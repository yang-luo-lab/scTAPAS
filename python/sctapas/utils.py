import gzip
import os
import pandas as pd
import numpy as np

# vcf loading

def extract_info_from_col(info_f):
    all_info_f = info_f.split(";")
    val_dict = dict([val.split("=") for val in all_info_f if "=" in val])
    [val_dict.update({stat: True}) for stat in all_info_f if "=" not in stat]
    return pd.Series(val_dict)

def get_chunk_resolver(out_dir):
    q2_vcfs = [vcf  for vcf in os.listdir(out_dir) if ".vcf" in vcf]
    q2_chunks = pd.DataFrame(data=np.array([vcf.split(".") for vcf in q2_vcfs])[:,1:3].astype(int), index=q2_vcfs, columns=["start", "end"]).sort_values("start")
    b_ovlp = q2_chunks.values.flatten()[1:-1].reshape(-1, 2)
    b_halfways = (b_ovlp[:,0] - b_ovlp[:,1])/2
    q2_chunks["start_sub"] = np.concatenate(([0], b_halfways))
    q2_chunks["end_sub"] = np.concatenate((b_halfways, [0]))
    q2_chunks["start_bound"] = q2_chunks["start"] + q2_chunks["start_sub"]
    q2_chunks["end_bound"] = q2_chunks["end"] - q2_chunks["end_sub"]
    return q2_chunks

def read_vcf(fp, field, data_extract_func, info_inds=[0, 1,2,3,4,7], field_ind=8, data_start_ind=9, data_delim=":", info_only_mode=False, expand_info=True):
    ext = fp.split(".")[-1]
    if ext == "vcf":
        f = open(fp)
    else:
        f = gzip.open(fp, "rt")
    lines = f.readlines()
    var_mat = []
    f.close()
    header = None
    for l in lines:
        if l.startswith("##"):
            continue
        elif l.startswith("#"):
            header = l
        else:
            # line is data
            var = l.replace("\n", "").split("\t")
            var_info = [var[i] for i in info_inds]
            if not info_only_mode:
                fields = np.array(var[field_ind].split(data_delim))
                field_j = np.squeeze(np.argwhere(fields == field))
                var_data = [data_extract_func(v, field_j, data_delim) for v in var[data_start_ind:]]
                var_mat.append([*var_info, *var_data])
            else:
                var_mat.append(var_info)
    var_df = pd.DataFrame(var_mat)
    header = header.replace("\n", "")
    hd = header.split("\t")
    if info_only_mode:
        var_df.columns = [hd[i] for i in info_inds]
    else:
        var_df.columns = [*[hd[i] for i in info_inds], *hd[data_start_ind:]]
    # now separate info
    if expand_info:
        info = var_df["INFO"].apply(extract_info_from_col)
        var_df = var_df.join(info)
    return var_df

def read_chunked_vcf(out_dir, field, field_func, var_inds=[0, 1, 2, 3, 4, 7], field_info_ind=8, sample_ind_start=9, expand_info=True):
    chunks = get_chunk_resolver(out_dir)
    all_vcf = []
    for vcf_name in chunks.index:
        vcf = read_vcf(fp=os.path.join(out_dir, vcf_name), field=field, data_extract_func=field_func, info_inds=var_inds, field_ind=field_info_ind, data_start_ind=sample_ind_start, expand_info=expand_info)
        chunk_pos = vcf["POS"].astype(int)
        sb = chunks["start_bound"].loc[vcf_name]
        eb = chunks["end_bound"].loc[vcf_name]
        chunk_keep = np.logical_and(chunk_pos>sb, chunk_pos<=eb)
        all_vcf.append(vcf[chunk_keep])
    return pd.concat(all_vcf)

def extr_vcf_float(v, j, delim):
    parts = v.split(delim)
    if j >= len(parts):
        return np.nan
    value = parts[j]
    return float(value) if value != "." else np.nan

def extr_vcf_gt(v, j, delim, phased=False):
    parts = v.split(delim)
    if j >= len(parts):
        return np.nan
    gt = parts[j]
    if phased:
        gt_delim = "|"
    else:
        gt_delim="/"
    return sum([float(a) if a.isdigit() else np.nan for a in gt.split(gt_delim)])

def extr_phased(v, j, delim):
    return extr_vcf_gt(v, j, delim, phased=True)

def get_vcf_like_dose(field, vcf_path=None, vcf_dir=None, phased=False, var_inds=[0, 1, 2, 3, 4, 7], field_info_ind=8, sample_ind_start=9, expand_info=True):
    if field == "DS":
        field_func = extr_vcf_float
    elif field == "GT":
        if phased:
            field_func = extr_phased
        else:
            field_func = extr_vcf_gt
    if vcf_dir is None:# fp, field, data_extract_func, info_inds=[0, 1,2,3,4,7], field_ind=8, data_start_ind=9, data_delim=":", info_only_mode=False, expand_info=True
        return read_vcf(fp=vcf_path, field=field, data_extract_func=field_func, info_inds=var_inds, field_ind=field_info_ind, data_start_ind=sample_ind_start, expand_info=expand_info)
    elif vcf_path is None:# out_dir, field, field_func,
        return read_chunked_vcf(vcf_dir, field, field_func, var_inds, field_info_ind, sample_ind_start, expand_info=expand_info)
    else:
        print("must specify path to vcf or directory containing .vcf or .vcf.gz")
        return None

def split_vcf_like_dose(field, vcf_path=None, vcf_dir=None, phased=False, var_inds=[0, 1,2,3,4,7], field_info_ind=8, sample_ind_start=9, col_func=None, col_func_kwargs=None, var_id_func=None):
    vcf = get_vcf_like_dose(field, vcf_path=vcf_path, vcf_dir=vcf_dir, phased=phased, var_inds=var_inds, field_info_ind=field_info_ind, sample_ind_start=sample_ind_start, expand_info=True)
    if 'chr' not in vcf["#CHROM"].values[0]:
        add_chr = True
    else:
        add_chr = False
    if var_id_func is None:
        var_id_func = get_myvarid
    vcf["my_index"] = var_id_func(vcf, add_chr=add_chr)
    vcf["na"] = vcf.isna().sum(axis=1)
    vcf["id_len"] = vcf["ID"].str.len()
    ddup_vcf = vcf.groupby("my_index").apply(lambda y: y.sort_values(by=["na", "id_len"]).iloc[0])
    return split_sample_info(ddup_vcf.drop(columns=["na", "id_len", "my_index"]), n_info=len(var_inds), col_func=col_func, col_func_kwargs=col_func_kwargs)

def get_myvarid(vcf_df, add_chr=False, headings=["#CHROM", "POS", "REF", "ALT"]):
    vid = vcf_df[headings[0]].astype(str) +":"+ vcf_df[headings[1]].astype(str)+"_"+ vcf_df[headings[2]]+"_"+vcf_df[headings[3]]
    if add_chr:
        return "chr" + vid
    else:
        return vid

def split_sample_info(vcf_df, n_info=6, col_func=None, col_func_kwargs=None):
    info_fields = []
    for i in vcf_df["INFO"]:
        info_fields.extend([ii.split("=")[0] for ii in i.split(";")])
    info_info_n = len(np.unique(info_fields))#np.amax([i.count(";") for i in vcf_df["INFO"]]) + 1
    vcf_mat = vcf_df[vcf_df.columns[n_info:-info_info_n]]
    if col_func is not None:
        if col_func_kwargs is None:
            col_func_kwargs = {}
        vcf_mat.columns = [col_func(c, **col_func_kwargs) for c in vcf_mat.columns]
    # assumes standard ordering of columns!
    info_df = vcf_df[np.concatenate((vcf_df.columns[:n_info], vcf_df.columns[-info_info_n:])).reshape(-1)]
    return vcf_mat, info_df

# hla gene metadata

# palettes
classical_genes = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
hla_def_cdict = dict(zip(classical_genes, ["tab:red", "tab:blue", "tab:green", "tab:pink", "y", "tab:brown", "tab:purple", "tab:orange"]))

# positions based on hg38 ensembl
classical_gene_pos = {
    "A": [29941260, 29949606], # ENSG00000206503.17
    "B": [31353195, 31367067], # ENSG00000234745.16
    "C": [31268746, 31272171], # ENSG00000204525.20
    "DRB1": [32577902, 32589848], # ENSG00000196126.14
    "DQA1": [32628179, 32659533], # ENSG00000196735.15
    "DQB1": [32640622, 32668383], # ENSG00000179344.19
    "DPA1": [33064569, 33080781], # ENSG00000231389.10
    "DPB1": [33073786, 33091655] # ENSG00000223865.14
    }
