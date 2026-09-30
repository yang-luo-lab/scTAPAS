import pandas as pd
import scipy as sp
import numpy as np

def remove_gs_allele(gs_col):
    return gs_col.str.split("*", expand=True)[0]


def select_by_col(tcr_tab, col_name, var_name):
    return tcr_tab[tcr_tab[col_name]==var_name]


def inverse_normal_rank(tcr_tab):
    nrank = (tcr_tab.rank(axis=0) - 0.5)/tcr_tab.count()
    return pd.DataFrame(data=sp.stats.norm.ppf(nrank), index=tcr_tab.index, columns=tcr_tab.columns)

def make_ann_tcr(tcr_path, gex_path, tcr_barcode_col_name="barcode_id", ab_col_name="chain_composition", ab_col_value="double_alpha_beta", sample_id_col_name="sample_id"):
    import scanpy as sc
    tcr = pd.read_table(tcr_path)
    gex = sc.read(gex_path)
    # merge with gex
    tcr_gex = tcr.join(gex.obs, on=tcr_barcode_col_name, rsuffix = "_gex", how="inner")
    # extract only ab
    abtcr = tcr_gex[tcr_gex[ab_col_name]==ab_col_value]
    abtcr = abtcr.dropna(subset=sample_id_col_name)
    return abtcr

def demux_qc_clones(tcr, sample_name_col="sample_id", pool_name="pool", trav_name="v_gene_TRA",
                     CDR3a_nt_name="cdr3_nt_TRA", traj_name="j_gene_TRA", trbv_name="v_gene_TRB",
                     CDR3b_nt_name="cdr3_nt_TRB", trbj_name="j_gene_TRB"):
    tcr["clone_id_nt"] = tcr[trav_name] + tcr[CDR3a_nt_name] + tcr[traj_name] + tcr[trbv_name] + tcr[CDR3b_nt_name] + tcr[trbj_name]
    flags = tcr.groupby("clone_id_nt").apply(flag_cell_sample_clone, sample_name_col=sample_name_col, pool_name=pool_name)
    tcr = tcr.join(flags.to_frame("clone_owner"), on="clone_id_nt")
    return tcr[tcr["clone_owner"].fillna(tcr[sample_name_col])==tcr[sample_name_col]]
    
def flag_cell_sample_clone(clone, sample_name_col="sample_id", pool_name="pool"):
    cpsc = cell_per_sample_clone(clone=clone, sample_name_col=sample_name_col)
    # clone originates from only 1 pool
    if len(cpsc)==1:
        return np.nan
    else:
        # clone originates from different samples but they are in different pools
        if len(clone[pool_name].unique())==len(cpsc):
            return np.nan
        else:
            # at least two samples contributing to clone are in same pool
            top_n_cells = cpsc.max()
            is_top_cells = cpsc == top_n_cells
            if sum(is_top_cells)>1:
                return "amb"
            else:
                return str(np.squeeze(list(cpsc[is_top_cells].index)))

def cell_per_sample_clone(clone, sample_name_col="sample_id"):
    sams, counts = np.unique(clone[sample_name_col], return_counts=True)
    return pd.Series(counts, sams)

def extract_VDJ_usage(tcr_tab, vdj, locus, sample_id_col, clone_id_col, tsubset=None, by_clone=True):
    if tsubset is not None:
        print(tsubset)
        tcr_tab = select_by_col(tcr_tab, tsubset[0], tsubset[1])
    # note we are taking the first instance of v gene when multiple modes
    gn = f"{vdj}_gene_{locus}"
    tcr_tab[gn] = remove_gs_allele(tcr_tab[gn])
    if by_clone:
        gs = tcr_tab.groupby([sample_id_col, clone_id_col])[gn].agg(lambda x: pd.Series.mode(x)[0]).to_frame()
        gs = gs.reset_index()
        gs[gn] = gs[gn] + "_clone"
    else:
        gs = tcr_tab.copy()
        gs[gn] = gs[gn] + "_cell"
    gs_counts = gs.groupby([sample_id_col, gn]).agg(n=(clone_id_col, "count"))
    wide_gs_counts = gs_counts.reset_index().pivot(columns=gn, index=sample_id_col, values="n").fillna(0)
    wide_gs_counts.index.name = None
    wide_gs_counts.columns.name = None
    return wide_gs_counts

def get_vdj_usg(tcr_tab, sams, gene, locus, sample_id_col, clone_id_col, tsub, min_cln=1, min_obs_prop=0.05, by_clone=True):
    print(tsub)
    gtab = extract_VDJ_usage(tcr_tab, gene, locus, sample_id_col, clone_id_col, tsubset=tsub, by_clone=by_clone)
    # remove samples not in sample list, if given
    gtab = gtab[np.isin(gtab.index, sams)]
    gtab = gtab.reindex(sams)
    # remove rows/samples with minimum number of clones/cells
    gtab = gtab[gtab.sum(axis=1)>=min_cln]
    prop_obs = (gtab>0).sum(axis=0)/(len(gtab) - gtab.isna().sum())
    gtab = gtab[gtab.columns[prop_obs >= min_obs_prop]]
    ts_neat = "_".join(tsub)
    ts_neat = ts_neat.replace("/", "-")
    if by_clone:
        n = tcr_tab.groupby([sample_id_col])[clone_id_col].nunique()
    else:
        n = tcr_tab.groupby([sample_id_col]).size()
    gtab_norm = (gtab.T/n.loc[gtab.index]).T
    gtab_irnt = inverse_normal_rank(gtab_norm)
    gtab_irnt = gtab_irnt.fillna("NA")
    vgene_names = gtab_irnt.columns
    gtab_irnt = gtab_irnt.reset_index(names=["IID"])
    gtab_irnt["FID"] = gtab_irnt["IID"]
    neat_vgene_names = [vg.replace("-", "_").replace("/", "_") for vg in vgene_names]
    gtab_irnt_sel = gtab_irnt[["FID", "IID", *vgene_names]]
    gtab_irnt_sel.columns = ["FID", "IID", *neat_vgene_names]
    return gtab_irnt_sel