import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from adjustText import adjust_text
from sctapas import plot_utils as pu
import glob
from sctapas import utils
import scipy as sp

class TCRQTLPlots(pu.AbsPlot):
    def __init__(self, mode="show", sv_dir=None, res=600, ext=".png"):
        super().__init__(mode, sv_dir, res, ext)
        
    def mhc_manhattan(self, gwas, ann_ids=None, ann_pos=None, ann_id_func=None, adjust=False, buffer=500, top=None, sig=None, figsize=(5,4), plt_add=""):
        fig, ax = plt.subplots(figsize=figsize)
        pu.mhc_manhattan(fig, ax, gwas, ann_ids=ann_ids, ann_pos=ann_pos, ann_id_func=ann_id_func, adjust=adjust, buffer=buffer)
        if sig is not None:
            ax.axhline(y=sig, c="k", linestyle="--", alpha=0.2, zorder=0)
        ax.set_ylim(bottom=0, top=top)
        self.handle(f"{plt_add}mhc_manhattan")

    def meta_mhc_manhattan(self, plink_df=None, fn_patt=None, top=None, ann_pos=None, ann_id_func=None, annot_top=False, sig=None, figsize=(5,4), plt_add=""):
        if fn_patt is not None:
            all_fn = glob.glob(fn_patt)
            all_f = []
            for fn in all_fn:
                po = pd.read_table(fn)
                po["ID"] = fn.split(".")[-3].replace("_", "-") + "_" + po["ID"].astype(str)
                all_f.append(po)
            meta = pd.concat(all_f, axis=0)
        elif plink_df is not None:
            meta = plink_df
        else:
            ValueError("must supply filename pattern or plink dataframe")
        if annot_top:
            ann_ids = [meta[meta["P"].values == meta["P"].min()]["ID"].values[0]]
        else:
            ann_ids = None
        self.mhc_manhattan(meta, ann_ids=ann_ids, ann_id_func=ann_id_func, top=top, ann_pos=ann_pos, sig=sig, figsize=figsize, plt_add=plt_add)

    def assoc_boxplot(self, vid, seg, data_long, lev, cdict, figsize=(4, 3), plt_add=""):
        fig, ax = plt.subplots(figsize=figsize)
        usg_name = f"{seg} usage"
        data_long[vid] = data_long[vid].astype(int)
        sns.swarmplot(data=data_long, x=vid, y=usg_name, hue=lev, ax=ax, palette=cdict, legend=False, dodge=True)
        sns.boxplot(data=data_long, x=vid, y=usg_name, hue=lev, ax=ax, palette=cdict, legend=False, fill=False, gap=0.1, fliersize=0)
        pu.add_legend(fig, {k: v for k, v in cdict.items() if k in data_long[lev].unique()}, bbox=pu.find_right_legend_loc(fig, [ax], [ax]), loc="center left")
        ax.set_xlabel(vid.replace("_" , " ") + " dosage")
        ax.set_ylabel(usg_name.replace("_" , " "))
        self.handle(f"{plt_add}{vid}_{seg}_{lev}_boxplot")

def filename_split(fn):
    fn = fn.replace("TEM_TEMRA", "TEMTEMRA")
    stop_split = fn.split(".")
    pheno_name = stop_split[-3]
    pheno_origin = ".".join(stop_split[:-3])
    und_split = pheno_origin.split("_")
    dataset, gene, locus, subset_level1, subset_level2, subset_level3, subset_name, usage_mode, pheno, hla, multi, assoc = und_split
    return gene, locus, "_".join([subset_level1, subset_level2, subset_level3]), subset_name, pheno_name

def list_plink_model_outs(po_dir=None, po_filenames=None, model_ext=".linear"):
    if po_dir is not None:
        plink_outs = os.listdir(po_dir)
    else:
        plink_outs = po_filenames
    return np.array(plink_outs)[[model_ext in po for po in plink_outs]]

def get_plink(po_dir=None, po_filenames=None, model_ext=".linear", md_query=None, var_query=None, fn_split_func=filename_split):
    plink_results = list_plink_model_outs(po_dir=po_dir, po_filenames=po_filenames, model_ext=model_ext)#np.array(plink_outs)[[".linear" in po for po in plink_outs]]
    if md_query is None:
        md_query = {}
    if var_query is None:
        var_query = {}
    if po_dir is None:
        po_dir = ""
    all_pr = []
    for pr in plink_results:
        md_dict = dict(zip(["Gene", "Locus", "SubsetLevel", "SubsetName", "Seg"], fn_split_func(pr)))
        md_dict["FileName"] = pr
        open_plink = True
        for qcol in md_query.keys():
            if md_dict[qcol] not in np.array(md_query[qcol]).flatten():
                open_plink = False
            if not open_plink:
                break
        if open_plink:
            plink_res = pd.read_csv(os.path.join(po_dir, pr), sep="\t")
            # variant query separate because requires opening file
            for var_field, var_func in var_query.items():
                var_match = var_func(plink_res[var_field])
                plink_res = plink_res[var_match]
            for mdcol in md_dict.keys():
                plink_res[mdcol] = md_dict[mdcol]
            all_pr.append(plink_res)
    if len(all_pr) == 0:
        all_pr_ord = None
    else:
        all_pr_df  = pd.concat(all_pr)
        all_pr_ord = all_pr_df.sort_values(by="P", ascending=True).reset_index(drop=True)
    return all_pr_ord

def select_plink(po_dir=None, plink_df=None, level_sel=None, name_sel=None, var_patt=None):
    if po_dir is not None:
        plink_df = get_plink(po_dir=po_dir)
    if var_patt is not None:
        plink_df = plink_df[plink_df["ID"].str.contains(var_patt)]
    if level_sel is not None:
        plink_df = plink_df[np.isin(plink_df["SubsetLevel"], level_sel)]
    if name_sel is not None:
        plink_df = plink_df[np.isin(plink_df["SubsetName"], name_sel)]
    return plink_df

