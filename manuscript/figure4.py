# reproduction code for figure 4

import argparse
import os
import sys

_script_dir = os.path.dirname(os.path.abspath(__file__))
_sctapas_python = os.path.join(os.path.dirname(_script_dir), "python")
sys.path.insert(0, _sctapas_python)
sys.path.insert(0, _script_dir)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sctapas import plot_utils as pu
from sctapas import utils
import tcrqtl_eval as tqe

HLA_EVAL_CDICT = {
    "scTAPAS all":          "turquoise",
    "scTAPAS common":       "teal",
    "arcasHLA all":         "lightcoral",
    "arcasHLA common":      "brown",
    "Null (COMBAT) all":    "lightgrey",
    "Null (COMBAT) common": "darkgrey",
}

TCR_CDICT = {
    "CD4": (89/255, 170/255, 67/255, 1),
    "CD8": (158/255, 69/255, 137/255, 1),
}
SIG = 5e-8
LEV = "Annotation_major_subset"
CT_SEL = ["CD4", "CD8"]

def prep_hla_perf(perf_tab, sv_dir):
    os.makedirs(sv_dir, exist_ok=True)
    tab = perf_tab[np.logical_and(perf_tab["resolution"] == "4d", perf_tab["measure"] == "genewise_acc")].rename(columns={"gene": "Gene", "value": "Acc", "run": "neat_run"})
    tab["neat_run"] = tab["neat_run"] + " all" 
    tab_common = perf_tab[np.logical_and(perf_tab["resolution"] == "4d_common", perf_tab["measure"] == "genewise_acc")].rename(columns={"gene": "Gene", "value": "Acc", "run": "neat_run"})
    tab_common["neat_run"] = tab_common["neat_run"] + " common"
    return pd.concat([tab, tab_common])

@pu.name_output
def figure_4a(perf_tab, ext, sv_dir):
    tab = prep_hla_perf(perf_tab, sv_dir)
    cdict = HLA_EVAL_CDICT
    run_order = dict(zip(cdict.keys(), range(len(cdict))))
    gene_data = tab[np.isin(tab["Gene"], utils.classical_genes)]
    gene_order = dict(zip(utils.classical_genes, range(len(utils.classical_genes))))
    gene_data = gene_data[gene_data["Gene"].isin(utils.classical_genes)]
    gene_data = gene_data.sort_values(by=["neat_run", "Gene"], key=lambda y: [run_order.get(yi, 0) if yi in run_order else gene_order.get(yi, 0) for yi in y],)
    fig, ax = plt.subplots(figsize=(4.5, 3))
    pu.bars_either(fig, ax, gene_data, "Gene", "Acc", "neat_run", cdict)
    pu.AbsPlot(mode="save", sv_dir=sv_dir, ext=ext).handle("figure_4a")
    plt.close()

@pu.name_output
def figure_4b(summ_data, ext, sv_dir):
    cdict = HLA_EVAL_CDICT
    run_order = dict(zip(cdict.keys(), range(len(cdict))))
    summ_order = {"All": 0, "class 1": 1, "class 2": 2}
    summ_data = summ_data.sort_values(by=["neat_run", "Gene"], key=lambda y: [run_order.get(yi, 0) if yi in run_order else summ_order.get(yi, 0) for yi in y])
    fig, ax = plt.subplots(figsize=(2, 3))
    pu.bars_either(fig, ax, summ_data, "Gene", "Acc", "neat_run", cdict)
    pu.AbsPlot(mode="save", sv_dir=sv_dir, ext=ext).handle("figure_4b")
    plt.close()

def figure_4cd(cell_type, sctapas_df, ext, sv_dir):
    os.makedirs(sv_dir, exist_ok=True)
    tcr = tqe.TCRQTLPlots(ext=ext, mode="save", sv_dir=sv_dir)
    top_nlog10p = -np.log10(tqe.select_plink(plink_df=sctapas_df, name_sel=["CD4", "CD8"])["P"].min())
    tcr.meta_mhc_manhattan(tqe.select_plink(plink_df=sctapas_df, name_sel=cell_type), top=top_nlog10p, sig=-np.log10(SIG), figsize=(4, 3), plt_add=f"{cell_type}_")

@pu.name_output
def figure_4c(sctapas_df, ext, sv_dir):
    figure_4cd("CD4", sctapas_df=sctapas_df, ext=ext, sv_dir=sv_dir)

@pu.name_output
def figure_4d(sctapas_df, ext, sv_dir):
    figure_4cd("CD8", sctapas_df=sctapas_df, ext=ext, sv_dir=sv_dir)

@pu.name_output
def figure_4e(drb103_trav12_2_tab, ext, sv_dir):
    os.makedirs(sv_dir, exist_ok=True)
    tcr = tqe.TCRQTLPlots(ext=ext, mode="save", sv_dir=sv_dir)
    tcr.assoc_boxplot("HLA_DRB1*03", "TRAV12_2_clone", drb103_trav12_2_tab, "Annotation_major_subset", TCR_CDICT)

@pu.name_output
def figure_4f(dqb102_trav21_tab, ext, sv_dir):
    os.makedirs(sv_dir, exist_ok=True)
    tcr = tqe.TCRQTLPlots(ext=ext, mode="save", sv_dir=sv_dir)
    tcr.assoc_boxplot("HLA_DQB1*02", "TRAV21_clone", dqb102_trav21_tab, "Annotation_major_subset", TCR_CDICT)

@pu.name_output
def figure_4g(both_2d, ext, sv_dir):
    os.makedirs(sv_dir, exist_ok=True)
    tcr_plots = tqe.TCRQTLPlots(ext=ext, mode="save", sv_dir=sv_dir)
    sig = SIG
    cdict = TCR_CDICT
    fig, ax = plt.subplots(figsize=(6, 2))
    w, s = 0.2, 0
    n_id = len(both_2d["ID"].unique())
    for ct in both_2d["SubsetName"].unique():
        s += 0.1
        color = cdict[ct]
        res_ct = both_2d[both_2d["SubsetName"] == ct]
        edgec_sct = "k" if (res_ct["P_scTAPAS"] < sig).all() else None
        edgec_gsa = "k" if (res_ct["P_GSA"] < sig).all() else None
        ax.errorbar(x=res_ct["BETA_scTAPAS"], y=np.arange(n_id) - w - s, xerr=res_ct[["L95_err_scTAPAS", "U95_err_scTAPAS"]].values.T, fmt="o", c=color, capsize=3, markeredgecolor=edgec_sct)
        ax.errorbar(x=res_ct["BETA_GSA"],     y=np.arange(n_id) + w - s, xerr=res_ct[["L95_err_GSA",     "U95_err_GSA"]].values.T,     fmt="v", c=color, capsize=3, markeredgecolor=edgec_gsa)
    ax.set_xlabel(r"$\beta$")
    ax.set_yticks(np.arange(n_id) - 0.05 * len(both_2d["SubsetName"].unique()), both_2d["ID"].str.replace("_", "-").unique())
    ax.axvline(x=0, c="darkgrey", linestyle="--")
    cleg_pos = pu.find_right_legend_loc(fig, [ax], [ax], buffer=0.02)
    marker_handles = [
        mpl.lines.Line2D([0], [0], marker="o", linestyle="none", label="scTAPAS", c="k"),
        mpl.lines.Line2D([0], [0], marker="v", linestyle="none", label="array +\nHLA imputation", c="k"),
        mpl.lines.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="darkgray", markeredgecolor="k", label="significant"),
        mpl.lines.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="darkgray", markeredgecolor="darkgray", label="not significant"),
    ]
    fig.legend(handles=marker_handles, loc="center left", bbox_to_anchor=(cleg_pos[0], cleg_pos[1] - 0.1), frameon=False)
    pu.add_legend(fig, cdict, loc="center left", bbox=(cleg_pos[0], cleg_pos[1] + 0.3))
    tcr_plots.handle("scTAPAS_GSA_forest")


def parse_args():
    p = argparse.ArgumentParser(description="Reproduce Figure 4 (HLA eval) and Figure S4 (TCR QTL).")
    p.add_argument("--sv_dir", required=True, help="Output root directory")
    p.add_argument("--ext", default=".png", help="file extension for plots")
    p.add_argument("--hla_perf_table", help="Pre-computed performance table (Supplementary Table 7)")
    p.add_argument("--hla_summ_perf_table", help="Pre-computed performance summary table")
    p.add_argument("--sctapas_res_path", help="scTAPAS HLA-TCR association study summary statistics (Supplementary Table 8)")
    p.add_argument("--drb103_trav12_2_path", help="scTAPAS HLA-TCR input table for HLA-DRB1*03 and TRAV12-2")
    p.add_argument("--dqb102_trav21_path", help="scTAPAS HLA-TCR input table for HLA-DQB1*02 and TRAV21")
    p.add_argument("--scTAPAS_array_compare_path", help="scTAPAS and array HLA-TCR result comparison")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    perf_tab = pd.read_table(args.hla_perf_table, comment="#")
    figure_4a(perf_tab, ext=args.ext, sv_dir=args.sv_dir)
    summ_perf_tab = pd.read_table(args.hla_summ_perf_table)
    figure_4b(summ_perf_tab, ext=args.ext, sv_dir=args.sv_dir)
    sctapas_df = pd.read_table(args.sctapas_res_path, comment="#")
    figure_4c(sctapas_df=sctapas_df, ext=args.ext, sv_dir=args.sv_dir)
    figure_4d(sctapas_df=sctapas_df, ext=args.ext, sv_dir=args.sv_dir)
    drb103_tab = pd.read_table(args.drb103_trav12_2_path, comment="#")
    figure_4e(drb103_trav12_2_tab=drb103_tab, ext=args.ext, sv_dir=args.sv_dir)
    dqb102_tab = pd.read_table(args.dqb102_trav21_path, comment="#")
    figure_4f(dqb102_trav21_tab=dqb102_tab, ext=args.ext, sv_dir=args.sv_dir)
    comp_tab = pd.read_table(args.scTAPAS_array_compare_path)
    figure_4g(both_2d=comp_tab, ext=args.ext, sv_dir=args.sv_dir)
