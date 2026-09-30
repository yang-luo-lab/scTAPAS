# reproduction code for figure 2 and figure s1

import argparse
import os
import sys

# put python package and this dir on path so imports resolve
_script_dir = os.path.dirname(os.path.abspath(__file__))
_sctapas_python = os.path.join(os.path.dirname(_script_dir), "python")
sys.path.insert(0, _sctapas_python)
sys.path.insert(0, _script_dir)

import pandas as pd
import quilt2_eval as q2e
from sctapas import plot_utils as pu

def make_r2_2d_hist(eval_name, r2_tab, ext, sv_dir):
    os.makedirs(sv_dir, exist_ok=True)
    plt_add = f"{eval_name}_"
    q2p = q2e.QUILT2_plots(ext=ext, mode="save", sv_dir=sv_dir)
    r2_tab["R2 with genotyping"] = r2_tab["R2_q2DS-gt"]
    q2p.panels_2d_1dhist(r2_tab, "INFO_SCORE", "R2 with genotyping", 50,
                                      cmap="Blues", xmin=0, ymin=0, plt_add=plt_add)


def make_r2_maf_boxplot(eval_name, r2_tab, maf_name, ext, sv_dir):
    os.makedirs(sv_dir, exist_ok=True)
    plt_add = f"{eval_name}_"
    q2p = q2e.QUILT2_plots(ext=ext, mode="save", sv_dir=sv_dir)
    r2_tab["R2 with genotyping"] = r2_tab["R2_q2DS-gt"]
    q2p.maf_v_r2_boxplot_info_cond(
                    r2_tab, "R2 with genotyping", "INFO_SCORE", maf_name, "R2 with genotyping",
                    maf_name.replace("_", " from "), info_threshs=[0, 0.5, 0.8], plt_add=plt_add)


def make_r2_info_boxplot(eval_name, r2_tab, ext, sv_dir):
    os.makedirs(sv_dir, exist_ok=True)
    plt_add = f"{eval_name}_"
    q2p = q2e.QUILT2_plots(ext=ext, mode="save", sv_dir=sv_dir)
    r2_tab["R2 with genotyping"] = r2_tab["R2_q2DS-gt"]
    q2p.r2_boxplot_info_cond(r2_tab, "R2 with genotyping", "INFO_SCORE", "R2 with genotyping",
                info_threshs=[0, 0.5, 0.8], plt_add=plt_add)

def make_r2_hwe_boxplot(eval_name, r2_tab, ext, sv_dir):
    os.makedirs(sv_dir, exist_ok=True)
    plt_add = f"{eval_name}_"
    q2p = q2e.QUILT2_plots(ext=ext, mode="save", sv_dir=sv_dir)
    r2_tab["R2 with genotyping"] = r2_tab["R2_q2DS-gt"]
    q2p.r2_boxplot_hwe_cond(r2_tab, "R2 with genotyping", "HWE", "R2 with genotyping",
    hwe_thresh=1e-6, plt_add=plt_add)

def var_bar(info_summ, sv_dir, ext=".png", prop=True, plt_add=""):
    os.makedirs(sv_dir, exist_ok=True)
    q2p = q2e.QUILT2_plots(ext=ext, mode="save", sv_dir=sv_dir)
    if prop:
        q2p.var_prop_bar(info_summ, cond="pass_qc", cond_false_name="fail_qc", plt_add=plt_add)
    else:
        q2p.var_bar_w_cutoff(info_summ, cond="pass_qc", cond_false_name="fail_qc", plt_add=plt_add)

@pu.name_output
def fig_2a(r2_tab_gsa, ext, sv_dir):
    make_r2_2d_hist("GSA", r2_tab_gsa, ext, sv_dir)

@pu.name_output
def fig_2b(r2_tab_topmed, ext, sv_dir):
    make_r2_2d_hist("topmed", r2_tab_topmed, ext, sv_dir)

@pu.name_output
def fig_2c(r2_tab_gsa, ext, sv_dir):
    make_r2_info_boxplot("GSA", r2_tab_gsa, ext, sv_dir)
    
@pu.name_output
def fig_2d(r2_tab_topmed, ext, sv_dir):
    make_r2_info_boxplot("topmed", r2_tab_topmed, ext, sv_dir)

@pu.name_output
def fig_2e(r2_tab_gsa, ext, sv_dir):
    make_r2_maf_boxplot("GSA", r2_tab_gsa, "MAF_gt", ext, sv_dir)

@pu.name_output
def fig_2f(r2_tab_topmed, ext, sv_dir):
    make_r2_maf_boxplot("topmed", r2_tab_topmed, "MAF_gt", ext, sv_dir)

@pu.name_output
def fig_2g(info_summ_var_annot, sv_dir, ext):
    var_bar(info_summ_var_annot, sv_dir, ext=ext, prop=False, plt_add="var_annot_")

@pu.name_output
def fig_2h(info_summ_var_annot, sv_dir, ext):
    var_bar(info_summ_var_annot, sv_dir, ext=ext, prop=True, plt_add="var_annot_")


@pu.name_output
def fig_s1a(r2_tab_gsa, ext, sv_dir):
    make_r2_maf_boxplot("GSA", r2_tab_gsa, "MAF_1000G", ext, sv_dir)

@pu.name_output
def fig_s1b(r2_tab_topmed, ext, sv_dir):
    make_r2_maf_boxplot("topmed", r2_tab_topmed, "MAF_1000G", ext, sv_dir)

@pu.name_output
def fig_s1c(r2_tab_gsa, ext, sv_dir):
    make_r2_hwe_boxplot("GSA", r2_tab_gsa, ext, sv_dir)

@pu.name_output
def fig_s1d(r2_tab_topmed, ext, sv_dir):
    make_r2_hwe_boxplot("topmed", r2_tab_topmed, ext, sv_dir)

@pu.name_output
def fig_s1e(r2_sam_tab_path, sv_dir=None, ext=".png"):
    # per-sample r2 (quilt2 DS vs genotyping GT) against cell count
    r2_sam = pd.read_table(r2_sam_tab_path)
    q2p = q2e.QUILT2_plots(ext=ext, mode="save", sv_dir=sv_dir)
    os.makedirs(sv_dir, exist_ok=True)
    q2p.r2_sam_tab_plot(r2_sam, "n_cells", "R2_sam")

@pu.name_output
def fig_s1f(info_summ_gene_status, sv_dir, ext):
    var_bar(info_summ_gene_status, sv_dir, ext=ext, prop=False, plt_add="gene_status_")

@pu.name_output
def fig_s1g(info_summ_gene_status, sv_dir, ext):
    var_bar(info_summ_gene_status, sv_dir, ext=ext, prop=True, plt_add="gene_status_")


def parse_args():
    p = argparse.ArgumentParser(
        description="Reproduce autosomal QUILT2 evaluation plots for the manuscript."
    )
    p.add_argument("--r2_gsa", required=True,
                   help="Merged GSA R2 table TSV")
    p.add_argument("--r2_topmed", required=True,
                   help="Merged TopMed R2 table TSV")
    p.add_argument("--info_summ_var_annot", required=True,
                   help="Variant info summarised over variant annotations")
    p.add_argument("--info_summ_gene_status", required=True,
                   help="Variant info summarised over variant gene status")
    p.add_argument("--ext", default=".png", help="file extension for plots")
    p.add_argument("--sv_dir", help="Output directory for plots")
    # sample-wise R2 vs. cell count plot
    p.add_argument("--r2_sam_tab_path", required=True, help="table for sample-wise R2 plot with n_cells, R2")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    r2_tab_gsa = pd.read_table(args.r2_gsa)
    r2_tab_topmed = pd.read_table(args.r2_topmed)
    fig_2a(r2_tab_gsa=r2_tab_gsa, ext=args.ext, sv_dir=args.sv_dir)
    fig_2b(r2_tab_topmed=r2_tab_topmed, ext=args.ext, sv_dir=args.sv_dir)
    fig_2c(r2_tab_gsa=r2_tab_gsa, ext=args.ext, sv_dir=args.sv_dir)
    fig_2d(r2_tab_topmed=r2_tab_topmed, ext=args.ext, sv_dir=args.sv_dir)
    fig_2e(r2_tab_gsa=r2_tab_gsa, ext=args.ext, sv_dir=args.sv_dir)
    fig_2f(r2_tab_topmed=r2_tab_topmed, ext=args.ext, sv_dir=args.sv_dir)
    info_summ_var_annot = pd.read_table(args.info_summ_var_annot, index_col=0)
    fig_2g(info_summ_var_annot=info_summ_var_annot, ext=args.ext, sv_dir=args.sv_dir)
    fig_2h(info_summ_var_annot=info_summ_var_annot, ext=args.ext, sv_dir=args.sv_dir)
    fig_s1a(r2_tab_gsa=r2_tab_gsa, ext=args.ext, sv_dir=args.sv_dir)
    fig_s1b(r2_tab_topmed=r2_tab_topmed, ext=args.ext, sv_dir=args.sv_dir)
    fig_s1c(r2_tab_gsa=r2_tab_gsa, ext=args.ext, sv_dir=args.sv_dir)
    fig_s1d(r2_tab_topmed=r2_tab_topmed, ext=args.ext, sv_dir=args.sv_dir)
    fig_s1e(args.r2_sam_tab_path, sv_dir=args.sv_dir, ext=args.ext,)
    info_summ_gene_status = pd.read_table(args.info_summ_gene_status, index_col=0)
    fig_s1f(info_summ_gene_status=info_summ_gene_status, ext=args.ext, sv_dir=args.sv_dir)
    fig_s1g(info_summ_gene_status=info_summ_gene_status, ext=args.ext, sv_dir=args.sv_dir)
