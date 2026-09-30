# reproduction code for figure 3 and figure s2

import argparse
import os
import sys

_script_dir = os.path.dirname(os.path.abspath(__file__))
_sctapas_python = os.path.join(os.path.dirname(_script_dir), "python")
sys.path.insert(0, _sctapas_python)

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib_venn
import numpy as np
import pandas as pd
from sctapas import plot_utils as pu
from sctapas.ancestry import AncestryPCA

def neaten_name(name):
    replacements = {"gt": "array + TOPMed", "scrna": "scTAPAS"}
    neat = " ".join(replacements.get(p, p) for p in name.split("_"))
    if "pval nominal" in neat:
        if "neglog10" in neat:
            neat = neat.replace("pval nominal", r"$-\log_{10}P_{\mathrm{nominal}}$").replace("neglog10", "")
        else:
            neat = neat.replace("pval nominal", "$P_{nominal}$")
    return neat

def my_compute_colours(color_a, color_b, color_ab):
    ccv = mpl.colors.ColorConverter()
    return [np.array(ccv.to_rgb(c)) for c in [color_a, color_b, color_ab]]

def add_egene_cols(matches):
    # add helper columns the plots need
    matches = matches.copy()
    in_sct = matches["egene_in_scTAPAS"]
    in_gt = matches["egene_in_array_gt"]
    matches["egene_in_which"] = np.where(
        in_sct & in_gt, "both",
        np.where(in_sct, "scrna_only", "gt_only"))
    # both_lead: lead SNV identical in both datasets
    matches["both_lead"] = matches["scTAPAS_lead_SNV"] == matches["array_gt_lead_SNV"]
    return matches

class scRNAeQTLPlots(pu.AbsPlot):
    def __init__(self, mode="show", sv_dir=None, res=600, ext=".png"):
        super().__init__(mode=mode, sv_dir=sv_dir, res=res, ext=ext)

    def eqtl_scatter(self, xvals, xname, yvals, yname, cvals=None, cname=None, cdict=None, cmap=None, mvals=None, ms=4, mode_dict=None, xy=False, alpha=0.7, cvmax=1, cvmin=-1, cbar_label="sR2", mleg_colour="silver", figsize=(3, 3), plt_add=""):
        fig, ax = plt.subplots(figsize=figsize)
        if mode_dict is not None:
            mode = list(mode_dict.keys())[0]
            xvals = mode_dict[mode](xvals)
            yvals = mode_dict[mode](yvals)
            xname = f"{mode}_{xname}"
            yname = f"{mode}_{yname}"
        use_cdict = use_cbar = False
        if cdict is not None and cvals is not None:
            colours = [cdict[w] for w in np.atleast_1d(cvals.values)]
            use_cdict = True
        elif cvals is not None:
            if cmap is None:
                cmap = "viridis"
            norm = mpl.colors.Normalize(vmin=cvmin, vmax=cvmax)
            colours = [mpl.cm.get_cmap(cmap)(norm(c)) for c in cvals.values.reshape(-1)]
            use_cbar = True
        else:
            colours = np.repeat("k", len(xvals))
        xvals, yvals = np.asarray(xvals), np.asarray(yvals)
        mvals_arr = np.atleast_1d(mvals) if mvals is not None else None
        for i in range(len(xvals)):
            mval = mvals_arr[i] if mvals_arr is not None else None
            if mval == "lead_snp_pos_match":
                ax.scatter(xvals[i], yvals[i], facecolors="none", edgecolors=colours[i], alpha=alpha, marker="o", s=ms**2)
            elif mval in ("scrna_lead", "gt_lead", "single_lead"):
                ax.scatter(xvals[i], yvals[i], c=colours[i], alpha=alpha, marker="o", edgecolors=colours[i], s=ms**2)
            else:
                ax.scatter(xvals[i], yvals[i], c=colours[i], alpha=alpha, marker="o", edgecolors="k", s=ms**2)
        ax.set_xlabel(neaten_name(xname))
        ax.set_ylabel(neaten_name(yname))
        if xy:
            plt.gca().set_aspect('equal')
            ax = pu.add_xy(ax)
        single_legend_pos = pu.find_right_legend_loc(fig, [ax], [ax])
        has_pos_match = mvals_arr is not None and np.any(mvals_arr == "lead_snp_pos_match")
        has_single_lead = mvals_arr is not None and np.any(np.isin(mvals_arr, ["scrna_lead", "gt_lead", "single_lead"]))
        has_both_lead = mvals_arr is not None and np.any(mvals_arr == "both_lead")
        show_match_legend = has_pos_match or has_single_lead or has_both_lead
        if use_cdict:
            colour_legend_pos = (single_legend_pos[0], single_legend_pos[1] * (1.4 if show_match_legend else 1))
            pu.add_legend(fig, {neaten_name(k): v for k, v in cdict.items()}, loc="center left", bbox=colour_legend_pos, title="eGene in")
        if use_cbar:
            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.get_cmap(cmap)), ax=ax, location="bottom", shrink=0.6, label=cbar_label)
        if show_match_legend:
            match_handles = []
            if has_both_lead:
                match_handles.append(mpl.lines.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=mleg_colour, markeredgecolor="k", markersize=ms*1.5, label="both lead"))
            if has_single_lead:
                match_handles.append(mpl.lines.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=mleg_colour, markeredgecolor=mleg_colour, markersize=ms*1.5, label="single lead"))
            if has_pos_match:
                match_handles.append(mpl.lines.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="none", markeredgecolor=mleg_colour, markersize=ms*1.5, label="positional match"))
            fig.legend(handles=match_handles, loc="center left", bbox_to_anchor=(single_legend_pos[0], single_legend_pos[1] * 0.5), title="Match type", frameon=False)
        self.handle(f"{plt_add}{xname}_{yname}_{cname}_scatter", tight=True)


    def eqtl_venn(self, matches, cdict, a=0.8, include_cell_type=False, figsize=(4, 2), plt_add=""):
        matplotlib_venn._venn2._compute_colors = my_compute_colours
        matches = matches.copy()
        eqtl_colname = "phenotype_id_ct" if include_cell_type else "phenotype_id"
        if include_cell_type:
            matches["phenotype_id_ct"] = matches["phenotype_id"] + matches["cell_type"]
        in_sct = matches["egene_in_scTAPAS"]
        in_gt = matches["egene_in_array_gt"]
        eqtl_sets = (set(matches[in_sct][eqtl_colname]), set(matches[in_gt][eqtl_colname]))
        scrna = matches[in_sct][eqtl_colname].nunique()
        gt = matches[in_gt][eqtl_colname].nunique()
        scrna_only = matches[in_sct & ~in_gt][eqtl_colname].nunique()
        gt_only = matches[in_gt & ~in_sct][eqtl_colname].nunique()
        fig, ax = plt.subplots(figsize=figsize)
        matplotlib_venn._venn2.venn2(
            eqtl_sets,
            set_labels=(
                f"scTAPAS-only: {scrna_only/scrna:.1%}\nscTAPAS-shared: {1 - (scrna_only/scrna):.1%}",
                f"array+TOPMed-only: {gt_only/gt:.1%}\narray+TOPMed-shared: {1 - (gt_only/gt):.1%}\n",
            ),
            set_colors=(cdict["scrna_only"], cdict["gt_only"], cdict["both"]),
            alpha=a, ax=ax,
        )
        self.handle(f"{plt_add}eGene_venn_ct{include_cell_type}")
        plt.close()

    def nvar_bar(self, nvar, which_cdict):
        nvar_bar = nvar.sum()[["array+TOPMed", "scTAPAS"]]
        fig, ax = plt.subplots(figsize=(2, 3))
        ax.bar([i.replace("+", "\n+\n") for i in nvar_bar.index], nvar_bar.values, color=[which_cdict["gt_only"], which_cdict["scrna_only"]])
        ax.set_ylabel("number of variants")
        fmt = mpl.ticker.ScalarFormatter(useMathText=True)
        fmt.set_scientific(True)
        fmt.set_powerlimits((0, 0))
        fmt.set_useOffset(False)
        ax.yaxis.set_major_formatter(fmt)
        ax.set_axisbelow(True)
        self.handle("nvar_bar")

    def onek1k_bar(self, jtab, col_name, count_per_egene=True, figsize=(2, 4)):
        if count_per_egene:
            jtab = jtab.drop_duplicates("phenotype_id")
            ylab = "proportion egenes"
        else:
            ylab = "proportion cell type egenes"
        n_egenes = jtab.groupby(["egene_in_which", col_name]).size().to_frame(name="n_var").reset_index()
        total_egenes = jtab.groupby("egene_in_which").size().to_frame(name="total")
        n_egenes = n_egenes.join(total_egenes, on="egene_in_which")
        n_egenes["prop"] = n_egenes["n_var"] / n_egenes["total"]
        neg_wide = n_egenes.pivot(index="egene_in_which", columns=col_name, values="prop")
        fig, ax = plt.subplots(figsize=figsize)
        cdict = {True: "gold", False: "palevioletred"}
        neg_wide.index = [neaten_name(i) for i in neg_wide.index]
        ax.bar(neg_wide.index, neg_wide[True], color=cdict[True])
        ax.bar(neg_wide.index, neg_wide[False], bottom=neg_wide[True], color=cdict[False])
        ax.margins(y=0)
        plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
        for x, (true_val, false_val) in enumerate(zip(neg_wide[True], neg_wide[False])):
            if true_val > 0.03:
                ax.text(x, true_val / 2, f"{true_val*100:.3g}%", ha="center", va="center", fontsize=8)
            if false_val > 0.03:
                ax.text(x, true_val + false_val / 2, f"{false_val*100:.3g}%", ha="center", va="center", fontsize=8)
        pu.add_legend(fig, cdict={"in OneK1K": cdict[True], "not in OneK1K": cdict[False]}, bbox=pu.find_right_legend_loc(fig, [ax], [ax]), loc="center left")
        ax.set_ylabel(ylab)
        self.handle(f"{col_name}_per_egene{count_per_egene}_onek1k_prop_bar")

WHICH_CDICT = {"both": "mediumseagreen", "scrna_only": "dodgerblue", "gt_only": "tomato"}

@pu.name_output
def figure_3a(anc_pc_tab, ext, sv_dir):
    os.makedirs(sv_dir, exist_ok=True)
    sct_anc = AncestryPCA(ext=ext, mode="save", sv_dir=sv_dir)
    pcs = [f"PC{i+1}_AVG" for i in range(2)]
    target = anc_pc_tab[anc_pc_tab["dataset"]=="COMBAT"]
    ref = anc_pc_tab[anc_pc_tab["dataset"]=="1000G"]
    cdict = dict(zip(['EUR', 'SAS', 'EAS', 'AFR', 'AMR'], ['r', 'b', 'g', 'y', 'm']))
    sct_anc.scatter_target_ref_pca_split(target[pcs], ref[pcs], pcs, cdict, target["scTAPAS"], target["Genotyping"], ref["1000G"], ref_marker="o", target_marker="o", ref_alpha=0.1, target_edgecolor="k", figsize=(5, 4), plt_add="")

@pu.name_output
def figure_3b(nvar_tab, ext, sv_dir):
    plots = scRNAeQTLPlots(ext=ext, mode="save", sv_dir=sv_dir)
    plots.nvar_bar(nvar=nvar_tab, which_cdict=WHICH_CDICT)

@pu.name_output
def figure_3c(matches, ext, sv_dir):
    plots = scRNAeQTLPlots(ext=ext, mode="save", sv_dir=sv_dir)
    plots.eqtl_venn(matches=matches, cdict=WHICH_CDICT, include_cell_type=True)

@pu.name_output
def figure_3d(matches, ext, sv_dir):
    plots = scRNAeQTLPlots(ext=ext, mode="save", sv_dir=sv_dir)
    matches_exact = matches[matches["SNV_match_in_scTAPAS_array"]]
    cvals = matches_exact["egene_in_which"]
    mvals = np.where(matches_exact["both_lead"], "both_lead", "single_lead")
    plots.eqtl_scatter(matches_exact["pval_nominal_array_gt"].astype(float), "gt_pval_nominal", matches_exact["pval_nominal_scTAPAS"].astype(float), "scrna_pval_nominal", cvals, "egene_in_which", cdict=WHICH_CDICT, mvals=mvals, mode_dict={"neglog10": lambda y: -np.log10(y)}, xy=True, alpha=0.8)

@pu.name_output
def figure_3e(matches, ext, sv_dir):
    plots = scRNAeQTLPlots(ext=ext, mode="save", sv_dir=sv_dir)
    matches_exact = matches[matches["SNV_match_in_scTAPAS_array"]]
    cvals = matches_exact["egene_in_which"]
    mvals = np.where(matches_exact["both_lead"], "both_lead", "single_lead")
    plots.eqtl_scatter(matches_exact["slope_array_gt"], "gt_slope", matches_exact["slope_scTAPAS"], "scrna_slope", cvals, "egene_in_which", cdict=WHICH_CDICT, mvals=mvals, xy=True, alpha=0.8)

@pu.name_output
def figure_3f(matches, ext, sv_dir):
    plots = scRNAeQTLPlots(ext=ext, mode="save", sv_dir=sv_dir)
    plots.onek1k_bar(matches, "egene_in_1k1k", count_per_egene=True)

@pu.name_output
def figure_s2(matches, ext, sv_dir):
    plots = scRNAeQTLPlots(ext=ext, mode="save", sv_dir=sv_dir)
    matches_mismatch = matches.groupby("phenotype_id").apply(lambda y: True if y["egene_in_1k1k"].all() and not y["egene_in_1k1k_same_cell"].any() else False)
    matches = matches.join(matches_mismatch.to_frame(name="egene_in_1k1k_all_ct_mismatch"), on="phenotype_id")
    matches["egene_in_1k1k_any_same_cell"] = np.logical_and(matches["egene_in_1k1k"], np.logical_not(matches["egene_in_1k1k_all_ct_mismatch"]))
    plots.onek1k_bar(matches, "egene_in_1k1k_any_same_cell", count_per_egene=True)

def parse_args():
    p = argparse.ArgumentParser(description="Reproduce Figure 3 (ancestry and eQTL) and Figure S2.")
    p.add_argument("--ext", default=".png", help="file extension for plots")
    p.add_argument("--sv_dir", help="output directory for plots")
    p.add_argument("--matches_path", help="Pre-computed scrna_gt_eqtl_matches.tsv (with OneK1K columns)")
    p.add_argument("--nvar_path", help="Table with number of variants in each eQTL analysis, supp table 6")
    p.add_argument("--anc_pca_path", help="Pre-computed table for ancestry")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    anc_pca_tab = pd.read_table(args.anc_pca_path)
    figure_3a(anc_pc_tab=anc_pca_tab, ext=args.ext, sv_dir=args.sv_dir)
    nvar_tab = pd.read_table(args.nvar_path, comment="#")
    figure_3b(nvar_tab=nvar_tab, ext=args.ext, sv_dir=args.sv_dir)
    matches_all = add_egene_cols(pd.read_table(args.matches_path, comment="#"))
    figure_3c(matches=matches_all, ext=args.ext, sv_dir=args.sv_dir)
    figure_3d(matches=matches_all, ext=args.ext, sv_dir=args.sv_dir)
    figure_3e(matches=matches_all, ext=args.ext, sv_dir=args.sv_dir)
    figure_3f(matches=matches_all, ext=args.ext, sv_dir=args.sv_dir)
    figure_s2(matches=matches_all, ext=args.ext, sv_dir=args.sv_dir)

