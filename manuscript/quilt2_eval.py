import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
import matplotlib as mpl
import upsetplot
from scipy.stats import chi2_contingency
import statsmodels.api as sm

from sctapas import plot_utils as pu

maf_bins = [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]


class QUILT2_plots(pu.AbsPlot):
    def __init__(self, mode="show", sv_dir=None, res=600, ext=".png"):
        super().__init__(mode, sv_dir, res, ext)

    def panels_2d_1dhist(self, data, xname, yname, bins, cmap="Blues", xmin=-1, xmax=1, ymin=-1, ymax=1, hist_colour="k", figsize=(3,3), plt_add=""):
        data_nmiss = data[[xname, yname]].dropna()
        label_renamer = {xname: xname.replace("R2", "$R^2$"), yname: yname.replace("R2", "$R^2$")}
        pu.panels_2d_1dhist(data_nmiss[xname], xname, data_nmiss[yname], yname, bins, mode="hist", cmap=cmap, cname="number of variants", xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, hist1d_colour=hist_colour, label_renamer=label_renamer, figsize=figsize)
        res = sp.stats.pearsonr(data_nmiss[xname], data_nmiss[yname])
        self.handle(f"{plt_add}{xname}_{yname}_panel_hist", tight=True)
        return pd.Series([res.statistic, res.pvalue], index=["r", "p"])

    def maf_v_r2_boxplot_info_cond(self, data, r2_name, info_name, maf_name, neat_r2_name="R2", neat_maf_name="MAF",
                                   info_threshs=[0, 0.5, 0.8, 0.9], figsize=(5, 3), ylabel_offset=None, plt_add=""):
        got_nans = data[[r2_name, info_name, maf_name]].isna().any(axis=1)
        print(f"removing {sum(got_nans)} rows with missing data")
        neat_r2_name = neat_r2_name.replace("R2", "$R^2$")
        conds = [data[info_name]>t for t in info_threshs]
        digitized = [np.digitize(data[maf_name][c], maf_bins) for c in conds]
        ys = [data[r2_name][c] for c in conds]
        range_cons = range(len(conds))
        cond_names = [rf"{info_name}$ > {t}$" if t>0 else "all" for t in info_threshs]
        maf_bins_neat = [f"{int(100*maf_bins[i])}%-{int(100*maf_bins[i+1])}%" for i in range(len(maf_bins)-1)]
        used_bins = [j for j in range(len(maf_bins_neat))
                     if len(pd.concat([ys[i][digitized[i]==(j+1)] for i in range_cons])) > 0]
        fig, ax = plt.subplots(1, len(used_bins), figsize=figsize, sharey=True)
        if len(used_bins) == 1:
            ax = [ax]
        colours = list(mpl.colormaps['Dark2'].colors)
        for ax_idx, j in enumerate(used_bins):
            mb = maf_bins_neat[j]
            y_con = [ys[i][digitized[i]==(j+1)] for i in range_cons]
            box = ax[ax_idx].boxplot(y_con, positions=range_cons, showfliers=False)
            for item in ['boxes', 'medians']:
                [plt.setp(box[item][i], color=colours[i]) for i in range_cons]
            for item in ['whiskers', 'fliers', 'caps']:
                start = 0
                for i in np.array(range_cons):
                    end = start + 2
                    plt.setp(box[item][start:end], color=colours[i])
                    start = end
            ax[ax_idx].set_ylim(-0.1, 1.1)
            ax[ax_idx].set_xticks([np.mean(range_cons)], [mb])
        ax[round(len(used_bins)/2)].set_xlabel(neat_maf_name)
        fig.supylabel(neat_r2_name, **({} if ylabel_offset is None else {"x": ylabel_offset}))
        pu.add_legend(fig, dict(zip(cond_names, colours)), loc="center left", bbox=pu.find_right_legend_loc(fig, ax, ax))
        self.handle(f"{plt_add}{r2_name}_{maf_name}_info_box")

    def r2_boxplot_info_cond(self, data, r2_name, info_name, neat_r2_name="R2", info_threshs=[0, 0.5, 0.8, 0.9], figsize=(3, 2), plt_add=""):
        got_nans = data[[r2_name, info_name]].isna().any(axis=1)
        print(f"removing {sum(got_nans)} rows with missing data")
        neat_r2_name = neat_r2_name.replace("R2", "$R^2$")
        conds = [data[info_name]>t for t in info_threshs]
        ys = [data[r2_name][c] for c in conds]
        range_cons = range(len(conds))
        cond_names = [rf"{info_name}$ > {t}$" if t>0 else "all" for t in info_threshs]
        fig, ax = plt.subplots(1, figsize=figsize, sharey=True)
        colours = list(mpl.colormaps['Dark2'].colors)
        y_con = [ys[i] for i in range_cons]
        box = ax.boxplot(y_con, positions=range_cons, showfliers=False)
        for item in ['boxes', 'medians']:
            [plt.setp(box[item][i], color=colours[i]) for i in range_cons]
        for item in ['whiskers', 'fliers', 'caps']:
            start = 0
            for i in np.array(range_cons):
                end = start + 2
                plt.setp(box[item][start:end], color=colours[i])
                start = end
            ax.set_ylim(-0.1, 1.1)
        ax.set_xticks(range_cons, cond_names)
        plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
        ax.set_ylabel(neat_r2_name)
        self.handle(f"{plt_add}{r2_name}_info_box")

    def r2_boxplot_hwe_cond(self, data, r2_name, hwe_name, neat_r2_name="R2", hwe_thresh=1e-6, figsize=(3, 3), ylabel_offset=None, plt_add=""):
        got_nans = data[[r2_name, hwe_name]].isna().any(axis=1)
        print(f"removing {sum(got_nans)} rows with missing data")
        neat_r2_name = neat_r2_name.replace("R2", "$R^2$")
        conds = [data[hwe_name]>=hwe_thresh, data[hwe_name]<hwe_thresh]
        ys = [data[r2_name][c] for c in conds]
        range_cons = range(len(conds))
        cond_names = [rf"{hwe_name} $\geq {hwe_thresh}$", f"{hwe_name} $< {hwe_thresh}$"]
        fig, ax = plt.subplots(1, figsize=figsize, sharey=True)
        colours = ["black", "firebrick"]
        y_con = [ys[i] for i in range_cons]
        box = ax.boxplot(y_con, positions=range_cons, showfliers=False)
        for item in ['boxes', 'medians']:
            [plt.setp(box[item][i], color=colours[i]) for i in range_cons]
        for item in ['whiskers', 'fliers', 'caps']:
            start = 0
            for i in np.array(range_cons):
                end = start + 2
                plt.setp(box[item][start:end], color=colours[i])
                start = end
            ax.set_ylim(-0.1, 1.1)
        ax.set_xticks(range_cons, cond_names)
        ax.set_ylabel(neat_r2_name)
        self.handle(f"{plt_add}{r2_name}_hwe_box")


    def var_bar_w_cutoff(self, all_bar_data, cond="pass_qc", cond_false_name="fail_qc", var_cat="my_var_annot", figsize=(5, 4), plt_add=""):
        top_y1 = all_bar_data[cond].max()
        bottom_y2 = all_bar_data[all_bar_data["all"].astype(int) > top_y1]["all"].astype(int).min()
        buffer = np.abs(top_y1 - bottom_y2)/2
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=figsize)
        fig.subplots_adjust(hspace=0.05, left=0.15)
        ax1.bar(all_bar_data.index, all_bar_data[cond].values, color="lightskyblue")
        ax1.bar(all_bar_data.index, all_bar_data[cond_false_name].values, bottom=all_bar_data[cond].values, color="peachpuff")
        ax2.bar(all_bar_data.index, all_bar_data[cond].values, color="lightskyblue")
        ax2.bar(all_bar_data.index, all_bar_data[cond_false_name].values, bottom=all_bar_data[cond].values, color="peachpuff")
        ax1.set_ylim(bottom_y2 - buffer*1.5, None)
        ax2.set_ylim(0, top_y1 + buffer/3)
        ax1.spines.bottom.set_visible(False)
        ax2.spines.top.set_visible(False)
        ax1.xaxis.tick_top()
        ax1.tick_params(labeltop=False)
        ax2.xaxis.tick_bottom()
        order = int(np.floor(np.log10(max(top_y1, all_bar_data["all"].astype(int).max()))))
        ax1.yaxis.set_major_formatter(FixedOrderFormatter(order))
        ax2.yaxis.set_major_formatter(FixedOrderFormatter(order))
        fig.supylabel("count")
        d = .5
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                    linestyle="none", color='k', mec='k', mew=1, clip_on=False)
        ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
        ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
        ax1.yaxis.get_offset_text().set_fontsize(7)
        ax2.yaxis.get_offset_text().set_visible(False)
        plt.xticks(rotation=90)
        pu.add_legend(fig, cdict={cond.replace("_", " "): "lightskyblue", cond_false_name.replace("_", " "): "peachpuff"}, bbox=pu.find_right_legend_loc(fig, [ax1, ax2], [ax1, ax2]), loc="center left")
        self.handle(f"{plt_add}{cond}_{var_cat}_bar_cutoff")

    def var_prop_bar(self, all_bar_data, cond="pass_qc", cond_false_name="fail_qc", figsize=(5, 3), plt_add=""):
        fig, ax = plt.subplots(figsize=figsize)
        r1 = ax.barh(y=all_bar_data.index, width=all_bar_data[f"{cond}_prop"], color="lightskyblue")
        r2 = ax.barh(y=all_bar_data.index, width=all_bar_data[f"{cond_false_name}_prop"], left=all_bar_data[f"{cond}_prop"], color="peachpuff")
        ax.bar_label(r1, labels=all_bar_data[cond].values.astype(int), label_type='edge', color="mediumblue", fmt="{:.0f}")
        ax.bar_label(r2, labels=all_bar_data[cond].values.astype(int), label_type='center', color="sienna", fmt="{:.0f}")
        ax.invert_yaxis()
        ax.margins(x=0)
        ax.set_xlabel("Proportion")
        pu.add_legend(fig, cdict={cond.replace("_", " "): "lightskyblue", cond_false_name.replace("_", " "): "peachpuff"}, bbox=pu.find_right_legend_loc(fig, [ax], [ax]), loc="center left")
        self.handle(f"{plt_add}{cond}_bar_props")

    def r2_sam_tab_plot(self, data, n_cell_name, r2_name):
        fig, ax = plt.subplots()
        ax.scatter(data[n_cell_name], data[r2_name], s=10, c="y", edgecolor="k")
        ax.set_xscale("log")
        ax.set_xlabel("Number of cells")
        ax.set_ylabel("Sample-wise R2 with genotyping")
        self.handle("n_cell_samwise_r2_scatter")


class FixedOrderFormatter(mpl.ticker.ScalarFormatter):
    # scalarformatter pinned to a fixed power of ten
    def __init__(self, order, useMathText=True):
        super().__init__(useMathText=useMathText)
        self._fixed_order = order
        self.set_scientific(True)
        self.set_useOffset(False)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self._fixed_order

def consequence_to_my_annot(cons_field, va_mapper, va_priority):
    cons = cons_field.split(",")
    vas = [va_mapper.get(c, "unknown") for c in cons]
    return str(np.array(va_priority)[np.isin(va_priority, vas)][0])


variant_collapse = {
    "splice_acceptor_variant": "essential splice",
    "splice_donor_variant": "essential splice",
    "stop_gained": "nonsense",
    "missense_variant": "missense",
    "synonymous_variant": "synonymous",
    "5_prime_UTR_variant": "5'UTR",
    "3_prime_UTR_variant": "3'UTR",
    "intron_variant": "intron",
    "intergenic_variant": "intergenic",
    "not_transcribed": "intergenic",
    "frameshift_variant": "other coding",
    "stop_lost": "other coding",
    "start_lost": "other coding",
    "transcript_ablation": "other coding",
    "transcript_amplification": "other coding",
    "feature_elongation": "other coding",
    "feature_truncation": "other coding",
    "inframe_insertion": "other coding",
    "inframe_deletion": "other coding",
    "protein_altering_variant": "other coding",
    "incomplete_terminal_codon_variant": "other coding",
    "start_retained_variant": "other coding",
    "stop_retained_variant": "other coding",
    "coding_sequence_variant": "other coding",
    "coding_transcript_variant": "other coding",
    "NMD_transcript_variant": "other coding",
    "splice_donor_5th_base_variant": "intron",
    "splice_region_variant": "intron",
    "splice_donor_region_variant": "intron",
    "splice_polypyrimidine_tract_variant": "intron",
    "mature_miRNA_variant": "other noncoding",
    "non_coding_transcript_exon_variant": "other noncoding",
    "non_coding_transcript_variant": "other noncoding",
    "upstream_gene_variant": "up/down gene",
    "downstream_gene_variant": "up/down gene",
    "TFBS_ablation": "regulatory",
    "TFBS_amplification": "regulatory",
    "TF_binding_site_variant": "regulatory",
    "regulatory_region_ablation": "regulatory",
    "regulatory_region_amplification": "regulatory",
    "regulatory_region_variant": "regulatory",
    "sequence_variant": "unknown",
}

va_priority = [
    "essential splice",
    "nonsense",
    "missense",
    "synonymous",
    "other coding",
    "5'UTR",
    "3'UTR",
    "intron",
    "other noncoding",
    "up/down gene",
    "regulatory",
    "intergenic",
    "unknown",
]
