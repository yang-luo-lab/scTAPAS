import os
import functools
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from . import utils

# name_output sets this so handle() saves under the function name, not the name passed in
_forced_fname = None

def name_output(func):
    # save any plot made inside func under func's name, e.g. fig_3a -> fig_3a<ext>
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        global _forced_fname
        prev = _forced_fname
        _forced_fname = func.__name__
        try:
            return func(*args, **kwargs)
        finally:
            _forced_fname = prev
    return wrapper

class AbsPlot:
    def __init__(self, mode="show", sv_dir=None, res=600, ext=".png"):
        self.mode = mode
        self.ext = ext
        if self.mode == "save":
            if sv_dir is None:
                self.sv_dir = ""
            else:
                self.sv_dir = sv_dir
        else:
            self.sv_dir = None
        self.res = res

    def _get_sv_path(self, fn):
        return os.path.join(self.sv_dir, fn + self.ext)
    
    def handle(self, fn, tight=True):
        if _forced_fname is not None:
            fn = _forced_fname
        if self.mode == "show":
            # skip tight_layout if a layout engine is set, it errors once there's a colorbar
            if tight and plt.gcf().get_layout_engine() is None:
                plt.tight_layout()
            plt.show()
        elif self.mode == "save":
            fn = fn.replace(" ", "_")
            if tight:
                bbox_dict = {"bbox_inches": "tight"}
            else:
                bbox_dict = {}
            plt.savefig(self._get_sv_path(fn), dpi=self.res, transparent=True, **bbox_dict)

def add_legend(fig, cdict, bbox=(1.1, 0.5), loc='center right', frameon=False, title=""):
    rect = [mpl.patches.Patch(color=v, label=k) for k, v in cdict.items()]
    fig.legend(handles=rect, labels=cdict.keys(), loc=loc, borderaxespad=0, ncol=1, frameon=frameon, bbox_to_anchor=bbox, title=title)

def add_xy(ax, colour="lightgrey", linestyle="--", pad=0.05):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    lo = min(xmin, ymin)
    hi = max(xmax, ymax)
    pad = pad * (hi - lo)
    plt.plot([lo-pad, hi+pad], [lo-pad, hi+pad], c=colour, linestyle=linestyle, zorder=0)
    ax.set_xlim([lo, hi])
    ax.set_ylim([lo, hi])
    return ax

def find_right_legend_loc(fig, x_axes, y_axes, buffer=0.02):
    fig.canvas.draw()
    right_edge = max(ax.get_position().x1 for ax in x_axes)
    top_bottom = [[ax.get_position().y0, ax.get_position().y1] for ax in y_axes]
    scatter_ycenter = 0.5 * (np.amin(top_bottom) + np.amax(top_bottom))
    return right_edge + buffer, scatter_ycenter


def panels_2d_1dhist(x, xname, y, yname, bins, mode="scatter", colours=None, cname=None, cmap=None, cdict=None, xmin=-1, xmax=1, ymin=-1, ymax=1, s=2, alpha=0.8, hist1d_colour="k", label_renamer=None, figsize=(6,6)):
    fig, ax_scatter = plt.subplots(figsize=figsize, layout='constrained')
    divider = make_axes_locatable(ax_scatter)
    ax_histx = divider.append_axes("top", size="25%", pad=0.2, sharex=ax_scatter)
    ax_histy = divider.append_axes("right", size="25%", pad=0.2, sharey=ax_scatter)
    if cmap is None:
        if cdict is not None:
            colours = [cdict[c] for c in colours]
        else:
            colours = np.repeat(hist1d_colour, len(x))
    else:
        ax_cbar = divider.append_axes("bottom", size="6%", pad=0.6)
    if mode == "scatter":
        sc = ax_scatter.scatter(x, y, c=colours, cmap=cmap, s=s, alpha=alpha)
    elif mode == "hist":
        _, _, _, sc = ax_scatter.hist2d(x, y, cmap=cmap, bins=bins, norm=mpl.colors.LogNorm(vmin=1))
        if cname is None:
            cname = "count"
    if label_renamer is not None:
        xname = label_renamer.get(xname, xname)
        yname = label_renamer.get(yname, yname)
        cname = label_renamer.get(cname, cname)
    ax_scatter.set_ylabel(yname)
    ax_scatter.set_xlabel(xname)
    ax_histx.hist(x, bins=bins, color=hist1d_colour)
    ax_histx.set_ylabel("count")
    ax_histy.hist(y, bins=bins, orientation='horizontal', color=hist1d_colour)
    ax_histy.set_xlabel("count")
    if cmap is None and cdict is not None:
        add_legend(fig, cdict, title=cname, bbox=find_right_legend_loc(fig, [ax_scatter, ax_histy], [ax_scatter]), loc="center left")
    else:
        fig.colorbar(sc, cax=ax_cbar, label=cname, orientation="horizontal")
    ax_histx.tick_params(labelbottom=False)
    ax_histy.tick_params(labelleft=False)
    ax_histx.set_xlim(xmin, xmax)
    ax_histy.set_ylim(ymin, ymax)
    ax_scatter.set_xlim(xmin, xmax)
    ax_scatter.set_ylim(ymin, ymax)

def bin_func(y, digitized, bins, func):
    y_binned = [func(y[digitized==i]) for i in range(1, len(bins))]
    return y_binned

# hla/mhc manhattan plot helpers

def add_hla_fill_ax(hla_gene, a, ax):
    ax.axvspan(*utils.classical_gene_pos[hla_gene], color=utils.hla_def_cdict[hla_gene], alpha=a, label=f"HLA-{hla_gene}")
    return ax


def _manhattan(ax, gwas, xname, yname, hla_tag, palette, colour_default="lightgrey", s=2, a=0.8):
    ax.scatter(gwas[xname].values, gwas[yname], color=[palette.get(gn, colour_default) for gn in gwas[hla_tag].values], s=s, alpha=a)
    return ax


def accum_pos_over_chrom(plink_out, col_name="POScum"):
    chrom_mm = plink_out.groupby(by="CHROM")["POS"].agg(["min", "max"])
    chrom_mm["mm_diff"] = chrom_mm["max"] - chrom_mm["min"] + 1
    chrom_add = pd.Series(index=chrom_mm.index, data=np.cumsum([0, *chrom_mm["mm_diff"].values[:-1]]))
    plink_out[col_name] = np.nan
    for chrom in chrom_add.index:
        plink_out.loc[plink_out["CHROM"] == chrom, col_name] = plink_out["POS"][plink_out["CHROM"] == chrom] + chrom_add.loc[chrom]
    return plink_out


def add_hla_plink(plink_out, buffer=500, col_name="Gene"):
    plink_out["POS"] = plink_out["POS"].astype(int)
    plink_out[col_name] = "other"
    for g in utils.classical_genes:
        g_patts = ["_".join([var_type, g]) for var_type in ["SNP", "AA", "HLA"]]
        g_patt_all = "|".join(g_patts)
        is_gene = plink_out[plink_out["ID"].str.contains(g_patt_all)]
        if len(is_gene) == 0:
            continue
        gene_window = [is_gene["POS"].min(), is_gene["POS"].max()]
        gene_window_buffer = [gene_window[0] - buffer, gene_window[1] + buffer]
        hla_chr = is_gene["CHROM"].unique()
        below_upper = np.logical_and(np.isin(plink_out["CHROM"], hla_chr), plink_out["POS"] <= gene_window_buffer[1])
        above_lower = np.logical_and(np.isin(plink_out["CHROM"], hla_chr), plink_out["POS"] >= gene_window_buffer[0])
        plink_out.loc[np.logical_and(above_lower, below_upper), col_name] = g
    return plink_out


def aug_plink_out(plink_out, buffer=500):
    plink_out = add_hla_plink(plink_out, buffer=buffer, col_name="Gene")
    plink_out["nlog10P"] = -np.log10(plink_out["P"].values)
    plink_out = accum_pos_over_chrom(plink_out, col_name="POScum")
    plink_out["POSMbp"] = plink_out["POS"] / 1e6
    plink_out["POScumMbp"] = plink_out["POScum"] / 1e6
    return plink_out


def mhc_manhattan(fig, ax, gwas, ann_ids=None, ann_pos=None, ann_id_func=None, adjust=False, buffer=500):
    if adjust is True:
        from adjustText import adjust_text
    neat_ann_ids = [ann_id_func(aid) for aid in ann_ids] if (ann_ids is not None and ann_id_func is not None) else ann_ids
    gwas = aug_plink_out(gwas)
    ax = _manhattan(ax, gwas, "POScumMbp", "nlog10P", "Gene", utils.hla_def_cdict)
    if ann_ids is not None:
        if ann_pos is not None:
            [ax.annotate(n_ann_id, xy=tuple(gwas[gwas["ID"] == ann_id][["POSMbp", "nlog10P"]].to_numpy()[0]),
                 xytext=(ann_pos[i, 0], ann_pos[i, 1]), arrowprops=dict(arrowstyle='->', color='k'))
             for i, (ann_id, n_ann_id) in enumerate(zip(ann_ids, neat_ann_ids))]
        else:
            texts = [plt.text(*gwas[gwas["ID"] == ann_id][["POSMbp", "nlog10P"]].to_numpy()[0], n_ann_id)
                     for ann_id, n_ann_id in zip(ann_ids, neat_ann_ids)]
            if adjust:
                adjust_text(texts, arrowprops=dict(arrowstyle='->', color='k'),
                            ensure_inside_axes=True, ax=ax, expand=(-2, -2))
    add_legend(fig, utils.hla_def_cdict, bbox=find_right_legend_loc(fig, [ax], [ax]), loc="center left")
    ax.set_ylabel("-log10(P)")
    ax.set_xlabel("Chromosome 6 Position (Mbp)")
    ax.margins(y=0, x=0)
    return ax

def bars_either(fig, ax, data, catname, contname, cname, palette, sns_legend=False, vertical=True):
    if vertical:
        sns.barplot(data, y=contname, hue=cname, x=catname, palette=palette, legend=sns_legend)
        ax.set_ylim(0, 1)
        ax.set_ylabel(contname)
        ax.set_xlabel(None)
    else:
        sns.barplot(data, x=contname, hue=cname, y=catname, palette=palette, legend=sns_legend)
        ax.set_xlim(0, 1)
        ax.set_xlabel(contname)
        ax.set_ylabel(None)
    add_legend(fig, palette, bbox=find_right_legend_loc(fig, [ax], [ax]), loc="center left")