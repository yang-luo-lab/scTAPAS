import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from . import plot_utils as pu

def neaten_pc(pc_name):
    return pc_name.split("_")[0]

def half_marker(side, n=50):
    # closing the arc along the diameter gives a filled semicircle
    a0, a1 = (90, 270) if side == 'left' else (-90, 90)
    t = np.radians(np.linspace(a0, a1, n))
    verts = np.column_stack([np.cos(t), np.sin(t)])
    verts = np.vstack([verts, verts[0]])          # close along the flat edge
    codes = [mpl.path.Path.MOVETO] + [mpl.path.Path.LINETO] * (n - 1) + [mpl.path.Path.CLOSEPOLY]
    return mpl.path.Path(verts, codes)
    
def circle_marker(n=120):
    t = np.radians(np.linspace(0, 360, n))
    verts = np.column_stack([np.cos(t), np.sin(t)])
    verts = np.vstack([verts, verts[0]])
    codes = [mpl.path.Path.MOVETO] + [mpl.path.Path.LINETO] * (n - 1) + [mpl.path.Path.CLOSEPOLY]
    return mpl.path.Path(verts, codes)

class AncestryPCA(pu.AbsPlot):
    def __init__(self, mode="show", sv_dir=None, res=600, ext=".png"):
        super().__init__(mode=mode, sv_dir=sv_dir, res=res, ext=ext)
        
    def scatter_target_ref_pca(self, pca_target, pca_ref, pca_colname_pair, cdict, anc_target, anc_ref, ref_marker="o", target_marker="o", ref_alpha=0.1, target_edgecolor="k", figsize=(5, 4), plt_add=""):
        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(pca_ref[pca_colname_pair[0]], pca_ref[pca_colname_pair[1]], c=[cdict.get(c, "k") for c in anc_ref], marker=ref_marker, alpha=ref_alpha)
        ax.scatter(pca_target[pca_colname_pair[0]],pca_target[pca_colname_pair[1]], c=[cdict.get(c, "k") for c in anc_target], marker=target_marker, edgecolor=target_edgecolor)
        ax.set_xlabel(neaten_pc(pca_colname_pair[0]))
        ax.set_ylabel(neaten_pc(pca_colname_pair[1]))
        pu.add_legend(fig, {key: val for key, val in cdict.items() if key in [*anc_target, *anc_ref]}, loc='center left', bbox=pu.find_right_legend_loc(fig, [ax], [ax], buffer=0.02))
        self.handle(f"{plt_add}target_ref{pca_colname_pair[0]}_{pca_colname_pair[1]}_scatter", tight=True)

    def scatter_target_ref_pca_split(self, pca_target, pca_ref, pca_colname_pair, cdict, anc_target_left, anc_target_right, anc_ref, ref_marker="o", target_marker="o", ref_alpha=0.1, target_edgecolor="k", figsize=(5, 4), plt_add=""):
        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(pca_ref[pca_colname_pair[0]], pca_ref[pca_colname_pair[1]], c=[cdict.get(c, "k") for c in anc_ref], marker=ref_marker, alpha=ref_alpha, s=36)
        ms = 36
        common = dict(s=ms, linewidths=0)
        for i in range(pca_target.shape[0]):
            s0 = ax.scatter(pca_target[pca_colname_pair[0]].iloc[i],pca_target[pca_colname_pair[1]].iloc[i], s=ms*1.4, marker=circle_marker(), facecolors=target_edgecolor, linewidths=0, zorder=2)
            s1 = ax.scatter(pca_target[pca_colname_pair[0]].iloc[i],pca_target[pca_colname_pair[1]].iloc[i], marker=half_marker('left'),  c=cdict.get(anc_target_left[i], "k"), edgecolor=target_edgecolor, **common, zorder=2)
            s2 = ax.scatter(pca_target[pca_colname_pair[0]].iloc[i],pca_target[pca_colname_pair[1]].iloc[i], marker=half_marker('right'), c=cdict.get(anc_target_right[i], "k"), edgecolor=target_edgecolor, **common, zorder=2)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        ax.set_xlabel(neaten_pc(pca_colname_pair[0]))
        ax.set_ylabel(neaten_pc(pca_colname_pair[1]))
        cleg_pos = pu.find_right_legend_loc(fig, [ax], [ax], buffer=0.02)
        pu.add_legend(fig, {key: val for key, val in cdict.items() if key in [*anc_target_left, *anc_target_right, *anc_ref]}, loc='center left', bbox=[cleg_pos[0], cleg_pos[1]+0.2])
        marker_handles = [mpl.lines.Line2D([0], [0], marker=target_marker, linestyle="none", markerfacecolor="darkgray", markeredgecolor=target_edgecolor, markersize=np.sqrt(ms), label="scTAPAS"),
                          mpl.lines.Line2D([0], [0], marker=ref_marker, linestyle="none", markerfacecolor="darkgray", markeredgecolor="darkgray", markersize=np.sqrt(ms), label="1000 Genomes"),
                         mpl.lines.Line2D([0], [0], marker=half_marker('left'), linestyle="none", markerfacecolor="darkgray", markeredgecolor="darkgray",  markersize=np.sqrt(ms), label="scTAPAS ancestry"),
                         mpl.lines.Line2D([0], [0], marker=half_marker('right'), linestyle="none", markerfacecolor="darkgray", markeredgecolor="darkgray", markersize=np.sqrt(ms), label="Array-based ancestry"),]
        fig.legend(handles=marker_handles, loc="center left", bbox_to_anchor=(cleg_pos[0], cleg_pos[1]-0.2), title=None, frameon=False)
        self.handle(f"{plt_add}target_ref{pca_colname_pair[0]}_{pca_colname_pair[1]}_split_scatter", tight=True)

def anc_assign(pca_ref, pca_target, pc_colnames, ref_pop_colname, sam_colname):
    sd_ref = pca_ref[pc_colnames].std()
    ref_pcs = pca_ref[pc_colnames]/sd_ref
    ref_pcs = ref_pcs.join(pca_ref[[ref_pop_colname]])
    target_pcs = pca_target[pc_colnames]/sd_ref
    centroids = ref_pcs.groupby(ref_pop_colname)[pc_colnames].mean()
    ds = []
    for s in target_pcs.index:
        ds.append(np.abs(target_pcs.loc[s]-centroids).sum(axis=1))
    d = pd.concat(ds, axis=1).T
    d = d.join(pca_target[sam_colname]).set_index(sam_colname)
    return d.idxmin(axis=1)

class scTAPASAncestry:
    def __init__(self, ref_path, target_path, n_pcs_plink=10, ref_pop_colname="SuperPop", id_colname="#IID", cdict=None, mode="show", sv_dir=None, res=600, ext=".png"):
        self.pca_ref = pd.read_table(ref_path)
        self.pca_target = pd.read_table(target_path)
        self.pcs = self.pca_ref.columns[-n_pcs_plink:]
        self.ref_pop_colname = ref_pop_colname
        self.id_colname = id_colname
        self.plots = AncestryPCA(mode=mode, sv_dir=sv_dir, res=res, ext=ext)
        self.cdict = cdict if cdict is not None else dict(zip(['EUR', 'SAS', 'EAS', 'AFR', 'AMR'], ['r', 'b', 'g', 'y', 'm']))
        self.anc = None
        self.comp_anc = None
        
    def assign_anc(self, n_pcs):
        anc = anc_assign(self.pca_ref, self.pca_target, self.pcs[:n_pcs], self.ref_pop_colname, self.id_colname)
        if self.anc is None:
            self.anc = anc
        return anc

    def pca_plots(self, plot_upto_pc=8, ref_marker="o", target_marker="o", ref_alpha=0.1, target_edgecolor="k", figsize=(5, 4), plt_name_add=""):
        for i in range(plot_upto_pc-1):
            self.plots.scatter_target_ref_pca(self.pca_target, self.pca_ref, [self.pcs[i], self.pcs[i+1]], self.cdict, self.anc, self.pca_ref[self.ref_pop_colname], ref_marker=ref_marker, target_marker=target_marker, ref_alpha=ref_alpha, target_edgecolor=target_edgecolor, figsize=figsize, plt_add=plt_name_add)
