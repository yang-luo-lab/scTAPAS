process InferAncestry {
    label 'sctapas_env'
    stageInMode 'symlink'
    publishDir "${params.outdir}/sctapas_ancestry", mode: 'copy'
    input:
        path sctapas_pcs
        path onekg_pcs
        val  pcs_to_use

    output:
        path "sctapas_pop_assign.tsv", emit: assignments
        path "*.png",                  emit: plots

    script:
    """
    python3 - <<'EOF'
import sys
import os
sys.path.append(os.environ['SCTAPAS_PYTHON'])
import pandas as pd
import sctapas.ancestry
sct_anc = sctapas.ancestry.scTAPASAncestry("${onekg_pcs}", "${sctapas_pcs}", mode="save", sv_dir=".")
anc = sct_anc.assign_anc(${pcs_to_use})
sct_anc.pca_plots(plot_upto_pc=${pcs_to_use})
anc.to_frame(name="scTAPAS_pop_assign").to_csv("sctapas_pop_assign.tsv", sep="\t")
EOF
    """
}