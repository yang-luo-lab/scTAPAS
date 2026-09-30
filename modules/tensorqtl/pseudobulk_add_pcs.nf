process PseudobulkCovAddPCs {
    label 'sctapas_env'
    stageInMode 'symlink'
    publishDir "${params.outdir}/sctapas_tensorqtl/${identifier}", mode: 'copy'
    tag "${identifier}"

    input:
        tuple val(identifier), path(covs)
        path pcs
        val n_pcs
        val insert_position

    output:
        tuple val(identifier), path("covs_pcs_${identifier}.txt")

    script:
    """
python3 - <<'EOF'
import sys
import os
import pandas as pd
cov = pd.read_table("${covs}", index_col=0).T
pcs = pd.read_table("${pcs}", index_col=0)
pc_names = [c for c in pcs.columns if "PC" in c]
pcs_keep = pcs[pc_names[:${n_pcs}]]
cov_w_pcs = cov.join(pcs_keep)
cov_w_pcs_ord = cov_w_pcs[[*cov.columns[:${insert_position}], *pcs_keep.columns, *cov.columns[${insert_position}:]]]
cov_w_pcs_ord.T.to_csv("covs_pcs_${identifier}.txt", sep="\t")
EOF
    """


}