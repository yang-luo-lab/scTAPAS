process RemovePublicTCRMultiplexArtifacts {
    label 'sctapas_env'
    stageInMode 'copy'
    input:
        path tcr_tab
        val sample_id_col_name
        val pool_id_col_name
        val trav_name
        val cdr3a_nt_name
        val traj_name
        val trbv_name
        val cdr3b_nt_name
        val trbj_name

    output:
        path "abtcr_annot_tab_nopubtcrmplx.tsv"

    script:
        """
        python3 - <<'EOF'
import sys
import os
sys.path.append(os.environ["SCTAPAS_PYTHON"])
import sctapas.tcr_feat as tcr_feat
import pandas as pd
tcr_tab = pd.read_table("${tcr_tab}")
tcr_tab = tcr_feat.demux_qc_clones(
    tcr_tab,
    sample_name_col="${sample_id_col_name}",
    pool_name="${pool_id_col_name}",
    trav_name="${trav_name}",
    CDR3a_nt_name="${cdr3a_nt_name}",
    traj_name="${traj_name}",
    trbv_name="${trbv_name}",
    CDR3b_nt_name="${cdr3b_nt_name}",
    trbj_name="${trbj_name}")
tcr_tab.to_csv("abtcr_annot_tab_nopubtcrmplx.tsv", sep="\t")
EOF
        """
}
