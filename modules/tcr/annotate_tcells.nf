process AnnotateTCells {
    label 'sctapas_env'
    stageInMode 'copy'
    input:
        path tcr_tab
        path gex
        val tcr_barcode_col_name
        val ab_col_name
        val ab_col_value
        val sample_id_col_name

    output:
        path "abtcr_annot_tab.tsv"

    script:
        """
        python3 - <<'EOF'
import sys
import os
sys.path.append(os.environ["SCTAPAS_PYTHON"])
import sctapas.tcr_feat as tcr_feat
tcr_tab = tcr_feat.make_ann_tcr(
    "${tcr_tab}",
    "${gex}",
    tcr_barcode_col_name="${tcr_barcode_col_name}",
    ab_col_name="${ab_col_name}",
    ab_col_value="${ab_col_value}",
    sample_id_col_name="${sample_id_col_name}")
tcr_tab.to_csv("abtcr_annot_tab.tsv", sep="\t")
EOF
        """
}
