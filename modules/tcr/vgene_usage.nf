process VGeneUsage {
    label 'sctapas_env'
    publishDir "${params.outdir}/sctapas_tcrqtl/pheno", mode: "copy"
    input:
        tuple val(subset_level), val(subset_name), val(gene), val(locus), path(tcr_tab)
        val min_cln
        val min_obs_prop
        val sample_id_col
        val clone_id_col

    output:
        path "${gene}_${locus}_${subset_level}_${subset_name.replace("/", "_")}_pheno.txt"

    script:
    """
    python3 - <<'EOF'
import sys
import os
sys.path.append(os.environ["SCTAPAS_PYTHON"])
import sctapas.tcr_feat as tcr_feat
import pandas as pd
tcr_tab = pd.read_table("${tcr_tab}", index_col=0)
tsub = list(["${subset_level}", "${subset_name}"])
by_clone = True
sams = tcr_tab["${sample_id_col}"].unique()
gtab = tcr_feat.get_vdj_usg(tcr_tab, sams, "${gene}", "${locus}", "${sample_id_col}", "${clone_id_col}", tsub, ${min_cln}, ${min_obs_prop}, by_clone)
gtab.to_csv("${gene}_${locus}_${subset_level}_${subset_name.replace("/", "_")}_pheno.txt", sep="\t", index=False)
EOF
    """
}
