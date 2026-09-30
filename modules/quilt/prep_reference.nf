process QUILT2_PrepRef {
    stageInMode 'symlink'
    label 'quilt2_env'
    input:
        path genetic_map
        path reference_vcf
        path reference_vcf_index
        val start
        val end
        val chr

    output:
        tuple path("reference_panel/RData/*.RData"), val(start), val(end), val(chr)

    script:
    """
    QUILT2_prepare_reference.R \
    --genetic_map_file="${genetic_map}" \
    --reference_vcf_file="${reference_vcf}" \
    --chr="${chr}" \
    --regionStart="${start}" \
    --regionEnd="${end}" \
    --nGen=100 \
    --buffer=500000 \
    --outputdir="reference_panel"
    """
}
