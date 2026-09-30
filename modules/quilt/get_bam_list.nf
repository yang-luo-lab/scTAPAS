process BamList {
    label 'sctapas_env'
    input:
        val chrom

    output:
        tuple val(chrom), path("*.txt")

    script:
    def pattern = params.bam_pattern.replace("CHROM", "${chrom}")
    """
    ls ${pattern} > bam_list_${chrom}.txt
    """
}
