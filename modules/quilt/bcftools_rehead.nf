process ReheadVCF {
    label 'sctapas_env'
    input:
        path vcf
        path bam_list
        val chrom

    output:
        tuple val(chrom), path("${vcf.name.replace('.output.vcf.gz', '.vcf.gz')}"), emit: chrom_vcf
        path "${vcf.name.replace('.output.vcf.gz', '.vcf.gz')}.tbi",                emit: index

    script:
    def outname = vcf.name.replace('.output.vcf.gz', '.vcf.gz')
    """
    awk '${params.sample_name_awk}' ${bam_list} > sample_names.txt
    bcftools reheader -s sample_names.txt ${vcf} -o "${outname}"
    tabix -p vcf "${outname}"
    """
}
