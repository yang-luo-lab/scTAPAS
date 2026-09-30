process ExtractPositions {
    label 'sctapas_env'
    input:
        tuple path(vcf), val(start), val(end), val(chrom)

    output:
        path("${chrom}-${start}-${end}_${vcf.name}"),     emit: vcf
        path("${chrom}-${start}-${end}_${vcf.name}.tbi"), emit: tbi

    script:
    """
    tabix -p vcf ${vcf.name}
    bcftools view -r ${chrom}:${start}-${end} ${vcf.name} -Oz -o "${chrom}-${start}-${end}_${vcf.name}"
    tabix -p vcf "${chrom}-${start}-${end}_${vcf.name}"
    """
}
