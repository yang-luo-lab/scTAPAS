process FilterVariants {
    label 'sctapas_env'
    publishDir "${params.outdir}/sctapas_output", mode: 'copy'
    input:
        path vcf
        val infoscore
        val maf
        val hwe

    output:
        tuple path("${vcf.getSimpleName()}_info${infoscore}_maf${maf}_hwe${hwe}.vcf.gz"), path("${vcf.getSimpleName()}_info${infoscore}_maf${maf}_hwe${hwe}.vcf.gz.tbi")

    script:
    """
    bcftools view -i 'INFO/INFO_SCORE>${infoscore} & INFO/EAF>${maf} & INFO/EAF<${1-maf.toFloat()} & INFO/HWE>${hwe}' -Oz -o "${vcf.getSimpleName()}_info${infoscore}_maf${maf}_hwe${hwe}.vcf.gz" ${vcf.name}
    tabix -p vcf "${vcf.getSimpleName()}_info${infoscore}_maf${maf}_hwe${hwe}.vcf.gz"
    """
}
