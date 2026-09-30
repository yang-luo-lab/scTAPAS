process PrepHLA {
    label 'sctapas_env'
    publishDir "${params.outdir}/sctapas_tcrqtl/geno_cov", mode: "copy"
    input:
        path vcf

    output:
        tuple path("${vcf.baseName - '.vcf'}.for_tcrqtl.bim"), path("${vcf.baseName - '.vcf'}.for_tcrqtl.bed"), path("${vcf.baseName - '.vcf'}.for_tcrqtl.fam")

    script:
    """
    plink2 --vcf ${vcf} --double-id --maf 0.05 --make-bed --out "${vcf.baseName - '.vcf'}.for_tcrqtl"
    """
}
