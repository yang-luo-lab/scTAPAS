process Prep1000G {
    stageInMode 'symlink'
    label 'sctapas_env'
    input:
        path vcf
        path vcf_ind
        val chrom

    output:
        tuple val(chrom), path("*.noNA12878.vcf.gz"), path("*.noNA12878.vcf.gz.*")

    script:
    def rename = params.num_chrs
        ? """plink2 --vcf "1000GP.${chrom}.noNA12878.vcf.gz" --output-chr M --keep-allele-order --export vcf bgz --out "1000GP.${chrom}.noNA12878_tmp"
             mv "1000GP.${chrom}.noNA12878_tmp.vcf.gz" "1000GP.${chrom}.noNA12878.vcf.gz" """
        : ""
    """
    bcftools norm -m -any ${vcf} -Ou | bcftools view -m 2 -M 2 -v snps -s ^NA12878,NA12891,NA12892 -Oz -o "1000GP.${chrom}.noNA12878.vcf.gz"
    ${rename}
    bcftools index -f "1000GP.${chrom}.noNA12878.vcf.gz"
    """
}
