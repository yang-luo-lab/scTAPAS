def setids    = '@:#_\\$r_\\$a'
def in_dosage  = 'dosage=DS'
def out_dosage = 'vcf-dosage=DS'

process NormVariants {
    label 'sctapas_env'
    publishDir "${params.outdir}/quilt2_output_merged", mode: 'copy'
    input:
        path vcf
        path ref_fasta

    output:
        path "${vcf.name.replace('.vcf.gz', '')}_norm.vcf.gz"

    script:
    def base       = vcf.name.replace('.vcf.gz', '')
    def outchr_norm = params.num_chrs ? 'M' : 'chrM'
    """
    plink2 --vcf ${vcf.name} ${in_dosage} --output-chr ${outchr_norm} --keep-allele-order --export vcf bgz ${out_dosage} --out "${base}_chrfix"
    bcftools norm -m - -f "${ref_fasta}" -c s "${base}_chrfix.vcf.gz" -Oz -o "${base}_bcfnorm.vcf.gz"
    plink2 --vcf "${base}_bcfnorm.vcf.gz" ${in_dosage} --output-chr chrM --set-all-var-ids ${setids} --export vcf bgz ${out_dosage} --out "${base}_norm"
    """
}
