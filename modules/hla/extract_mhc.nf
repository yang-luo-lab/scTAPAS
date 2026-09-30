def setids    = '@:#_\\$r_\\$a'
def in_dosage  = 'dosage=DS'
def out_dosage = 'vcf-dosage=DS'

process FilterMHCVars {
    label 'sctapas_env'
    publishDir "${params.outdir}/sctapas_mhc", mode: 'copy'
    input:
        path vcf
        path mhc_ref_panel_vars

    output:
        path "${vcf.getSimpleName()}_mhc_ref_vars.vcf.gz",     emit: vcf

    script:
    """
    plink2 --vcf ${vcf} ${in_dosage} --output-chr chrM --set-all-var-ids ${setids} --extract ${mhc_ref_panel_vars} --export vcf bgz ${out_dosage} --out "${vcf.getSimpleName()}_mhc_ref_vars"
    """
}
