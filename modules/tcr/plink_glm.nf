process RunGLM {
    label 'sctapas_env'
    publishDir "${params.outdir}/sctapas_tcrqtl/plink_out", mode: "copy"
    errorStrategy 'ignore'
    input:
        tuple path(bim), path(bed), path(fam), path(pheno), path(cov)
        val ci

    output:
        path "*_HLA_multi_assoc.*"

    script:
    """
    plink2 --bfile ${bim.baseName} \
          --pheno ${pheno} \
          --covar ${cov} \
          --no-parents \
          --glm omit-ref hide-covar cols=chrom,pos,ref,alt,test,nobs,beta,se,ci,tz,p,a1freqcc,a1freq \
          --ci ${ci} \
          --out "${pheno.baseName}_HLA_multi_assoc" \
          --covar-variance-standardize
    """
}
