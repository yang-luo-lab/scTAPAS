process FindBatchHits {
    label 'sctapas_env'
    stageInMode 'symlink'
    publishDir "${params.outdir}/sctapas_ancestry", mode: 'copy'
    input:
        tuple path(sctapas_pvar), path(sctapas_psam), path(sctapas_pgen)
        tuple path(onekg_pvar), path(onekg_psam), path(onekg_pgen)

    output:
        path "assoc.source.glm.logistic.hybrid", emit: assoc
        path "assoc.log", emit: log

    script:
    """
        plink2 --pfile "${sctapas_pvar.baseName}" --output-chr chrM --make-bed --out sctapas
        plink2 --pfile "${onekg_pvar.baseName}" --output-chr chrM --make-bed --out onekg
        plink --bfile onekg --bmerge sctapas --output-chr chrM --make-bed --out combined
        cat <(cat onekg.fam | awk -v OFS='\t' 'BEGIN{print "#IID", "source"} {print \$2, 1}')  <(cat sctapas.fam | awk -v OFS='\t' '{print \$2, 2}') > pheno.txt
        plink2 --bfile combined --ci 0.95 --glm omit-ref allow-no-covars cols=chrom,pos,ref,alt,test,nobs,beta,se,ci,tz,p,a1freqcc,a1freq --pheno pheno.txt --output-chr chrM --out assoc
    """
}
