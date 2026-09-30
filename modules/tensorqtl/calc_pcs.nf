process CalcPCs {
    label 'sctapas_env'
    stageInMode 'symlink'
    publishDir "${params.outdir}/sctapas_tensorqtl", mode: 'copy'
    input:
        path sctapas_vcf
        val mhc_from
        val mhc_to

    output:
        path "sctapas_pcs.sscore"

    script:
    """
    printf "chr6\\t${mhc_from}\\t${mhc_to}\\n" > mhc.range
    plink2 --vcf ${sctapas_vcf} --output-chr chrM --exclude bed1 mhc.range --make-pgen --out sctapas_nomhc
    plink2 --pfile sctapas_nomhc --indep-pairwise 50 5 0.2 --bad-ld --output-chr chrM --out sctapas.prune
    plink2 --pfile sctapas_nomhc --extract sctapas.prune.prune.in --output-chr chrM --make-pgen --out sctapas_pruned
    plink2 --pfile sctapas_pruned --freq --bad-freqs --pca biallelic-var-wts --output-chr chrM --out pca_sctapas
    plink2 --pfile sctapas_pruned --read-freq pca_sctapas.afreq --score pca_sctapas.eigenvec.var 2 3 header-read no-mean-imputation variance-standardize --score-col-nums 5-14 --out sctapas_pcs
    """
}
