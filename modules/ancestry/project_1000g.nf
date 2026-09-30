process ProjectPCs {
    label 'sctapas_env'
    stageInMode 'symlink'
    publishDir "${params.outdir}/sctapas_ancestry", mode: 'copy'
    input:
        tuple path(sctapas_pvar), path(sctapas_psam), path(sctapas_pgen)
        tuple path(onekg_pvar), path(onekg_psam), path(onekg_pgen)

    output:
        path "sctapas_proj.sscore", emit: sctapas
        path "onekg_proj.sscore", emit: onekg
        path "*.log", emit: log

    script:
    """
    plink2 --pfile "${onekg_pvar.baseName}" --indep-pairwise 50 5 0.2 --output-chr chrM --out onekg.prune
    plink2 --pfile "${onekg_pvar.baseName}" --extract onekg.prune.prune.in --output-chr chrM --make-pgen --out onekg_pruned
    plink2 --pfile onekg_pruned --freq --pca biallelic-var-wts --output-chr chrM --out pca_onekg
    plink2 --pfile "${sctapas_pvar.baseName}" --extract onekg.prune.prune.in --read-freq pca_onekg.afreq --score pca_onekg.eigenvec.var 2 3 header-read no-mean-imputation variance-standardize --score-col-nums 5-14 --out sctapas_proj
    plink2 --pfile onekg_pruned --read-freq pca_onekg.afreq --score pca_onekg.eigenvec.var 2 3 header-read no-mean-imputation variance-standardize --score-col-nums 5-14 --out onekg_proj
    """
}
