process CisQTL {
    label 'sctapas_env'
    stageInMode 'symlink'
    publishDir "${params.outdir}/sctapas_tensorqtl/${identifier}", mode: 'copy'
    errorStrategy 'ignore'
    tag "${identifier}"
    input:
        path vcf
        tuple val(identifier), path(pseudobulk), path(covariates)
        val qval_lambda

    output:
        path "${identifier}.cis_qtl_pairs.*.parquet"
        path "${identifier}.tensorQTL.cis_nominal.log"
        path "${identifier}.cis_qtl.txt.gz", optional: true
        path "${identifier}.tensorQTL.cis.log", optional: true
        path "${identifier}.cis_independent_qtl.txt.gz", optional: true
        path "${identifier}.tensorQTL.cis_independent.log", optional: true

    script:
    """
    plink2 --vcf ${vcf} dosage=DS --output-chr chrM --keep-allele-order --make-bed --out input
    python3 -m tensorqtl input ${pseudobulk} ${identifier} --covariates ${covariates} --mode cis_nominal --maf_threshold 0.05
    python3 -m tensorqtl input ${pseudobulk} ${identifier} --covariates ${covariates} --mode cis --maf_threshold 0.05 --qvalue_lambda ${qval_lambda} || echo "WARNING: tensorQTL cis failed for ${identifier}" >&2
    [ -f "${identifier}.cis_qtl.txt.gz" ] && python3 -m tensorqtl input ${pseudobulk} ${identifier} --covariates ${covariates} --cis_output "${identifier}.cis_qtl.txt.gz" --mode cis_independent --qvalue_lambda ${qval_lambda} || echo "WARNING: tensorQTL cis_independent failed or skipped for ${identifier}" >&2
    """
}
