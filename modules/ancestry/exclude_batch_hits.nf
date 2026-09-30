process ExcludeBatchHits {
    label 'sctapas_env'
    stageInMode 'symlink'
    input:
        tuple path(sctapas_pvar), path(sctapas_psam), path(sctapas_pgen)
        tuple path(onekg_pvar), path(onekg_psam), path(onekg_pgen)
        path assoc
        val thresh

    output:
        path "sig_ids.txt", emit: exbatch_vars
        tuple path("sctapas_exbatch.pvar"), path("sctapas_exbatch.psam"), path("sctapas_exbatch.pgen"), emit: sctapas
        tuple path("onekg_exbatch.pvar"), path("onekg_exbatch.psam"), path("onekg_exbatch.pgen"), emit: onekg
        path "*.log", emit: logs

    script:
    """
    awk 'NR==1{for(i=1;i<=NF;i++) c[\$i]=i; next}
       \$c["P"] != "NA" && \$c["P"] < ${thresh} {print \$c["ID"]}' ${assoc} > sig_ids.txt
    plink2 --pfile "${sctapas_pvar.baseName}" --output-chr chrM --exclude sig_ids.txt --make-pgen --out sctapas_exbatch
    plink2 --pfile "${onekg_pvar.baseName}" --output-chr chrM --exclude sig_ids.txt --make-pgen --out onekg_exbatch
    """
}
