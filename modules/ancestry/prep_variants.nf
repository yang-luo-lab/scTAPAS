def setids = '@:#_\\$r_\\$a'

process PrepVariants {
    label 'sctapas_env'
    stageInMode 'symlink'
    input:
        path sctapas_vcf
        tuple path(onekg_pvar), path(onekg_psam), path(onekg_pgen)
        path rel_samples
        val mhc_from
        val mhc_to

    output:
        tuple path("sctapas_ovlp.pvar"), path("sctapas_ovlp.psam"), path("sctapas_ovlp.pgen"), emit: sctapas
        tuple path("onekg_unrel_ovlp.pvar"), path("onekg_unrel_ovlp.psam"), path("onekg_unrel_ovlp.pgen"), emit: onekg
        path "*.log", emit: logs

    script:
    """
    printf "chr6\\t${mhc_from}\\t${mhc_to}\\n" > mhc.range
    plink2 --vcf ${sctapas_vcf} --output-chr chrM --max-alleles 2 --min-alleles 2 --snps-only just-acgt --set-all-var-ids ${setids} --exclude bed1 mhc.range --make-pgen --out sctapas
    plink2 --pfile ${onekg_pvar.baseName} --output-chr chrM --set-all-var-ids ${setids} --allow-extra-chr --chr 1-22 --new-id-max-allele-len 487 --remove ${rel_samples} --make-pgen --out onekg_unrel
    comm -12 <( cut -f3 sctapas.pvar | sort ) <( cut -f3 onekg_unrel.pvar | sort ) > ovlp.txt
    plink2 --pfile sctapas --output-chr chrM --extract ovlp.txt --make-pgen --out sctapas_ovlp
    plink2 --pfile onekg_unrel --output-chr chrM --extract ovlp.txt --make-pgen --out onekg_unrel_ovlp
    """
}
