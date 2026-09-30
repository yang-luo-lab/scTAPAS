process ConcatVCF {
    label 'sctapas_env'
    tag { name_add }
    input:
        tuple val(name_add), path(sample_vcfs), path(sample_vcf_indexes)

    output:
        path "${name_add}_merged_sorted.vcf.gz", emit: vcf
        path "${name_add}_merged_sorted.vcf.gz.tbi", emit: index

    script:
    """
    bcftools concat -Oz -a -o "${name_add}_merged.vcf.gz" ${sample_vcfs.join(' ')}
    bcftools sort "${name_add}_merged.vcf.gz" -Oz -o "${name_add}_merged_sorted.vcf.gz"
    tabix -p vcf "${name_add}_merged_sorted.vcf.gz"
    """
}
