process QUILT2 {
    stageInMode 'symlink'
    publishDir "${params.outdir}/quilt2_output_raw", mode: "copy"
    label 'quilt2_env'
    input:
        path bam_list
        tuple path(ref_panel), val(start), val(end), val(chr)

    output:
        tuple val(chr), val(start), val(end), path("*.vcf.gz"), emit: vcf

    script:
    """
    QUILT2.R \
    --prepared_reference_filename="${ref_panel}" \
    --bamlist="${bam_list}" \
    --method=diploid \
    --chr="${chr}" \
    --regionStart="${start}" \
    --regionEnd="${end}" \
    --nGen=100 \
    --buffer=500000 \
    --output_filename="${chr}.${start}.${end}.quilt2.diploid.output.vcf.gz"
    """
}
