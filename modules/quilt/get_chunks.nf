process QUILT_Chunks {
    publishDir params.outdir, mode: "copy"
    label 'quilt2_env'
    input:
        path genetic_map
        val chrom

    output:
        tuple val(chrom), path("chunks_${chrom}.csv")

    script:
    """
    zcat $genetic_map > "${genetic_map.baseName}"
    Rscript - "${chrom}" "${genetic_map.baseName}" <<'RSCRIPT'
    library('QUILT')
    args <- commandArgs(trailingOnly = TRUE)
    chr <- args[1]
    genetic_map_file <- args[2]

    chunks <- QUILT::quilt_chunk_map(chr = chr, genetic_map_file = genetic_map_file)
    write.csv(chunks\$region,
      file = paste0('fromR_chunks_', chr, '.csv'),
      row.names = FALSE,
      quote = FALSE)
    RSCRIPT
    cat "fromR_chunks_${chrom}.csv" | tail -n +2 | tr ':' '-' | tr '-' ',' > "chunks_${chrom}.csv"
    """
}
