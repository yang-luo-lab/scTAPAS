include { CisQTL              } from '../modules/tensorqtl/cis_qtl'
include { CalcPCs             } from '../modules/tensorqtl/calc_pcs'
include { PseudobulkCovAddPCs } from '../modules/tensorqtl/pseudobulk_add_pcs'

workflow TensorQTL {
    take:
        vcf

    main:
        CalcPCs(vcf, params.mhc_from, params.mhc_to)
        id_pb_channel = Channel.fromList(params.identifiers).map { id ->
            tuple(id, file(params.pseudobulk_pattern.replace("IDENTIFIER", id)))}
        id_cov_channel = Channel.fromList(params.identifiers).map { id ->
            tuple(id, file(params.covariates_pattern.replace("IDENTIFIER", id)))}
        PseudobulkCovAddPCs(id_cov_channel, CalcPCs.out.first(), params.n_pcs_tqtl, params.insert_position_pcs_tqtl)
        id_channel = id_pb_channel.combine(PseudobulkCovAddPCs.out, by: 0)
        CisQTL(vcf.first(), id_channel, params.qval_lambda)
}

// Standalone entry: point params.vcf at any post-QC VCF
workflow {
    TensorQTL(Channel.fromPath(params.vcf))
}
