include { FilterMHCVars } from '../modules/hla/extract_mhc'

workflow ExtractMHC {
    take:
        vcf

    main:
        FilterMHCVars(vcf, file(params.mhc_ref_panel_vars))

    emit:
        vcf = FilterMHCVars.out.vcf
}

// Standalone entry: point params.vcf at any VCF from sctapas_quilt2
workflow {
    ExtractMHC(Channel.fromPath(params.vcf))
}
