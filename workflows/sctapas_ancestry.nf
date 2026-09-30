include { PrepVariants     } from '../modules/ancestry/prep_variants'
include { FindBatchHits    } from '../modules/ancestry/find_batch_hits'
include { ExcludeBatchHits } from '../modules/ancestry/exclude_batch_hits'
include { ProjectPCs       } from '../modules/ancestry/project_1000g'
include { InferAncestry    } from '../modules/ancestry/infer_ancestry'

workflow AncestryInference {
    take:
        vcf

    main:
        onekg    = Channel.of(params.onek_gen_patt).map { p -> tuple(file("${p}.pvar"), file("${p}.psam"), file("${p}.pgen")) }
        rel_sams = Channel.of(file(params.onek_rel_sams))
        PrepVariants(vcf, onekg, rel_sams, params.mhc_from, params.mhc_to)
        FindBatchHits(PrepVariants.out.sctapas, PrepVariants.out.onekg)
        ExcludeBatchHits(PrepVariants.out.sctapas, PrepVariants.out.onekg, FindBatchHits.out.assoc, 1e-6)
        ProjectPCs(ExcludeBatchHits.out.sctapas, ExcludeBatchHits.out.onekg)
        InferAncestry(ProjectPCs.out.sctapas, ProjectPCs.out.onekg, params.pcs_to_use)

    emit:
        ancestry = InferAncestry.out.assignments
}

// Standalone entry: point params.vcf at any VCF
workflow {
    AncestryInference(Channel.fromPath(params.vcf))
}
