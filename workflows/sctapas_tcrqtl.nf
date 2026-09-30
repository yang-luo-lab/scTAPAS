include { PrepHLA                           } from '../modules/tcr/prep_hla_vcf'
include { VGeneUsage                        } from '../modules/tcr/vgene_usage'
include { RunGLM                            } from '../modules/tcr/plink_glm'
include { AnnotateTCells                    } from '../modules/tcr/annotate_tcells'
include { RemovePublicTCRMultiplexArtifacts } from '../modules/tcr/remove_public_tcr_mplex'

workflow TCRQTLs {
    take:
        sample_id_col_name
        clone_id_col_name
        min_cln
        min_obs_prop_tcr
        glm_ci
        tcr_barcode_col_name
        ab_col_name
        ab_col_value
        rem_pub_tcr_mplex_arts
        pool_id_col_name
        trav_name
        cdr3a_nt_name
        traj_name
        trbv_name
        cdr3b_nt_name
        trbj_name

    main:
        options_tuple = Channel.fromList(params.subset_options)
            .flatMap { level, names -> names.collect { name -> tuple(level, name) } }

        hla  = Channel.fromPath(params.hla_vcf)
        gex = Channel.fromPath(params.anndata)
        tcr_covs = Channel.value(file(params.tcr_covs))
        tcr  = Channel.fromPath(params.tcr_path)
        vdj  = Channel.of("v")
        loci = Channel.of("TRA", "TRB")

        PrepHLA(hla)

        AnnotateTCells(tcr, gex, tcr_barcode_col_name, ab_col_name, ab_col_value, sample_id_col_name)

        processed_tcr = rem_pub_tcr_mplex_arts
            ? RemovePublicTCRMultiplexArtifacts(AnnotateTCells.out, sample_id_col_name, pool_id_col_name, trav_name, cdr3a_nt_name, traj_name, trbv_name, cdr3b_nt_name, trbj_name)
            : AnnotateTCells.out

        usg_input = options_tuple.combine(vdj).combine(loci).combine(processed_tcr)
        VGeneUsage(usg_input, min_cln, min_obs_prop_tcr, sample_id_col_name, clone_id_col_name)

        // bfiles × pheno × tcr_covs
        glm_input = PrepHLA.out.combine(VGeneUsage.out).combine(tcr_covs)
        RunGLM(glm_input, glm_ci)
}

workflow {
    TCRQTLs(
        params.sample_id_col_name,
        params.clone_id_col_name,
        params.min_cln,
        params.min_obs_prop_tcr,
        params.glm_ci,
        params.tcr_barcode_col_name,
        params.ab_col_name,
        params.ab_col_value,
        params.rem_pub_tcr_mplex_arts,
        params.pool_id_col_name,
        params.trav_name,
        params.cdr3a_nt_name,
        params.traj_name,
        params.trbv_name,
        params.cdr3b_nt_name,
        params.trbj_name)
}
