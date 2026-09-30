include { QUILT2_PrepRef    } from '../modules/quilt/prep_reference'
include { QUILT2            } from '../modules/quilt/impute'
include { BamList           } from '../modules/quilt/get_bam_list'
include { ReheadVCF         } from '../modules/quilt/bcftools_rehead'
include { QUILT_Chunks      } from '../modules/quilt/get_chunks'
include { Prep1000G         } from '../modules/quilt/prep_1000g_vcf'
include { ChunkLimits       } from '../modules/quilt/chunk_limits'
include { ConcatVCF         } from '../modules/quilt/concat_vcf'
include { ExtractPositions  } from '../modules/quilt/extract_positions'
include { NormVariants      } from '../modules/quilt/fix_variants'
include { FilterVariants    } from '../modules/quilt/choose_variants'
include { ExtractMHC        } from './sctapas_hla'
include { AncestryInference } from './sctapas_ancestry'
include { TensorQTL         } from './sctapas_tensorqtl'

// reference files always use chr prefix, convert numerical chr for lookups
def chrM(chr)     { params.num_chrs ? "chr${chr}" : "${chr}" }
def gmap(chr)     { file(params.genetic_map_pattern.replace("CHROM", chrM(chr))) }
def thouG(chr)    { file(params.thouG_vcf_pattern.replace("CHROM", chrM(chr))) }
def thouGInd(chr) { file(params.thouG_vcf_pattern.replace("CHROM", chrM(chr)) + ".tbi") }

workflow QUILT2_impute {
    take:
        start_pos
        end_pos
        chrom
        ref_vcf      // value channel, callers pass .first() so it broadcasts to all chunks
        ref_vcf_ind
        tagged_bam   // (chr, start, end, bam), join by key so chunks don't mis-pair

    main:
        genetic_map = chrom.map { chr -> gmap(chr) }
        QUILT2_PrepRef(genetic_map, ref_vcf, ref_vcf_ind, start_pos, end_pos, chrom)
        // join ref panel output with bam by (chr, start, end)
        quilt2_in = QUILT2_PrepRef.out
            .map { rd, s, e, chr -> tuple(chr, s, e, rd) }
            .join(tagged_bam, by: [0, 1, 2])
        // quilt2_in: (chr, start, end, rdata, bam)
        QUILT2(
            quilt2_in.map { chr, s, e, rd, bam -> bam },
            quilt2_in.map { chr, s, e, rd, bam -> tuple(rd, s, e, chr) }
        )
        // QUILT2.out.vcf: (chr, start, end, vcf), join with quilt2_in's bam by same key
        rehead_in = QUILT2.out.vcf
            .join(quilt2_in.map { chr, s, e, rd, bam -> tuple(chr, s, e, bam) }, by: [0, 1, 2])
        ReheadVCF(
            rehead_in.map { chr, s, e, vcf, bam -> vcf },
            rehead_in.map { chr, s, e, vcf, bam -> bam },
            rehead_in.map { chr, s, e, vcf, bam -> chr }
        )

    emit:
        chrom_vcf = ReheadVCF.out.chrom_vcf
}

workflow QUILT2_merge {
    // single input: (chrom, csv, vcf) tuples joined by chrom at the call site
    take:
        chrom_csv_vcf

    main:
        // dedup to one (chrom, csv) per chromosome for ChunkLimits
        chrom_csv = chrom_csv_vcf.map { chr, csv, vcf -> tuple(chr, csv) }
            .groupTuple()
            .map { chr, csvs -> tuple(chr, csvs[0]) }
        ChunkLimits(chrom_csv)
        bounds = ChunkLimits.out.splitCsv(header: true)
            .map { row -> tuple(row.vcf_name, row.start_bound.toInteger(), row.end_bound.toInteger(), row.chrom) }
        // join VCFs to bounds by filename
        merge_bounds = chrom_csv_vcf.map { chr, csv, vcf -> tuple(vcf.name, vcf) }
            .join(bounds)
            .map { name, vcf, start, end, chr -> tuple(vcf, start, end, chr) }
        ExtractPositions(merge_bounds)
        concat_input = ExtractPositions.out.vcf.collect()
            .map  { vcfs -> tuple("all", vcfs) }
            .join (ExtractPositions.out.tbi.collect().map { tbis -> tuple("all", tbis) })
        ConcatVCF(concat_input)

    emit:
        vcf = ConcatVCF.out.vcf
}

workflow QUILT2_QC {
    take:
        vcf

    main:
        NormVariants(vcf, params.ref_fasta)
        FilterVariants(NormVariants.out, 0.8, 0.05, 1e-6)

    emit:
        vcf = FilterVariants.out
}

workflow scTAPAS_region {
    take:
        chrom
        bp_from
        bp_to

    main:
        Prep1000G(chrom.map { thouG(it) }, chrom.map { thouGInd(it) }, chrom)
        bam_list = BamList(chrom)
        // tagged_bam: (chr, start, end, bam) so QUILT2_impute can join by key
        tagged_bam = chrom.combine(bp_from).combine(bp_to)
            .combine(bam_list.map { chr, bam -> bam }.first())
        QUILT2_impute(
            bp_from, bp_to, chrom,
            Prep1000G.out.map { chr, vcf, ind -> vcf }.first(),
            Prep1000G.out.map { chr, vcf, ind -> ind }.first(),
            tagged_bam
        )
        QUILT2_QC(QUILT2_impute.out.chrom_vcf.map { chr, vcf -> vcf })

    emit:
        vcf = QUILT2_QC.out.vcf
}

workflow scTAPAS_MHC {
    scTAPAS_region(
        Channel.of(params.num_chrs ? "6" : "chr6"),
        params.mhc_from,
        params.mhc_to
    )
    if (params.prep_mhc_imp) ExtractMHC(scTAPAS_region.out.vcf.map { vcf, tbi -> vcf })
}

workflow scTAPAS_REGION {
    scTAPAS_region(
        Channel.of(params.chrom),
        params.region_from,
        params.region_to
    )
}

workflow scTAPAS_CHROM {
    chrom_ch = Channel.value(params.chrom)
    Prep1000G(chrom_ch.map { thouG(it) }, chrom_ch.map { thouGInd(it) }, chrom_ch)
    bam_list = BamList(chrom_ch)
    QUILT_Chunks(gmap(params.chrom), params.chrom)
    rows = QUILT_Chunks.out.flatMap { chr, csv -> csv.splitCsv(header: false).collect { row -> tuple(chr, row[1].toInteger(), row[2].toInteger()) } }
    QUILT2_impute(
        rows.map { it[1] }, rows.map { it[2] }, rows.map { it[0] },
        Prep1000G.out.map { chr, vcf, ind -> vcf }.first(),
        Prep1000G.out.map { chr, vcf, ind -> ind }.first(),
        rows.combine(bam_list.map { chr, bam -> bam }.first())
    )
    chrom_csv_vcf = QUILT_Chunks.out.combine(QUILT2_impute.out.chrom_vcf, by: 0)
    QUILT2_merge(chrom_csv_vcf)
    QUILT2_QC(QUILT2_merge.out.vcf)
    if (params.prep_mhc_imp)   ExtractMHC(QUILT2_QC.out.vcf.map { vcf, tbi -> vcf })
    if (params.infer_ancestry) AncestryInference(QUILT2_QC.out.vcf.map { vcf, tbi -> vcf })
    if (params.run_tensorqtl)  TensorQTL(QUILT2_QC.out.vcf.map { vcf, tbi -> vcf })
}

workflow scTAPAS_chroms {
    take:
        chroms  // channel of chromosome names to process

    main:
        // one Prep1000G and BamList job per chromosome
        Prep1000G(chroms.map { chr -> thouG(chr) }, chroms.map { chr -> thouGInd(chr) }, chroms)
        BamList(chroms)
        QUILT_Chunks(chroms.map { chr -> gmap(chr) }, chroms)
        // fan out per-chunk rows for imputation
        rows = QUILT_Chunks.out.flatMap { chr, csv ->
            csv.splitCsv(header: false).collect { row -> tuple(chr, row[1].toInteger(), row[2].toInteger()) }
        }
        // combine per-chunk rows with per-chromosome ref prep by chrom key
        // rows: (chrom, start, end)  Prep1000G.out: (chrom, ref_vcf, ref_vcf_ind)
        chunk_prep = rows.combine(Prep1000G.out, by: 0)
        // tagged_bam: (chr, start, end, bam), pairs each chunk with its chromosome's bam list by key
        tagged_bam = rows.combine(BamList.out, by: 0)
        QUILT2_impute(
            chunk_prep.map { chr, s, e, rv, ri -> s   },
            chunk_prep.map { chr, s, e, rv, ri -> e   },
            chunk_prep.map { chr, s, e, rv, ri -> chr },
            chunk_prep.map { chr, s, e, rv, ri -> rv  },
            chunk_prep.map { chr, s, e, rv, ri -> ri  },
            tagged_bam
        )
        // join chrom_csv with chrom_vcf by chrom to get (chrom, csv, vcf) tuples
        chrom_csv_vcf = QUILT_Chunks.out.combine(QUILT2_impute.out.chrom_vcf, by: 0)
        QUILT2_merge(chrom_csv_vcf)
        QUILT2_QC(QUILT2_merge.out.vcf)
        if (params.prep_mhc_imp)   ExtractMHC(QUILT2_QC.out.vcf.map { vcf, tbi -> vcf })
        if (params.infer_ancestry) AncestryInference(QUILT2_QC.out.vcf.map { vcf, tbi -> vcf })
        if (params.run_tensorqtl)  TensorQTL(QUILT2_QC.out.vcf.map { vcf, tbi -> vcf })

    emit:
        vcf = QUILT2_QC.out.vcf
}

workflow scTAPAS_CHROMS {
    // chroms passed in params.chroms, e.g. ["chr6", "chr22"]
    scTAPAS_chroms(Channel.fromList(params.chroms))
}

workflow scTAPAS_AUTO {
    scTAPAS_chroms(Channel.of(1..22).map { n -> "${params.num_chrs ? '' : 'chr'}${n}" })
}
