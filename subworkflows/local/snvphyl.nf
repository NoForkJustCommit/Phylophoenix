//
// Subworkflow: Running SNVPhyl from --> https://github.com/DHQP/SNVPhyl_Nextflow/tree/nf-core
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_INPUT                } from '../../modules/local/convert_input'
include { INDEXING                     } from '../../modules/local/indexing'
include { FIND_REPEATS                 } from '../../modules/local/find_repeats'
include { SMALT_MAP                    } from '../../modules/local/smalt_map'
include { SORT_INDEX_BAMS              } from '../../modules/local/sort_index_bams'
include { GENERATE_LINE_1              } from '../../modules/local/generate_line_1'
include { VERIFYING_MAP_Q              } from '../../modules/local/verifying_map_q'
include { FREEBAYES                    } from '../../modules/local/freebayes'
include { FILTER_FREEBAYES             } from '../../modules/local/filter_freebayes'
include { BGZIP_FREEBAYES_VCF          } from '../../modules/local/bgzip_freebayes_vcf'
include { FREEBAYES_VCF_TO_BCF         } from '../../modules/local/freebayes_vcf_to_bcf'
include { MPILEUP                      } from '../../modules/local/mpileup'
include { BGZIP_MPILEUP_VCF            } from '../../modules/local/bgzip_mpileup_vcf'
include { BCFTOOLS_CALL                } from '../../modules/local/bcftools_call'
include { CONSOLIDATE_BCFS             } from '../../modules/local/consolidate_bcfs'
include { CONSOLIDATE_FILTERED_DENSITY } from '../../modules/local/consolidate_filtered_density'
include { GENERATE_LINE_2              } from '../../modules/local/generate_line_2'
include { FILTER_STATS                 } from '../../modules/local/filter_stats'
include { VCF2SNV_ALIGNMENT            } from '../../modules/local/vcf2snv_alignment'
include { PHYML                        } from '../../modules/local/phyml'
include { MAKE_SNV                     } from '../../modules/local/make_snv'

/*
========================================================================================
    IMPORT NF-CORE SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK                  } from './input_check'
include { FASTQC                       } from '../../modules/nf-core/fastqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def add_st(input) {
    // create meta map
    def meta = [:]
    meta.seq_type = input[0]
    files = [ meta, input[1], input[2], input[3] ]
    return files
}

def map_join(channel_a, channel_b, key){
    channel_a
        .map{ it -> [it[key], it] }
        //.cross(channel_b.map{it -> [it[key], it]})
        //.map { it[0][1] + it[1][1] }
}

def join_by_st(ch_1, ch_2){ 
    // Use for cases where ch_1 = [ [meta.id, meta.seq_type] file(s)] and ch_2 = [ [meta.seq_type] file(s)]
    id = ch_1.map{ meta, file -> meta.id}
    new_1 = ch_1.map{ meta_old, file -> 
        def meta = [:]
        meta.seq_type = meta_old.seq_type
        return tuple( meta, file)} // Need to do this to keep [ [meta.seq_type], file] format. This just drops meta.id from channel
    new_1.view()
    ch_2.view()
    ch_3 = ch_1.map{ meta, file -> [ meta, file ]}.join(ch_2, by: [0,0])
    ch_3.view()
    ch_3.map{ meta_old, file -> 
        def meta = [:]
        meta.seq_type = meta_old.seq_type
        meta.id = id
        return tuple( meta, mash_dists)}
    ch_3.view()
    return ch_3
}

workflow SNVPHYL {
    take:
        samplesheet // GRIPHIN.out.griphin_samplesheet;l
        reference   // GET_REFERENCE_SEQ.out.reference

    main:
        ch_versions = Channel.empty() // Used to collect the software versions

        //This helps convert what was given to griphin to a sample sheet with the paths to fastq files
        CONVERT_INPUT (
            samplesheet
        )

        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        INPUT_CHECK (
            CONVERT_INPUT.out.updated_samplesheet
        )
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        // Run FastQC
        FASTQC (
            INPUT_CHECK.out.reads
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())

        //1. index process takes 1 input channel as a argument
        INDEXING(
            reference
        )

        //2. find repeats process takes 1 input channel as a argument
        FIND_REPEATS(
            reference
        )
        ch_versions = ch_versions.mix(FIND_REPEATS.out.versions)

        //3. smalt map process takes 2 input channels as arguments
        //INPUT_CHECK.out.reads.view()
        index_ch = INDEXING.out.ref_indexes.map{ add_st(it) }
        //index_ch.view()

        smalt_ch = join_by_st(INPUT_CHECK.out.reads, INDEXING.out.ref_indexes)
        //INDEXING.out.ref_indexes.view()
        //smalt_ch = INPUT_CHECK.out.reads.join(INDEXING.out.ref_indexes, by: [0,0])
        smalt_ch.view()
        //map_join(INPUT_CHECK.out.reads, index_ch, key: "seq_type" )

        // Function to get list of [ meta, [ scaffolds_1, scaffolds_2 ] ]
        //INDEXING.out.ref_indexes.map{ add_st(it) }.view()
        //INDEXING.out.ref_indexes.map{ seq_type, ref_fai, ref_sma, ref_smi -> [[seq_type], [ref_fai, ref_sma, ref_smi]]}.view()

        SMALT_MAP(
            INPUT_CHECK.out.reads, smalt_ch
        )
        ch_versions = ch_versions.mix(SMALT_MAP.out.versions)

        //4. sorting and indexing bam files from smalt process takes 1 input channel as an arguments
        SORT_INDEX_BAMS(
            SMALT_MAP.out.bams
        )
        ch_versions = ch_versions.mix(SORT_INDEX_BAMS.out.versions)

        //5. Generating mapping_quality.txt file
        GENERATE_LINE_1(
            SORT_INDEX_BAMS.out.sorted_bams.collect()
        )
        //ch_versions = ch_versions.mix(GENERATE_LINE_1.out.versions)

        VERIFYING_MAP_Q(
            SORT_INDEX_BAMS.out.sorted_bams.collect(), GENERATE_LINE_1.out.bam_lines_file.splitText()
        )
        ch_versions = ch_versions.mix(VERIFYING_MAP_Q.out.versions)

        //6. freebays variant calling process takes 2 input channels as arguments
        //SORT_INDEX_BAMS.out.sorted_bams_and_sampleID.view()
        //reference.view()
        FREEBAYES(
            SORT_INDEX_BAMS.out.sorted_bams_and_sampleID, reference
        )
        ch_versions = ch_versions.mix(FREEBAYES.out.versions)

        //7. filter freebays variant file process takes 1 input channel as an argument
        FILTER_FREEBAYES(
            FREEBAYES.out.vcf_files
        )
        ch_versions = ch_versions.mix(FILTER_FREEBAYES.out.versions)

        // Zip up the freebayes vcf
        BGZIP_FREEBAYES_VCF(
            FILTER_FREEBAYES.out.filtered_vcf
        )
        ch_versions = ch_versions.mix(BGZIP_FREEBAYES_VCF.out.versions)

        //8. Convert vcf freebays variant file to bcf process takes 1 input channel as an argument
        FREEBAYES_VCF_TO_BCF(
            BGZIP_FREEBAYES_VCF.out.filtered_zipped_vcf
        )
        ch_versions = ch_versions.mix(FREEBAYES_VCF_TO_BCF.out.versions)

        //9. mplileup process takes 1 input channel as argument
        MPILEUP(
            SORT_INDEX_BAMS.out.sorted_bams_and_sampleID, reference
        )
        ch_versions = ch_versions.mix(MPILEUP.out.versions)

        // Zip up the mpileup vcf
        BGZIP_MPILEUP_VCF(
            MPILEUP.out.mpileup
        )
        ch_versions = ch_versions.mix(BGZIP_MPILEUP_VCF.out.versions)

        //10. mplileup variant calls takes 1 input channel as an argument
        BCFTOOLS_CALL(
            BGZIP_MPILEUP_VCF.out.mpileup_zipped
        )
        ch_versions = ch_versions.mix(BCFTOOLS_CALL.out.versions)

        //Joining channels of multiple outputs
        combined_ch = BCFTOOLS_CALL.out.mpileup_bcf.join(FREEBAYES_VCF_TO_BCF.out.filtered_bcf)
        //11. consolidate variant calling files process takes 2 input channels as arguments
        CONSOLIDATE_BCFS(
            combined_ch
        )
        ch_versions = ch_versions.mix(CONSOLIDATE_BCFS.out.versions)

        // Concat filtered densities to make new invalid_postions
        CONSOLIDATE_FILTERED_DENSITY(
            CONSOLIDATE_BCFS.out.filtered_densities.collect(), FIND_REPEATS.out.repeats_bed_file
        )
        ch_versions = ch_versions.mix(CONSOLIDATE_FILTERED_DENSITY.out.versions)

        // Making string that looks like... this is needed for the next process
        //--consolidate_vcf 2021JQ-00457-WAPHL-M5130-211029=2021JQ-00457-WAPHL-M5130-211029_consolidated.bcf --consolidate_vcf 2021JQ-00459-WAPHL-M5130-211029=2021JQ-00459-WAPHL-M5130-211029_consolidated.bcf --consolidate_vcf 2021JQ-00460-WAPHL-M5130-211029=2021JQ-00460-WAPHL-M5130-211029_consolidated.bcf
        GENERATE_LINE_2(
            CONSOLIDATE_BCFS.out.consolidated_bcfs.collect()
        )
        //ch_versions = ch_versions.mix(GENERATE_LINE_2.out.versions)

        // Get line out of file we just made that has the --consolidate_vcf line...
        //13. consolidate variant calling files process takes 2 input channels as arguments
        VCF2SNV_ALIGNMENT(
            GENERATE_LINE_2.out.consolidation_line.splitText(), CONSOLIDATE_BCFS.out.consolidated_bcfs.collect(), CONSOLIDATE_FILTERED_DENSITY.out.new_invalid_positions, reference, CONSOLIDATE_BCFS.out.consolidated_bcf_index.collect()
        )
        ch_versions = ch_versions.mix(VCF2SNV_ALIGNMENT.out.versions)

        //14. Filter Stats
        FILTER_STATS(
            VCF2SNV_ALIGNMENT.out.snvTable
        )
        ch_versions = ch_versions.mix(FILTER_STATS.out.versions)

        //15. Using phyml to build tree process takes 1 input channel as an argument
        PHYML(
            VCF2SNV_ALIGNMENT.out.snvAlignment
        )
        ch_versions = ch_versions.mix(PHYML.out.versions)

        //16. Make SNVMatix.tsv
        MAKE_SNV(
            VCF2SNV_ALIGNMENT.out.snvAlignment
        )
        ch_versions = ch_versions.mix(MAKE_SNV.out.versions)

    emit:
        versions = ch_versions // channel: [ versions.yml ]


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
