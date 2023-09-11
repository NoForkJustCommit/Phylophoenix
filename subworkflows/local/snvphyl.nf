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

def join_index_by_st(reads_ch, reference_index_list){ 
    // Use for cases where ch_1 = [ [meta.id, meta.seq_type] file(s)] and ch_2 = [ [meta.seq_type] file(s)]
    seq_type_1 = reads_ch.get(0).seq_type // get the meta.seq_type for the reads
    // now loop through reference list to find a match
    for ( ref_indexes in reference_index_list ) {
        // get reference_indexes seq_type
        seq_type_2 = ref_indexes.get(0).seq_type // get the meta.seq_type for the reference indexes
        if ( seq_type_2 == seq_type_1 ) {// if the seq_type match create a new channel
            reads = reads_ch.get(1)
            matched_ch = [ [reads_ch.get(0), reads], [ref_indexes.get(0), ref_indexes.get(1), ref_indexes.get(2), ref_indexes.get(3) ]]
        } 
    }
    return matched_ch
}

def join_ref_by_st(bams_ch, reference_index_list, ref_genome_list){ 
    // Use for cases where ch_1 = [ [meta.id, meta.seq_type] file(s)] and ch_2 = [ [meta.seq_type] file(s)]
    seq_type_1 = bams_ch.get(0).seq_type // get the meta.seq_type for the reads
    // now loop through reference list to find a match
    for ( ref_indexes in reference_index_list ) {
        // get reference seq_type
        seq_type_2 = ref_indexes.get(0).seq_type // get the meta.seq_type for the reference indexes
        if ( seq_type_2 == seq_type_1 ) {// if the seq_type match create a new channel
            for ( ref_genome in ref_genome_list ) {
                seq_type_3 = ref_genome.get(0).seq_type // get the meta.seq_type for the reference genome
                if (seq_type_2 == seq_type_3) {
                    matched_ch = [ [bams_ch.get(0), bams_ch.get(1)], [ref_indexes.get(0), ref_indexes.get(1), ref_indexes.get(2), ref_indexes.get(3) ], [ref_genome.get(0), ref_genome.get(1)]]
                }
            }
        } 
    }
    return matched_ch
}

workflow SNVPHYL {
    take:
        samplesheet // GRIPHIN.out.griphin_samplesheet
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
        ch_versions = ch_versions.mix(INDEXING.out.versions)

        //2. find repeats process takes 1 input channel as a argument
        FIND_REPEATS(
            reference
        )
        ch_versions = ch_versions.mix(FIND_REPEATS.out.versions)

        //Doing some channel rearragement to be able to combine the correct reference indexes for the st with the reads of the same st
        smalt_ch = INPUT_CHECK.out.reads.map{meta, reads -> [[meta, reads]] } // just reformating to use in function
            .combine(INDEXING.out.ref_indexes.map{meta, fai, sma, smi -> [[meta, fai, sma, smi]] }.collect())// get all references to loop through
            .map{ reads, ref_index -> join_index_by_st(reads, [ref_index]) } // use custom function to combine by seq_type

        //3. smalt map process takes 2 input channels as arguments
        SMALT_MAP(
            smalt_ch.map{meta_and_reads, meta_2_and_ref_indexes -> meta_and_reads },
            smalt_ch.map{meta_and_reads, meta_2_and_ref_indexes -> meta_2_and_ref_indexes }
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
        ch_versions = ch_versions.mix(GENERATE_LINE_1.out.versions)

        VERIFYING_MAP_Q(
            SORT_INDEX_BAMS.out.sorted_bams.collect(), GENERATE_LINE_1.out.bam_lines_file.splitText()
        )
        ch_versions = ch_versions.mix(VERIFYING_MAP_Q.out.versions)

        //Doing some channel rearragement to be able to combine the correct reference and its indexes for the st with the bams from the same st
        sorted_bams_with_ch = SORT_INDEX_BAMS.out.sorted_bams_and_sampleID.map{meta, bams -> [[meta, bams]] } // just reformating to use in function
            .combine(INDEXING.out.ref_indexes.map{meta, fai, sma, smi -> [[meta, fai, sma, smi]] }.collect())// get all references to loop through
            .combine(reference.map{meta, reference -> [[meta, reference]] }.collect())// get all references index files to loop through
            .map{ bams, ref_indexes, reference -> join_ref_by_st(bams, [ref_indexes], [reference]) } // use custom function to combine by seq_type

        //6. freebays variant calling process takes 3 input channels as arguments. You can do without indexes, but then freebayes makes them so including
        FREEBAYES(
            sorted_bams_with_ch.map{meta_and_bams, meta_2_and_ref_indexes, meta_3_and_reference -> meta_and_bams },
            sorted_bams_with_ch.map{meta_and_bams, meta_2_and_ref_indexes, meta_3_and_reference -> meta_2_and_ref_indexes },
            sorted_bams_with_ch.map{meta_and_bams, meta_2_and_ref_indexes, meta_3_and_reference -> meta_3_and_reference }
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

        //9. mplileup process takes 1 input channel as argument. You can do without indexes, but then mpileup makes them so including for speed
        MPILEUP(
            sorted_bams_with_ch.map{meta_and_bams, meta_2_and_ref_indexes, meta_3_and_reference -> meta_and_bams },
            sorted_bams_with_ch.map{meta_and_bams, meta_2_and_ref_indexes, meta_3_and_reference -> meta_2_and_ref_indexes },
            sorted_bams_with_ch.map{meta_and_bams, meta_2_and_ref_indexes, meta_3_and_reference -> meta_3_and_reference }
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
