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
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/

def add_st(input) {
    // create meta map
    def meta = [:]
    meta.seq_type = input[0]
    files = [ meta, input[1], input[2], input[3] ]
    return files
}

def collect_and_format_indexes(input_ch){
    //If you use the .collect() alone you get --> [seq_type:ST398], 2023CB-00249.filtered.scaffolds.fa.fai, 2023CB-00249.filtered.scaffolds.fa.sma, 2023CB-00249.filtered.scaffolds.fa.smi
    //This function returns the same info like this --> [[seq_type:ST398], 2023CB-00249.filtered.scaffolds.fa.fai, 2023CB-00249.filtered.scaffolds.fa.sma, 2023CB-00249.filtered.scaffolds.fa.smi]
    //This allows you to keep the index information together with the ST type
    def sample_num = input_ch.size() - 1 // substract 1 to handle zero indexing
    ref_index_list = []
    count = 0 // inti count
    while (count <= sample_num) { // loop through each index (one per st)
        meta = input_ch.get(count)
        fai = input_ch.get(count+1)
        sma = input_ch.get(count+2)
        smi = input_ch.get(count+3)
        ref_index_list.add([meta, fai, sma, smi])
        count=count+4
    }
    return ref_index_list
}

def collect_and_format_refs(input_ch){
    //If you use the .collect() alone you get --> [seq_type:ST398], 2023CB-00249.filtered.scaffolds.fa.fai, 2023CB-00249.filtered.scaffolds.fa.sma, 2023CB-00249.filtered.scaffolds.fa.smi
    //This function returns the same info like this --> [[seq_type:ST398], 2023CB-00249.filtered.scaffolds.fa.fai, 2023CB-00249.filtered.scaffolds.fa.sma, 2023CB-00249.filtered.scaffolds.fa.smi]
    //This allows you to keep the index information together with the ST type
    def sample_num = input_ch.size() - 1 // substract 1 to handle zero indexing
    ref_seqs_list = []
    count = 0
    while (count <=sample_num) { // loop through each index (one per st)
        meta = input_ch.get(count)
        fasta = input_ch.get(count+1)
        ref_seqs_list.add([meta, fasta])
        count=count+2
    }
    return ref_seqs_list
}

def join_index_by_st(input_ch){
    // Use for cases where ch_1 = [ [meta.id, meta.seq_type] file(s)] and ch_2 = [ [meta.seq_type] file(s) ]
    reads_ch = input_ch.take(1).get(0) // getting meta information from the reads (lists inside of lists...)
    seq_type_1 = reads_ch.get(0).seq_type // get the meta.seq_type for the reads
    //size of list
    list_size = input_ch.size() - 1 // have to add one because of indexing 0 issues\
    // now loop through reference list to find a match
    for ( int i in (1..list_size) ) {
        ref_indexes = input_ch.get(i) // remove is actually pulling that one selected
        // get reference_indexes seq_type
        seq_type_2 = ref_indexes.get(0).seq_type // get the meta.seq_type for the reference indexes
        if ( seq_type_2 == seq_type_1 ) {// if the seq_type match create a new channel
            reads = reads_ch.get(1)
            matched_ch = [ [reads_ch.get(0), reads_ch.get(1)], [ref_indexes.get(0), ref_indexes.get(1), ref_indexes.get(2), ref_indexes.get(3) ] ]
        }
    }
    return matched_ch
}

def join_ref_by_st(input_ch){
    // Use for cases where ch_1 = [ [meta.id, meta.seq_type] file(s)] and ch_2 = [ [meta.seq_type] file(s)]
    bams_ch = input_ch.get(0) // getting meta information from the bams (lists inside of lists...)
    seq_type_1 = bams_ch.get(0).seq_type // get the meta.seq_type for the bams
    ref_indexes = input_ch.get(1) // getting meta information from the indexes (lists inside of lists...)
    seq_type_2 = ref_indexes.get(0).seq_type // get the meta.seq_type for the indexes
    // get size of list
    list_size = input_ch.size() - 1 // have to add one because of indexing 0 issues
    count = 0 // inti count
    // now loop through reference list to find a match
    for ( int i in (2..list_size) ) {
        ref_genome = input_ch.get(i) // get the reference
        // get reference seq_type
        seq_type_3 = ref_genome.get(0).seq_type // get the meta.seq_type for the reference sequence
        if ( seq_type_3 == seq_type_1 ) {// if the seq_type of reference sequence and the bams match create a new channel
            // double check the reference seq type matches the ref_indexes seq type
            if (seq_type_3 == seq_type_2) {
                matched_ch = [ [bams_ch.get(0), bams_ch.get(1)], [ref_indexes.get(0), ref_indexes.get(1), ref_indexes.get(2), ref_indexes.get(3)], [ref_genome.get(0), ref_genome.get(1)] ]
            } else {
                exit 1, "ERROR: the seq type of the refernce indexes does not match. This is a bug please open github issue!\n"
            }
        }
    }
    return matched_ch
}

def sort_collected_by_st(input_ch){
    // this function takes in a list of one file type with each one having a meta.seq_type and returns:
    // [ [ meta.seq_type ], file_1, file_2, file_3 ]
    count = 0
    sample_num = 1 // this will be used to get bam that matches the
    def st_list = []
    for (isolate in input_ch) { // loop through each st type
        sequence_type =isolate.get(0).seq_type // get st
        st_list.add(sequence_type.toString()) // add st to list
    }
    //get the unique st's
    def unique_st = st_list.unique()

    def complete_list = []
    for ( looping_st in unique_st) { // loop through each st
        def isolate_st_list = [] // create an empty list to store bams that are of the same ST
        def meta = [:]
        meta.seq_type = looping_st // set current looping st
        for (isolate in input_ch) { // loop through isolates
            checking_st = isolate.get(0).seq_type // get st of current isolate
            if ( looping_st == checking_st ) { // if the current st is the same as the looping st add isolate to list
                // add it to the tuple that was created before
                isolate_st_list.add(isolate.get(1)) // add bam file to the list
            }
        }
        complete_list.add([meta, isolate_st_list])
    }
    return complete_list.get(0)
}

def sort_collected_by_st_b(input_ch, input_st){
    def complete_list = []
    def isolate_st_list = [] // create an empty list to store bams that are of the same ST
    def meta = [:]
    meta.seq_type = input_st // set current looping st
    for (isolate in input_ch) { // loop through isolates
        checking_st = isolate.get(0).seq_type // get st of current isolate
        if ( checking_st == input_st ) { // if the current st is the same as the looping st add isolate to list
            // add it to the tuple that was created before
            isolate_st_list.add(isolate.get(1)) // add bam file to the list
        }
    }
    complete_list.add([meta, isolate_st_list])
    return complete_list.get(0)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

        // reformat to get all index files in the same list together --> helps with joining downstream
        index_ch = INDEXING.out.ref_indexes.collect().map{ collect_and_format_indexes(it) }
        // reformat to get reference sequence file and meta.seq_type in the same list together --> helps with joining downstream
        ref_ch = reference.collect().map{ collect_and_format_refs(it) }
        // get the st 

        //2. find repeats process takes 1 input channel as a argument
        FIND_REPEATS(
            reference
        )
        ch_versions = ch_versions.mix(FIND_REPEATS.out.versions)

        //Doing some channel rearragement to be able to combine the correct reference indexes for the st with the reads of the same st
        smalt_ch = INPUT_CHECK.out.reads.map{meta, reads -> [[meta, reads]] } // just reformating to use in function
            .combine(index_ch) // add in all the index files
            .map{ join_index_by_st(it) } // use custom function to combine by seq_type

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

        // collect all sorted bam files and separate them in their own channel by ST
        // We have to start with the list of STs and flatten or so we end up with the right number of output channels. 
        sort_indexed_bams_ch = INPUT_CHECK.out.st_list.flatten().combine( // take the list of STs in data set and flatten so one goes in at a time. 
            SORT_INDEX_BAMS.out.sorted_bams.map{ meta, bams -> [[meta, bams]]}// use map to put the meta and bam in same list so they are coupled
            .collect().map{ it -> [it] } // collect all and add [] so its all one channel for the next map. 
        ).map{ sts, meta_and_bams -> sort_collected_by_st_b(meta_and_bams, sts) } // get sorted bams afor the st

        //5. Generating mapping_quality.txt file
        GENERATE_LINE_1(
            sort_indexed_bams_ch
        )
        ch_versions = ch_versions.mix(GENERATE_LINE_1.out.versions)

        // Combine sorted bams and bam_lines by ST
        verifying_map_q_ch = sort_indexed_bams_ch.join(GENERATE_LINE_1.out.bam_lines_file, by: [0])

        VERIFYING_MAP_Q(
            verifying_map_q_ch.map{ meta, sorted_bams, bam_lines_file -> [ meta, sorted_bams ]},
            verifying_map_q_ch.map{ meta, sorted_bams, bam_lines_file -> bam_lines_file}.splitText()
        )
        ch_versions = ch_versions.mix(VERIFYING_MAP_Q.out.versions)

        // Doing some channel rearragement to be able to combine the correct reference and its indexes for the st with the bams from the same st
        sorted_bams_with_ch = SORT_INDEX_BAMS.out.sorted_bams.map{meta, bams -> [[meta, bams]] } // just reformating to use in function
            .combine(index_ch) // add in all the index files
            .map{ join_index_by_st(it) } // use custom function to combine indexes with bams by seq_type
            .combine(ref_ch) // add in all the reference files
            .map{ join_ref_by_st(it) }// use custom function to combine by seq_type

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
        combined_ch = BCFTOOLS_CALL.out.mpileup_bcf.join(FREEBAYES_VCF_TO_BCF.out.filtered_bcf, by: [0])

        //11. consolidate variant calling files process takes 2 input channels as arguments
        CONSOLIDATE_BCFS(
            combined_ch
        )
        ch_versions = ch_versions.mix(CONSOLIDATE_BCFS.out.versions)

        // collect all sorted filter densities files and separate them in their own channel by ST
        // We have to start with the list of STs and flatten or so we end up with the right number of output channels. 
        consolidate_filtered_densities_ch = INPUT_CHECK.out.st_list.flatten().combine( // take the list of STs in data set and flatten so one goes in at a time. 
            CONSOLIDATE_BCFS.out.filtered_densities.map{ meta, filtered_densities -> [[ meta, filtered_densities ]]} // use map to put the meta and filter densities in same list so they are coupled
            .collect().map{ it -> [it] } // collect all and add [] so its all one channel for the next map. 
        ).map{ sts, meta_and_filtered_densities -> sort_collected_by_st_b(meta_and_filtered_densities, sts) } // get sorted filter densities afor the st
        .join(FIND_REPEATS.out.repeats_bed_file, by: [0]) // add in the bed file for the st

        // Concat filtered densities to make new invalid_postions
        CONSOLIDATE_FILTERED_DENSITY(
            consolidate_filtered_densities_ch.map{ meta, filtered_densities, repeats_bed_file -> [ meta, filtered_densities ] },
            consolidate_filtered_densities_ch.map{ meta, filtered_densities, repeats_bed_file -> [ meta, repeats_bed_file ] }
        )
        ch_versions = ch_versions.mix(CONSOLIDATE_FILTERED_DENSITY.out.versions)

        // collect all bcf files and separate them in their own channel by ST
        // We have to start with the list of STs and flatten or so we end up with the right number of output channels. 
        consolidated_bcfs_ch = INPUT_CHECK.out.st_list.flatten().combine( // take the list of STs in data set and flatten so one goes in at a time. 
            CONSOLIDATE_BCFS.out.consolidated_bcfs.map{ meta, consolidated_bcfs -> [[ meta, consolidated_bcfs ]]} // use map to put the meta and bcfs in same list so they are coupled
            .collect().map{ it -> [it] } // collect all and add [] so its all one channel for the next map. 
        ).map{ sts, meta_and_consolidated_bcfs -> sort_collected_by_st_b(meta_and_consolidated_bcfs, sts) } // get sorted bcfs for the st

        // Making string that looks like... this is needed for the next process
        //--consolidate_vcf 2021JQ-00457-WAPHL-M5130-211029=2021JQ-00457-WAPHL-M5130-211029_consolidated.bcf --consolidate_vcf 2021JQ-00459-WAPHL-M5130-211029=2021JQ-00459-WAPHL-M5130-211029_consolidated.bcf --consolidate_vcf 2021JQ-00460-WAPHL-M5130-211029=2021JQ-00460-WAPHL-M5130-211029_consolidated.bcf
        GENERATE_LINE_2(
            consolidated_bcfs_ch
        )
        ch_versions = ch_versions.mix(GENERATE_LINE_2.out.versions)

        // collect all sorted bcf index files and separate them into their own channel by ST
        consolidate_bcf_indexes_ch = INPUT_CHECK.out.st_list.flatten().combine( // take the list of STs in data set and flatten so one goes in at a time. 
            CONSOLIDATE_BCFS.out.consolidated_bcf_index.map{ meta, consolidated_bcf_index -> [[ meta, consolidated_bcf_index ]]}
            .collect().map{ it -> [it] } // collect all and add [] so its all one channel for the next map. 
        ).map{ sts, meta_and_consolidated_bcf_index -> sort_collected_by_st_b(meta_and_consolidated_bcf_index, sts) } // get sorted bcfs for the st

        vcf2snv_alignment_ch = GENERATE_LINE_2.out.consolidation_line.join(consolidated_bcfs_ch, by: [0])
                                .join(CONSOLIDATE_FILTERED_DENSITY.out.new_invalid_positions, by: [0])
                                .join(reference, by: [0])
                                .join(consolidate_bcf_indexes_ch, by: [0])

        // Get line out of file we just made that has the --consolidate_vcf line...
        //13. consolidate variant calling files process takes 2 input channels as arguments
        VCF2SNV_ALIGNMENT(
            vcf2snv_alignment_ch.map{ meta, consolidation_line, consolidated_bcfs, new_invalid_positions, reference, consolidate_bcf_indexes -> [ meta, consolidation_line ] }.splitText(),
            vcf2snv_alignment_ch.map{ meta, consolidation_line, consolidated_bcfs, new_invalid_positions, reference, consolidate_bcf_indexes -> [ meta, consolidated_bcfs ] },
            vcf2snv_alignment_ch.map{ meta, consolidation_line, consolidated_bcfs, new_invalid_positions, reference, consolidate_bcf_indexes -> [ meta, new_invalid_positions ] },
            vcf2snv_alignment_ch.map{ meta, consolidation_line, consolidated_bcfs, new_invalid_positions, reference, consolidate_bcf_indexes -> [ meta, reference ] },
            vcf2snv_alignment_ch.map{ meta, consolidation_line, consolidated_bcfs, new_invalid_positions, reference, consolidate_bcf_indexes -> [ meta, consolidate_bcf_indexes ] }
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
