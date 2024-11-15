/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPhylophoenix.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Modules for running all vs all
//
include { GRIPHIN_WORKFLOW             } from '../subworkflows/local/griphin_workflow'
include { GET_COMPARISONS              } from '../modules/local/create_comparisons'
include { MASH_DIST    as MASH_DIST    } from '../modules/local/mash_distance'
include { GET_CENTROID as GET_CENTROID } from '../modules/local/get_centroid'
include { ASSET_PREP   as ASSET_PREP   } from '../modules/local/asset_prep'
include { CREATE_META  as CREATE_META  } from '../subworkflows/local/create_meta'
include { SNVPHYL                      } from '../subworkflows/local/snvphyl'
include { RENAME_REF_IN_OUTPUT         } from '../modules/local/rename_ref_in_output'
include { CLEAN_AND_CREATE_METADATA    } from '../modules/local/clean_and_create_metadata'
include { COMBINE_GRIPHIN_SNVPHYL      } from '../modules/local/combine_griphin_snvphyl'

//
// Modules for running snvphyl by st
//
include { GET_SEQUENCE_TYPES                                           } from '../modules/local/get_sequence_types'
include { CREATE_META               as CREATE_META_BY_ST               } from '../subworkflows/local/create_meta'
include { MASH_DIST                 as MASH_DIST_BY_ST                 } from '../modules/local/mash_distance'
include { GET_CENTROID              as GET_CENTROID_BY_ST              } from '../modules/local/get_centroid'
include { ASSET_PREP                as ASSET_PREP_BY_ST                } from '../modules/local/asset_prep'
include { SNVPHYL                   as SNVPHYL_BY_ST                   } from '../subworkflows/local/snvphyl'
include { RENAME_REF_IN_OUTPUT      as RENAME_REF_IN_OUTPUT_BY_ST      } from '../modules/local/rename_ref_in_output'
include { CLEAN_AND_CREATE_METADATA as CLEAN_AND_CREATE_METADATA_BY_ST } from '../modules/local/clean_and_create_metadata'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
//def multiqc_report = []

workflow PHYLOPHOENIX {
    take:
        input_samplesheet_path
        indir

    main:
        ch_versions = Channel.empty()

        // get geonames file into a channel. Coding it this way and not params so we can use a glob and not be verbose.
        geonames_ch = Channel.fromPath("${baseDir}/assets/databases/*_geolocation.txt.xz").collect()

        // Create report
        GRIPHIN_WORKFLOW (
            input_samplesheet_path, indir
        )
        ch_versions = ch_versions.mix(GRIPHIN_WORKFLOW.out.versions)

        // If you pass --no_all then samples will not be run all together
        if (params.no_all==false) {
            // Allow outdir to be relative
            //outdir_path = Channel.fromPath(params.outdir, relative: true, type: 'dir')

            // Creates samplesheets with sampleid,seq_type,path_to_assembly
            GET_COMPARISONS (
                GRIPHIN_WORKFLOW.out.directory_samplesheet
            )
            ch_versions = ch_versions.mix(GET_COMPARISONS.out.versions)

            if (params.metadata!=null) {
                // get metadata into channel
                metadata  = Channel.fromPath(params.metadata, relative: true)
                // creates channel: [ val(meta.id, meta.st), [ scaffolds_1, scaffolds_2 ] ]
                CREATE_META (
                    GET_COMPARISONS.out.samplesheet, GET_COMPARISONS.out.snv_samplesheet, metadata, false
                )
            } else {
                // creates channel: [ val(meta.id, meta.st), [ scaffolds_1, scaffolds_2 ] ]
                CREATE_META (
                    GET_COMPARISONS.out.samplesheet, GET_COMPARISONS.out.snv_samplesheet, null, false
                )
            }

            // Run Mash on groups of samples by seq type
            MASH_DIST (
                CREATE_META.out.st_scaffolds
            )
            ch_versions = ch_versions.mix(MASH_DIST.out.versions)

            // Creating channel [ ST, [distance_1, distance_2] ]
            dist_ch = MASH_DIST.out.dist.map{ meta, dist -> [ dist ] }.collect() // drop meta and collect all distance files
            centroid_ch = dist_ch.map{ mash_dist -> 
                        def meta = [:]
                        meta.seq_type = "All_STs"
                        return tuple (meta , mash_dist)} // add back "All_STs" as the meta value
            .combine(GRIPHIN_WORKFLOW.out.directory_samplesheet) // Add samplesheet to all mash distance channels

            // Take in all mash distance files then use the samplesheet to return the centroid assembly
            // Get centroid, by calculating the average mash distance
            GET_CENTROID (
                centroid_ch
            )
            ch_versions = ch_versions.mix(GET_CENTROID.out.versions)

            // Unzip centroid assembly: SNVPhyl requires it unzipped
            // also unzip the geoname files for cleaning metadata file. 
            ASSET_PREP (
                // Bring in centroid into channel
                GET_CENTROID.out.centroid_path.splitCsv( header:false, sep:',' ).map{meta, list -> 
                def scaffold = list[0] // extract the file from the list
                return [meta, scaffold]},  // get into format [[meta], scaffold]
                geonames_ch, CREATE_META.out.st_snv_samplesheets
            )
            ch_versions = ch_versions.mix(ASSET_PREP.out.versions)

            // Check and correct the metadata file if it was passed
            if (params.metadata!=null) {
                CLEAN_AND_CREATE_METADATA (
                    CREATE_META.out.split_metadata, GRIPHIN_WORKFLOW.out.griphin_tsv_report, ASSET_PREP.out.unzipped_geodata
                )
                ch_versions = ch_versions.mix(CLEAN_AND_CREATE_METADATA.out.versions)
            }

            // Make SNVPHYL channel by joining by seq type
            all_ch = ASSET_PREP.out.st_snv_samplesheets.join(ASSET_PREP.out.unzipped_fasta, by: [0])

            // Run snvphyl on each st type on its own input
            SNVPHYL (
                all_ch.map{ seq_type, samplesheet, unzipped_fasta -> [seq_type, samplesheet]}, all_ch.map{ seq_type, samplesheet, unzipped_fasta -> [seq_type, unzipped_fasta] } // reference
            )
            ch_versions = ch_versions.mix(SNVPHYL.out.versions)

            if (params.metadata!=null) {
                final_output_ch = GET_CENTROID.out.centroid_info.join(SNVPHYL.out.phylogeneticTree, by: [0]).join(SNVPHYL.out.snvMatrix, by: [0]).join(CLEAN_AND_CREATE_METADATA.out.updated_samplesheet, by: [0])
            } else {
                final_output_ch = GET_CENTROID.out.centroid_info.join(SNVPHYL.out.phylogeneticTree, by: [0]).join(SNVPHYL.out.snvMatrix, by: [0]).combine([])
            }

            // Rename reference to actual sample name
            RENAME_REF_IN_OUTPUT (
                final_output_ch
            )
            ch_versions = ch_versions.mix(RENAME_REF_IN_OUTPUT.out.versions)
        }

        // If you pass --by_st then samples will be broken up by st type and SNVPhyl run on each st on its own
        if (params.by_st==true) {

            // Creates samplesheets with sample,seq_type,path_to_assembly
            GET_SEQUENCE_TYPES (
                GRIPHIN_WORKFLOW.out.directory_samplesheet, GRIPHIN_WORKFLOW.out.griphin_report
            )
            ch_versions = ch_versions.mix(GET_SEQUENCE_TYPES.out.versions)

            if (params.metadata!=null) {
                // get metadata into channel
                metadata  = Channel.fromPath(params.metadata, relative: true)

                // creates channel: [ val(meta.id, meta.st), [ scaffolds_1, scaffolds_2 ] ]
                CREATE_META_BY_ST (
                    GET_SEQUENCE_TYPES.out.st_samplesheets, GET_SEQUENCE_TYPES.out.st_snv_samplesheets, metadata, true
                )
            } else {
                // creates channel: [ val(meta.id, meta.st), [ scaffolds_1, scaffolds_2 ] ]
                CREATE_META_BY_ST (
                    GET_SEQUENCE_TYPES.out.st_samplesheets, GET_SEQUENCE_TYPES.out.st_snv_samplesheets, [], false
                )
            }

            // Run Mash on groups of samples by seq type
            MASH_DIST_BY_ST (
                CREATE_META_BY_ST.out.st_scaffolds
            )
            ch_versions = ch_versions.mix(MASH_DIST_BY_ST.out.versions)

            // Creating channel [ ST, [distance_1, distance_2] ]
            st_mash_dists = MASH_DIST_BY_ST.out.dist.map{ meta, mash_dist -> 
                    def key = meta.seq_type 
                    return tuple (key, mash_dist)}
                .groupTuple() //group by st type. this returns [ST#, [ distance_1, distance_2 ]]
                .map{meta_old, mash_dists -> 
                    def meta = [:]
                    meta.seq_type = meta_old
                    return tuple( meta, mash_dists)} // Now returns [[meta.seq_type], [ distance_1, distance_2 ]]

            // Add samplesheet to all mash distance channels
            centroid_st_ch =  st_mash_dists.combine(GRIPHIN_WORKFLOW.out.directory_samplesheet)

            // Take in all mash distance files for a seq type and then use the samplesheet to return the centroid assembly
            // Get centroid for ST groups, by calculating the average mash distance
            GET_CENTROID_BY_ST (
                centroid_st_ch
            )
            ch_versions = ch_versions.mix(GET_CENTROID_BY_ST.out.versions)

            // Unzip centroid assembly: SNVPhyl requires it unzipped
            ASSET_PREP_BY_ST (
                GET_CENTROID_BY_ST.out.centroid_path.splitCsv( header:false, sep:',' ), // Bring in centroid into channel
                geonames_ch,
                CREATE_META_BY_ST.out.st_snv_samplesheets
            )
            ch_versions = ch_versions.mix(ASSET_PREP_BY_ST.out.versions)

            // Check and correct the metadata file
            if (params.metadata!=null) {
                // clean up metadata file, add geolocation information
                CLEAN_AND_CREATE_METADATA_BY_ST (
                    CREATE_META_BY_ST.out.split_metadata, GRIPHIN_WORKFLOW.out.griphin_tsv_report, ASSET_PREP.out.unzipped_geodata
                )
                ch_versions = ch_versions.mix(CLEAN_AND_CREATE_METADATA_BY_ST.out.versions)
            }

            // Make SNVPHYL channel by joining by seq type
            st_ch = CREATE_META_BY_ST.out.st_snv_samplesheets.join(ASSET_PREP_BY_ST.out.unzipped_fasta, by: [0])

            // Run snvphyl on each st type on its own input
            SNVPHYL_BY_ST (
                st_ch.map{ seq_type, samplesheet, unzipped_fasta -> [seq_type, samplesheet]}, 
                st_ch.map{ seq_type, samplesheet, unzipped_fasta -> [seq_type, unzipped_fasta] } // reference
            )
            ch_versions = ch_versions.mix(SNVPHYL_BY_ST.out.versions)

            if (params.metadata!=null) {
                final_st_output_by_st_ch = GET_CENTROID_BY_ST.out.centroid_info.join(SNVPHYL_BY_ST.out.phylogeneticTree, by: [0]).join(SNVPHYL_BY_ST.out.snvMatrix, by: [0]).join(CLEAN_AND_CREATE_METADATA_BY_ST.out.updated_samplesheet, by: [0])
            } else {
                final_st_output_by_st_ch = GET_CENTROID_BY_ST.out.centroid_info.join(SNVPHYL_BY_ST.out.phylogeneticTree, by: [0]).join(SNVPHYL_BY_ST.out.snvMatrix, by: [0]).combine([])
            }

            // Rename reference to actual sample name
            RENAME_REF_IN_OUTPUT_BY_ST (
                final_st_output_by_st_ch
            )
            ch_versions = ch_versions.mix(RENAME_REF_IN_OUTPUT_BY_ST.out.versions)
        }

        if (params.by_st==true) {
            // collect files to add to griphin summary
            snvMatrix_ch = RENAME_REF_IN_OUTPUT.out.snvMatrix.collect().combine(RENAME_REF_IN_OUTPUT_BY_ST.out.snvMatrix.collect())
            vcf2core_ch = SNVPHYL.out.vcf2core.map{ meta, vcf2core -> vcf2core }.collect().combine(SNVPHYL_BY_ST.out.vcf2core.map{ meta, vcf2core -> vcf2core }.collect())
        } else {
            // collect files to add to griphin summary
            snvMatrix_ch = RENAME_REF_IN_OUTPUT.out.snvMatrix
            vcf2core_ch = SNVPHYL.out.vcf2core.map{ meta, vcf2core -> vcf2core }
        }

        COMBINE_GRIPHIN_SNVPHYL (
            snvMatrix_ch, vcf2core_ch, GRIPHIN_WORKFLOW.out.griphin_report
        )
        ch_versions = ch_versions.mix(COMBINE_GRIPHIN_SNVPHYL.out.versions)

        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

    emit:
        griphin_report = GRIPHIN_WORKFLOW.out.griphin_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
