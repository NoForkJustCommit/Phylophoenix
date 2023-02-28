/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPhylophoenix.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
//if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check input path parameters to see if they exist
if (params.input != null ) {  // if a samplesheet is passed
    // allow input to be relative, turn into string and strip off the everything after the last backslash to have remainder of as the full path to the samplesheet. 
    //input_samplesheet_path = Channel.fromPath(params.input, relative: true).map{ [it.toString().replaceAll(/([^\/]+$)/, "").replaceAll(/\/$/, "") ] }
    input_samplesheet_path = Channel.fromPath(params.input, relative: true)
    if (params.input_dir != null ) { //if samplesheet is passed and an input directory exit
        exit 1, 'You need EITHER an input samplesheet or a directory! Just pick one.' 
    } else { // if only samplesheet is passed check to make sure input is an actual file
        def checkPathParamList = [ params.input ]
        for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
    }
} else {
    if (params.input_dir != null ) { // if no samplesheet is passed, but an input directory is given
        input_samplesheet_path = []
        def checkPathParamList = [ params.input_dir ]
        for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
    } else { // if no samplesheet is passed and no input directory is given
        exit 1, 'You need EITHER an input samplesheet or a directory!' 
    }
}

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
include { ASSET_CHECK  as ASSET_CHECK  } from '../modules/local/asset_check'
include { CREATE_META  as CREATE_META  } from '../subworkflows/local/create_meta'
include { SNVPHYL                      } from '../subworkflows/local/snvphyl'

//
// Modules for running snvphyl by st
//
include { GET_SEQUENCE_TYPES                 } from '../modules/local/get_sequence_types'
include { CREATE_META  as CREATE_META_BY_ST  } from '../subworkflows/local/create_meta'
include { MASH_DIST    as MASH_DIST_BY_ST    } from '../modules/local/mash_distance'
include { GET_CENTROID as GET_CENTROID_BY_ST } from '../modules/local/get_centroid'
include { ASSET_CHECK  as ASSET_CHECK_BY_ST  } from '../modules/local/asset_check'
include { SNVPHYL      as SNVPHYL_BY_ST      } from '../subworkflows/local/snvphyl'

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
def multiqc_report = []

workflow PHYLOPHOENIX {

        main:
        ch_versions = Channel.empty()
        // Allow outdir to be relative
        //outdir_path = Channel.fromPath(params.outdir, relative: true, type: 'dir')

        // Create report
        GRIPHIN_WORKFLOW (
            input_samplesheet_path, params.input_dir
        )
        ch_versions = ch_versions.mix(GRIPHIN_WORKFLOW.out.versions)

        /*/ Creates samplesheets with sampleid,seq_type,path_to_assembly
        GET_COMPARISONS (
            GRIPHIN_WORKFLOW.out.griphin_samplesheet
        )
        ch_versions = ch_versions.mix(GET_COMPARISONS.out.versions)

        // creates channel: [ val(meta.id, meta.st), [ scaffolds_1, scaffolds_2 ] ]
        CREATE_META (
            GET_COMPARISONS.out.samplesheet, GET_COMPARISONS.out.snv_samplesheet
        )

        // Run Mash on groups of samples by seq type
        MASH_DIST (
            CREATE_META.out.st_scaffolds
        )
        ch_versions = ch_versions.mix(MASH_DIST.out.versions)

        // Creating channel [ ST, [distance_1, distance_2] ]
        dist_ch = MASH_DIST.out.dist.map{ meta, dist -> [ dist ] }.collect() // drop meta and collect all distance files
        centroid_ch = dist_ch.map{ dist -> [ "All_STs", dist ]} // add back "All_STs" as the meta value
        .combine(GRIPHIN_WORKFLOW.out.griphin_samplesheet) // Add samplesheet to all mash distance channels

        // Take in all mash distance files  then use the samplesheet to return the centroid assembly
        // Get centroid, by calculating the average mash distance
        GET_CENTROID (
            centroid_ch
        )
        ch_versions = ch_versions.mix(GET_CENTROID.out.versions)

        // Unzip centroid assembly: SNVPhyl requires it unzipped
        ASSET_CHECK(
            GET_CENTROID.out.centroid_path.splitCsv( header:false, sep:',' ) // Bring in centroid into channel
        )

        // Make SNVPHYL channel by joining by seq type
        all_ch = CREATE_META.out.st_snv_samplesheets.join(ASSET_CHECK.out.unzipped_fasta, by: [0])

        // Run snvphyl on each st type on its own input
        SNVPHYL (
            all_ch.map{ seq_type, samplesheet, unzipped_fasta -> [seq_type, samplesheet]}, all_ch.map{ seq_type, samplesheet, unzipped_fasta -> [seq_type, unzipped_fasta] }
        )
        ch_versions = ch_versions.mix(SNVPHYL.out.versions)*/

        // If you pass --both then samples will be broken up by st type and SNVPhyl run on each st on its own
        if (params.both==true) {

            // Creates samplesheets with sampleid,seq_type,path_to_assembly
            GET_SEQUENCE_TYPES (
                GRIPHIN_WORKFLOW.out.griphin_samplesheet, GRIPHIN_WORKFLOW.out.griphin_report
            )
            ch_versions = ch_versions.mix(GET_SEQUENCE_TYPES.out.versions)

            // creates channel: [ val(meta.id, meta.st), [ scaffolds_1, scaffolds_2 ] ]
            CREATE_META_BY_ST (
                GET_SEQUENCE_TYPES.out.st_samplesheets, GET_SEQUENCE_TYPES.out.st_snv_samplesheets
            )

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
            centroid_st_ch =  st_mash_dists.combine(GRIPHIN_WORKFLOW.out.griphin_samplesheet)

            // Take in all mash distance files for a seq type and then use the samplesheet to return the centroid assembly
            // Get centroid for ST groups, by calculating the average mash distance
            GET_CENTROID_BY_ST (
                centroid_st_ch
            )
            ch_versions = ch_versions.mix(GET_CENTROID_BY_ST.out.versions)

            // Unzip centroid assembly: SNVPhyl requires it unzipped
            ASSET_CHECK_BY_ST(
                GET_CENTROID_BY_ST.out.centroid_path.splitCsv( header:false, sep:',' ) // Bring in centroid into channel
            )

            // Make SNVPHYL channel by joining by seq type
            st_ch = CREATE_META_BY_ST.out.st_snv_samplesheets.join(ASSET_CHECK_BY_ST.out.unzipped_fasta, by: [0,0])

            // Run snvphyl on each st type on its own input
            SNVPHYL_BY_ST (
                st_ch.map{ seq_type, samplesheet, unzipped_fasta -> [seq_type, samplesheet]}, st_ch.map{ seq_type, samplesheet, unzipped_fasta -> [seq_type, unzipped_fasta] }
            )
            ch_versions = ch_versions.mix(SNVPHYL_BY_ST.out.versions)
        }

        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

        /*/
        // MODULE: MultiQC
        //
        workflow_summary    = WorkflowGriphin.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)
        ch_workflow_summary.view()

        methods_description    = WorkflowGriphin.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
        ch_methods_description = Channel.value(methods_description)
        ch_methods_description.view()*/

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
