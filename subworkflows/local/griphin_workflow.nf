//
// Subworkflow: one of the workflows to generate a report
//

include { GRIPHIN            } from '../../modules/local/griphin'
include { CREATE_SAMPLESHEET } from '../../modules/local/create_samplesheet'


workflow GRIPHIN_WORKFLOW {
    take:
        input_samplesheet_path // channel: tuple val(meta), path('*.json'): FASTP_TRIMD.out.json --> PHOENIX_EXQC.out.paired_trmd_json
        input_dir

    main:
        ch_versions     = Channel.empty() // Used to collect the software versions

        // Check if input samplesheet or input directory has passed, then check if a control list was passed
        if (params.input != null) {
            if (params.control_list != null){
                // Allow control list to be relative
                control_path = Channel.fromPath(params.control_list, relative: true)
                // Create report
                GRIPHIN (
                    input_samplesheet_path, params.ardb, params.prefix, control_path, []
                )
                ch_versions = ch_versions.mix(GRIPHIN.out.versions)
            } else {
                // Create report
                GRIPHIN (
                    input_samplesheet_path, params.ardb, params.prefix, [], []
                )
                ch_versions = ch_versions.mix(GRIPHIN.out.versions)
            }
        } else {
            // allow input directory to be relative
            inputdir_path = Channel.fromPath(params.input_dir, relative: true, type: 'dir') // this same path is needed to make the samplesheet

            // Create samplesheet
            CREATE_SAMPLESHEET (
                inputdir_path
            )
            ch_versions = ch_versions.mix(CREATE_SAMPLESHEET.out.versions)

            if (params.control_list != null){
                // Allow control list to be relative
                control_path = Channel.fromPath(params.control_list, relative: true)
                // Create report
                GRIPHIN (
                    CREATE_SAMPLESHEET.out.samplesheet, params.ardb, params.prefix, control_path, inputdir_path
                )
                ch_versions = ch_versions.mix(GRIPHIN.out.versions)
            } else {
                // Create report
                GRIPHIN (
                    CREATE_SAMPLESHEET.out.samplesheet, params.ardb, params.prefix, [], inputdir_path
                )
                ch_versions = ch_versions.mix(GRIPHIN.out.versions)
            }
        }

    emit:
        griphin_report      = GRIPHIN.out.griphin_report
        griphin_samplesheet = GRIPHIN.out.griphin_samplesheet
        versions            = ch_versions // channel: [ versions.yml ]
}