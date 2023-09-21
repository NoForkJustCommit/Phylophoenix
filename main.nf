#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cdcgov/phylophoenix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/cdcgov/phylophoenix
    Slack  : https://staph-b-dev.slack.com/channels/phoenix-h-dev
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)


//if you use --no_all phylophoenix assumes you want to do it by st and will set to true
if (params.no_all==true) {
    params.by_st = true
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PHYLOPHOENIX } from './workflows/phylophoenix'

//
// WORKFLOW: Run main nf-core/phylophoenix analysis pipeline
//
workflow PHYLOPHOENIX_WF {
    // Check input path parameters to see if they exist
    if (params.input != null ) {  // if a samplesheet is passed
        // allow input to be relative, turn into string and strip off the everything after the last backslash to have remainder of as the full path to the samplesheet. 
        //input_samplesheet_path = Channel.fromPath(params.input, relative: true).map{ [it.toString().replaceAll(/([^\/]+$)/, "").replaceAll(/\/$/, "") ] }
        input_samplesheet_path = Channel.fromPath(params.input, relative: true)
        if (params.input_dir != null ) { //if samplesheet is passed and an input directory exit
            exit 1, 'You need EITHER an input samplesheet or a directory! Just pick one.' 
        } else { // if only samplesheet is passed check to make sure input is an actual file
            input_dir = []
            def checkPathParamList = [ params.input ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
        }
    } else { // if no samplesheet is passed
        if (params.input_dir != null ) { // if no samplesheet is passed, but an input directory is given
            input_samplesheet_path = []
            def checkPathParamList = [ params.input_dir ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            input_dir = params.input_dir
        } else { // if no samplesheet is passed and no input directory is given
            exit 1, 'You need EITHER an input samplesheet or a directory!' 
        }
    }

    main:
        PHYLOPHOENIX ( input_samplesheet_path, input_dir )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    PHYLOPHOENIX_WF ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
