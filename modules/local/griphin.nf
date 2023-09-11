process GRIPHIN {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    path(sample_sheet) // -s . This is an empty list when no samplesheet is passed
    path(db) // -a
    val(prefix) // -o
    path(control_list) // -c
    val(coverage)
    val(entry)

    output:
    path("*_Summary.xlsx"),            emit: griphin_report
    path("*_Summary.tsv"),             emit: griphin_tsv_report
    path("Directory_samplesheet.csv"), optional: true, emit: griphin_samplesheet
    path("versions.yml"),              emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/griphin/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) {
        ica = ""
    } else if (params.ica==true) {
        ica = "python ${workflow.launchDir}/bin/"
    } else {
        error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods."
    }
    def samplesheet   = sample_sheet ? "--samplesheet ${sample_sheet}" : ""
    def controls      = control_list ? "--control_list ${control_list}" : ""
    def report_prefix = prefix ? "--output ${prefix}" : ""
    def phoenix       = entry ? "" : "--phoenix" // tells griphin in the run was a CDC (i.e. -entry CDC_PHOENIX) one or standard
    // get container info
    def container     = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """

    ${ica}GRiPHin.py --ar_db ${db} --coverage ${coverage} ${samplesheet} ${controls} ${phoenix} ${report_prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}