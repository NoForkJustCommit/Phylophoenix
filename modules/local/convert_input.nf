/* in the case of a given directory we can get reads from that easy. If a samplesheet is given then the directory is there. */
process CONVERT_INPUT {
    tag "${meta.seq_type}"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(sample_sheet)
    //path(directory)

    output:
    path("updated_samplesheet.csv"), emit: updated_samplesheet
    path("versions.yml"),            emit: versions

    script:
    def samplesheet = sample_sheet ? "--samplesheet ${sample_sheet}" : ""
    //def input_dir   = directory ? "--input_dir ${directory}" : ""
    //def report_prefix = prefix ? "--output ${prefix}" : ""
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) {
        ica = ""
    } else if (params.ica==true) {
        ica = "pearl ${workflow.launchDir}/bin/"
    } else {
        error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods."
    }
    // get container info
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    ${ica}convert_samplesheet.py ${samplesheet} -t ${meta.seq_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}