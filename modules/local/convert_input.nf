/* in the case of a given directory we can get reads from that easy. If a samplesheet is given then the directory is there. */
process CONVERT_INPUT {
    tag "${meta.seq_type}"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

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
    """
    convert_samplesheet.py ${samplesheet} -t ${meta.seq_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}