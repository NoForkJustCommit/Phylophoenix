process GET_COMPARISONS {
    label 'process_low'
    //container 'quay.io/jvhagey/phoenix:base_v1.1.0' install openpyxl

    input:
    path(griphin_samplesheet) // -s

    output:
    path("All_STs_samplesheet.csv"),         emit: samplesheet     // headers: id,seq_type,assembly_1,assembly_2
    path("SNVPhyl_All_STs_samplesheet.csv"), emit: snv_samplesheet // headers: id,directory
    path("versions.yml"),                    emit: versions

    script: // This script is bundled with the pipeline, in dhqp/griphin/bin/
    """
    create_comparisions.py -s ${griphin_samplesheet}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}