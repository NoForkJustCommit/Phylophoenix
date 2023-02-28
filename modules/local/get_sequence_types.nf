process GET_SEQUENCE_TYPES {
    label 'process_low'
    //container 'quay.io/jvhagey/phoenix:base_v1.1.0' install openpyxl

    input:
    path(griphin_samplesheet) // -s
    path(griphin_report) // -g

    output:
    path("ST*_samplesheet.csv"),         emit: st_samplesheets     // headers: id,seq_type,assembly_1,assembly_2
    path("SNVPhyl_ST*_samplesheet.csv"), emit: st_snv_samplesheets // headers: id,directory
    path("versions.yml"),                emit: versions

    script: // This script is bundled with the pipeline, in dhqp/griphin/bin/
    """
    get_st_types_new.py -g ${griphin_report} -s ${griphin_samplesheet}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}