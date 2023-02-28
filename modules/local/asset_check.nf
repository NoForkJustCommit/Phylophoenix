process ASSET_CHECK {
    tag "${meta.seq_type}"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(zipped_fasta)

    output:
    tuple val(meta), path("$gunzip"), emit: unzipped_fasta

    script:
    gunzip = zipped_fasta.toString() - '.gz'
    """
    if [[ ${zipped_fasta} == *.gz ]]
    then
        gunzip --force ${zipped_fasta}
    else
        :
    fi
    """
}
