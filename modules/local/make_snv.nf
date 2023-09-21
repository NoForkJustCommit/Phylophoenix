/* Making SNV matrix */
process MAKE_SNV {
    tag "${meta.seq_type}"
    label 'process_low'
    container "staphb/snvphyl-tools:1.8.2"

    input:
    tuple val(meta), path(snvAlignment_phy)

    output:
    tuple val(meta), path('snvMatrix.tsv'), emit: snvMatrix
    path("versions.yml"),                   emit: versions

    script:
    def container = task.container.toString() - "staphb/snvphyl-tools:"
    """
    snv_matrix.pl ${snvAlignment_phy} -o snvMatrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snvphyl-tools: ${container}
        perl: \$(perl --version | grep "This is perl" | sed 's/.*(\\(.*\\))/\\1/' | cut -d " " -f1)
    END_VERSIONS
    """
}