/* Find repeats and create invaild_positions.bed */

def VERSION = '1.8.2' // Version information not provided by in container

process FIND_REPEATS {
    tag "${meta.seq_type}"
    label 'process_low'
    container "staphb/snvphyl-tools:1.8.2"

    input:
    tuple val(meta), path(refgenome)

    output:
    path("invalid_positions.bed"),    emit: repeats_bed_file
    path("versions.yml"),             emit: versions

    script:
    """
    find-repeats.pl ${refgenome} --min-length 150 --min-pid 90 > invalid_positions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snvphyl-tools: $VERSION
        perl: \$(perl --version | grep "This is perl" | sed 's/.*(\\(.*\\))/\\1/' | cut -d " " -f1)
    """
}