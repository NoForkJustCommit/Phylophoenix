/* Generate a line for the next process */
process GENERATE_LINE_1 {
    label 'process_low'
    stageInMode 'copy'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    path(sorted_bams)

    output:
    path("bam_line.txt"), emit: bam_lines_file
    path("versions.yml"), emit: versions

    script:
    // get container info
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    count=0
    for f in *_sorted.bam; do
        ((count++))
        if [[ count == 1 ]]; then 
            echo "--bam bam\$count=./\$f " | tr -d "\\n" > bam_line.txt
        else
            echo "--bam bam\$count=./\$f " | tr -d "\\n" >> bam_line.txt
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}