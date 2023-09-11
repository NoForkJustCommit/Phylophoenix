/* Generate a line for the next process */
process GENERATE_LINE_2 {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    path(consolidated_bcf)

    output:
    path("consolidation_line.txt"), emit: consolidation_line
    path("versions.yml"),           emit: versions

    script:
    // get container info
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    for f in *_consolidated.bcf
    do
    fname=\$(basename \$f _consolidated.bcf)
    echo "--consolidate_vcf \$fname=\$f " | tr -d "\\n" >> consolidation_line.txt
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}