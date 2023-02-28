process GET_CENTROID {
    tag "${meta.seq_type}"
    label 'process_low'
    //container 'quay.io/jvhagey/phoenix:base_v1.1.0' install openpyxl

    input:
    tuple val(meta), path(mash_distance), \
    path(griphin_samplesheet) // -s

    output:
    tuple val(meta), path("path_to_${meta.seq_type}_centroid.csv"), emit: centroid_path
    path("versions.yml"),                                           emit: versions

    script: // This script is bundled with the pipeline, in dhqp/griphin/bin/
    """
    # combine all lists
    cat *.txt > ${meta.seq_type}_dists.tsv

    get_centroid.py -i ${meta.seq_type}_dists.tsv -s ${griphin_samplesheet} -t ${meta.seq_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}