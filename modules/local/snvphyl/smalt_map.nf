/* Map reads to reference genome & create BAM file */
process SMALT_MAP {
    tag "${meta.id}_${meta.seq_type}"
    label 'process_medium'
    container "staphb/smalt:0.7.6"

    input:
    tuple val(meta), file(reads)
    tuple val(meta_2), path(ref_fai), path(ref_sma), path(ref_smi)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bams
    path("versions.yml"),                    emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def ref_name = ref_fai.toString() - '.fai'
    """
    #confirm that ST type for the query and reference match
    if [[ "${meta.seq_type}" != "${meta_2.seq_type}" ]]
    then
        echo "Yikes, the ST type of your reference and query do not match. This shouldn't happen, please report the bug by opening a github issue."
        exit 1
    fi

    smalt map -f bam -n 4 -l pe -i 1000 -j 20 -r 1 -y 0.5 -o ${prefix}.bam ${ref_name} ${reads[0]} ${reads[1]}
    smalt version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smalt: \$( smalt version | grep "Version:" | sed 's/Version://' )
    END_VERSIONS
    """
}