/* Freebayes variant calling */
process FREEBAYES {
    tag "${meta.id}_${meta.seq_type}"
    label 'process_low'
    container "staphb/freebayes:1.3.6"

    input:
    tuple val(meta), path(sorted_bams)
    tuple val(meta_2), path(ref_fai), path(ref_sma), path(ref_smi)
    tuple val(meta_3), path(ref_genome)

    output:
    tuple val(meta), path( "${meta.id}_freebayes.vcf" ), emit: vcf_files
    path("versions.yml"),                                emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #confirm that ST type for the query and reference match
    if [[ "${meta.seq_type}" != "${meta_2.seq_type}" ]] && [[ "${meta.seq_type}" != "${meta_3.seq_type}" ]]
    then
        echo "Yikes, the ST type of your reference and query do not match. This shouldn't happen, please report the bug by opening a github issue."
        exit 1
    fi

    freebayes --bam ${sorted_bams} --ploidy 1 --fasta-reference ${ref_genome} --vcf ${prefix}_freebayes.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}