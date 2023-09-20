/* Mpileup */
process MPILEUP {
    tag "${meta.id}_${meta.seq_type}"
    label 'process_low'
    container "staphb/bcftools:1.14"

    input:
    tuple val(meta), path(sorted_bams)
    tuple val(meta_2), path(ref_fai), path(ref_sma), path(ref_smi)
    tuple val(meta_3), path(ref_genome)

    output:
    tuple val(meta), path( "${meta.id}_mpileup.vcf" ), emit: mpileup
    path("versions.yml"),                              emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #confirm that ST types are the same
    if [[ "${meta.seq_type}" != "${meta_2.seq_type}" ]] && [[ "${meta.seq_type}" != "${meta_3.seq_type}" ]]
    then
        echo "Yikes, the ST type of your reference and query do not match. This shouldn't happen, please report the bug by opening a github issue."
        exit 1
    fi

    bcftools mpileup --threads 4 --fasta-ref ${ref_genome} -A -B -C 0 -d 1024 -q 0 -Q 0 --output-type v -I --output ${prefix}_mpileup.vcf ${sorted_bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}