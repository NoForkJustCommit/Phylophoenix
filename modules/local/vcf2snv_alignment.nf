/* VCF2SNV_ALIGNMENT Call variants */
process VCF2SNV_ALIGNMENT {
    tag "${meta.seq_type}"
    label 'process_medium'
    container "staphb/snvphyl-tools:1.8.2"

    input:
    tuple val(meta), val(consolidate_bcfs)
    tuple val(meta), path(bcf)
    tuple val(meta), path(new_invalid_positions)
    tuple val(meta), path(refgenome)
    tuple val(meta), path(consolidated_bcf_index)

    output:
    tuple val(meta), path("${meta.seq_type}_snvAlignment.phy"), emit: snvAlignment
    tuple val(meta), path("${meta.seq_type}_vcf2core.tsv"),     emit: vcf2core
    tuple val(meta), path("${meta.seq_type}_snvTable.tsv"),     emit: snvTable
    path("versions.yml"),                                       emit: versions

    script:
    def container = task.container.toString() - "staphb/snvphyl-tools:"
    """
    vcf2snv_alignment.pl --reference reference --invalid-pos ${new_invalid_positions} --format fasta --format phylip --numcpus 4 --output-base snvalign --fasta ${refgenome} ${consolidate_bcfs} 
    mv snvalign-positions.tsv ${meta.seq_type}_snvTable.tsv
    mv snvalign-stats.csv ${meta.seq_type}_vcf2core.tsv
    if [[ -f snvalign.phy ]]; then
        mv snvalign.phy ${meta.seq_type}_snvAlignment.phy
        sed -i "s/'//" ${meta.seq_type}_snvAlignment.phy
        sed -i "s/'//" ${meta.seq_type}_snvAlignment.phy
    else
        touch ${meta.seq_type}_snvAlignment.phy
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snvphyl-tools: ${container}
        perl: \$(perl --version | grep "This is perl" | sed 's/.*(\\(.*\\))/\\1/' | cut -d " " -f1)
    END_VERSIONS
    """
}