process GRIPHIN {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(sample_sheet) // -s
    path(db) // -a
    val(prefix) // -o
    path(control_list) // -c
    val(complete_path)

    output:
    path("*.xlsx"),                  emit: griphin_report
    path("GRiPHin_samplesheet.csv"), emit: griphin_samplesheet
    path("versions.yml"),            emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/griphin/bin/
    //def samplesheet = sample_sheet ? "--samplesheet ${sample_sheet}" : ""
    def controls    = control_list ? "--control_list ${control_list}" : ""
    //def input_dir   = input ? "--directory ${input}" : ""
    def report_prefix     = prefix ? "--output ${prefix}" : ""
    """
    #create a samplesheet to be passed to GRiPHin.py, this helps to make sure the paths are full. 
    if [ ! -f GRiPHin_samplesheet.csv ]
    then
        while IFS="" read -r line;
        do
            sample_name=\$(echo \$line | cut -d ',' -f 1)
            echo
            if [[ "\$sample_name" == "sample" ]]; then
                echo "sample,directory" > GRiPHin_samplesheet.csv
            else
                #get the full path for the samples rather than the working directory
                full_path=\$(readlink -f ${complete_path}/Phoenix_Output_Report.tsv)
                full_dir=\$(echo \$full_path | sed 's/\\/Phoenix_Output_Report.tsv//')
                echo \$sample_name,\$full_dir/\$sample_name >> GRiPHin_samplesheet.csv
            fi
        done < ${sample_sheet}
    fi

    GRiPHin.py --samplesheet GRiPHin_samplesheet.csv --ar_db ${db} ${controls} ${report_prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}