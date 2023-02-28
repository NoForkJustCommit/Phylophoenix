//
// Getting ST types in order
//

workflow CREATE_META {
    take:
        samplesheets     // headers: id,seq_type,assembly_1,assembly_2
        snv_samplesheets // headers: id,directory

    main:
        // flatten so some samplesheet is  to create [ meta.id, meta.st] then group by st type collect them all and pass again to create list for each ST type
        st_scaffolds = samplesheets.flatten() //flatten so one samplesheet goes through at a time
            .splitCsv( header:true, sep:',' ) // split incoming samplesheet one row at a time
            .map{ create_assembly_channel(it) } // creates create [ [ meta.id, meta.st], assembly_1, assembly_2 ] for reach sample

        st_samplesheets = samplesheets.flatten()
            .map{ add_st_to_samplesheet(it, "st") }

        st_snv_samplesheets = snv_samplesheets.flatten()
            .map{ add_st_to_samplesheet(it, "st_snv") }

    emit:
        st_scaffolds        // channel: [ val(meta), [ scaffolds_1, scaffolds_2 ] ]
        st_samplesheets     // channel: [ seq_type, samplesheet ]
        st_snv_samplesheets // channel: [ seq_type, samplesheet ]
}

def add_st_to_samplesheet(samplesheet, samplesheet_type) {

    if (samplesheet_type == "st") {
        seq_type = samplesheet.toString().replaceAll("_samplesheet.csv", "") // remove _samplesheet.csv from path 
        seq_type = seq_type.split('/')[-1] // get the last string after the last backslash
    }
    if (samplesheet_type == "st_snv") {
        seq_type = samplesheet.toString().replaceAll("SNVPhyl_", "") // remove _samplesheet.csv from path
        seq_type = seq_type.replaceAll("_samplesheet.csv", "") // remove _samplesheet.csv from path 
        seq_type = seq_type.split('/')[-1] // get the last string after the last backslash
    }
    def meta = [:]
    meta.seq_type = seq_type
    new_samplesheet = [ meta, file(samplesheet) ]
    return new_samplesheet
}

// Function to get list of [ meta, [ scaffolds_1, scaffolds_2 ] ]
def create_assembly_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id = row.id
    meta.seq_type = row.seq_type

    // add path(s) of the assembly file(s) to the meta map
    def assembly_meta = []
    if (!file(row.assembly_1).exists()) {
        println(row.assembly_1)
        exit 1, "ERROR: Please check st samplesheet -> Assembly scaffolds file does not exist!\n${row.assembly}"
    }
    if (!file(row.assembly_2).exists()) {
        println(row.assembly_2)
        exit 1, "ERROR: Please check st samplesheet -> Assembly scaffolds file does not exist!\n${row.assembly}"
    }
    assembly_meta = [ meta, file(row.assembly_1), file(row.assembly_2) ]
    return assembly_meta
}


/*def create_st_tuples(input) {
    // the output of this function is to get a tuple like this:
    // [[[id:sample1, st:ST1], $PATH/sample1.filtered.scaffolds.fa.gz], [id:sample2, st:ST1], $PATH/sample2.filtered.scaffolds.fa.gz], [[id:sample3, st:ST2], $PATH/sample3.filtered.scaffolds.fa.gz], [id:sample4, st:ST2], $PATH/sample4.filtered.scaffolds.fa.gz]]
    //println(input)
    def starting_st = "ST" // initializing st
    count = 0
    sample_num = 1
    for (item in input) { // loop through each st type
        if ( count == 0 || count % 2 == 0) { // even number or zero it will be meta information, odd number will be the sample info 
            // get the current st
            new_st = item.st
            if (starting_st == new_st) { // if the new st is the same as the last st
                // add it to the tuple that was created before
                new_st_tuple.add(item)
                new_st_tuple.add(input.get(sample_num)[0])
            } else {
                // if the incoming st is a new one
                if (count == 0 ) { // check if this is the first time through and create an empty list if so
                    final_tuple = []
                } else { // if this is not the first time through add the completed tuple to the final list
                    // add the last tuple the final tuple then
                    final_tuple.add(new_st_tuple)
                }
                // set new st type
                starting_st = new_st
                // create a new tuble
                new_st_tuple = [new_st]
                // add current sample to this table
                new_st_tuple.add(item)
                new_st_tuple.add(input.get(sample_num)[0])
            }
            sample_num = sample_num + 2
            count = count + 1
        } else { // if count is an even number then only just add to variables
            count = count + 1
        }
    }
    return final_tuple
}*/