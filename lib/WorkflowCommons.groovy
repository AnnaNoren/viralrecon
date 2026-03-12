//
// This file holds several functions common to the multiple workflows in the nf-core/viralrecon pipeline
//
import nextflow.Nextflow

class WorkflowCommons {

    //
    // Function to get column entries from a file
    //
    public static ArrayList getColFromFile(input_file, col=0, uniqify=false, sep='\t') {
        def vals = []
        input_file.eachLine { line ->
            def val = line.split(sep)[col]
            if (uniqify) {
                if (!vals.contains(val)) {
                    vals << val
                }
            } else {
                vals << val
            }
        }
        return vals
    }

    //
    // Function that returns the number of lines in a file
    //
    public static Integer getNumLinesInFile(input_file) {
        def num_lines = 0
        input_file.eachLine { line ->
            num_lines ++
        }
        return num_lines
    }

    //
    // Function to generate an error if contigs in BED file do not match those in reference genome
    //
    public static void checkContigsInBED(fai_contigs, bed_contigs, log) {
        def intersect = bed_contigs.intersect(fai_contigs)
        if (intersect.size() != bed_contigs.size()) {
            def diff = bed_contigs.minus(intersect).sort()
            Nextflow.error("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Contigs in primer BED file do not match those in the reference genome:\n\n" +
                "  ${diff.join('\n  ')}\n\n" +
                "  Please check:\n" +
                "    - Primer BED file supplied with --primer_bed\n" +
                "    - Genome FASTA file supplied with --fasta\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        }
    }

    //
    // Function to read in all fields into a Groovy Map from Nextclade CSV output file
    //
    // See: https://stackoverflow.com/a/67766919
    public static Map getNextcladeFieldMapFromCsv(nextclade_report) {
        def headers   = []
        def field_map = [:]
        nextclade_report.readLines().eachWithIndex { row, row_index ->
            def vals = row.split(';')
            if (row_index == 0) {
                headers = vals
            } else {
                def cells = headers.eachWithIndex { header, header_index ->
                    def val = (header_index <= vals.size()-1) ? vals[header_index] : ''
                    field_map[header] = val ?: 'NA'
                }
            }
        }
        return field_map
    }

    //
    // Function to get number of variants reported in BCFTools stats file
    //
    public static Integer getNumVariantsFromBCFToolsStats(bcftools_stats) {
        def num_vars = 0
        bcftools_stats.eachLine { line ->
            def matcher = line =~ /SN\s*0\s*number\sof\srecords:\s*([\d]+)/
            if (matcher) num_vars = matcher[0][1].toInteger()
        }
        return num_vars
    }
}
