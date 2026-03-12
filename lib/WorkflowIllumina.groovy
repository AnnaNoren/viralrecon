//
// This file holds several functions specific to the workflow/illumina.nf in the nf-core/viralrecon pipeline
//
import nextflow.Nextflow
import groovy.json.JsonSlurper

class WorkflowIllumina {

    //
    // Print warning if genome fasta has more than one sequence
    //
    public static void isMultiFasta(fasta_file, log) {
        def count = 0
        def line  = null
        fasta_file.withReader { reader ->
            while (line = reader.readLine()) {
                if (line.contains('>')) {
                    count++
                    if (count > 1) {
                        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                            "  This pipeline does not officially support multi-fasta genome files!\n\n" +
                            "  The parameters and processes are tailored for viral genome analysis.\n" +
                            "  Please amend the '--fasta' parameter.\n" +
                            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                        break
                    }
                }
            }
        }
    }

    //
    // Check if the primer BED file supplied to the pipeline is from the SWIFT/SNAP protocol
    //
    public static void checkIfSwiftProtocol(primer_bed_file, name_prefix, log) {
        def count = 0
        def line  = null
        primer_bed_file.withReader { reader ->
            while (line = reader.readLine()) {
                def name = line.split('\t')[3]
                if (name.contains(name_prefix)) {
                    count++
                    if (count > 1) {
                        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                            "  Found '${name_prefix}' in the name field of the primer BED file!\n" +
                            "  This suggests that you have used the SWIFT/SNAP protocol to prep your samples.\n" +
                            "  If so, please set '--ivar_trim_offset 5' as suggested in the issue below:\n" +
                            "  https://github.com/nf-core/viralrecon/issues/170\n" +
                            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                        break
                    }
                }
            }
        }
    }

    //
    // Function that parses fastp json output file to get total number of reads after trimming
    //
    public static Integer getFastpReadsAfterFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['after_filtering']['total_reads'].toInteger()
    }

    //
    // Function that parses fastp json output file to get total number of reads before trimming
    //
    public static Integer getFastpReadsBeforeFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['before_filtering']['total_reads'].toInteger()
    }
}
