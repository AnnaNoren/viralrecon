//
// Consensus calling QC
//

include { QUAST                  } from '../../../modules/nf-core/quast/main'
include { PANGOLIN_UPDATEDATA    } from '../../../modules/nf-core/pangolin/updatedata/main'
include { PANGOLIN_RUN           } from '../../../modules/nf-core/pangolin/run/main'
include { NEXTCLADE_RUN          } from '../../../modules/nf-core/nextclade/run/main'
include { PLOT_BASE_DENSITY      } from '../../../modules/local/plot_base_density'
include { UNTAR as UNTAR_PANGODB } from '../../../modules/nf-core/untar/main'

workflow CONSENSUS_QC {
    take:
    consensus    // channel: [ val(meta), [ consensus ] ]
    fasta        // channel: /path/to/genome.fasta
    gff          // channel: /path/to/genome.gff
    nextclade_db // channel: /path/to/nextclade_db/

    main:

    ch_versions = channel.empty()

    //
    // Consensus QC report across samples with QUAST
    //
    ch_quast_results = channel.empty()
    ch_quast_tsv     = channel.empty()
    if (!params.skip_variants_quast) {
    consensus
        .collect{ it[1] }
        .map { consensus_collect -> tuple([id: "quast"], consensus_collect) }
        .set { ch_to_quast }

        QUAST (
            ch_to_quast,
            fasta.map { [ [:], it ] },
            gff
        )
        ch_quast_results = QUAST.out.results
        ch_quast_tsv     = QUAST.out.tsv
        ch_versions      = ch_versions.mix(QUAST.out.versions)
    }

    //
    // Lineage analysis with Pangolin
    //
    ch_pangolin_report = channel.empty()
    ch_pango_database = channel.empty()

    if (!params.skip_pangolin) {
        if (!params.pango_database) {
            PANGOLIN_UPDATEDATA('pangolin_db')
            ch_pango_database = PANGOLIN_UPDATEDATA.out.db
            ch_versions       = ch_versions.mix(PANGOLIN_UPDATEDATA.out.versions.first())
        } else {
            if (params.pango_database.endsWith('.tar.gz')) {
                UNTAR_PANGODB (
                    [ [:], params.pango_database ]
                )
                ch_pango_database = UNTAR_PANGODB.out.untar.map { it[1] }
                ch_versions       = ch_versions.mix(UNTAR_PANGODB.out.versions)
            } else {
                ch_pango_database = channel.value(file(params.pango_database, type: 'dir'))
            }
        }
        PANGOLIN_RUN (
            consensus,
            ch_pango_database
        )
        ch_pangolin_report = PANGOLIN_RUN.out.report
        ch_versions        = ch_versions.mix(PANGOLIN_RUN.out.versions.first())
    }

    //
    // Lineage analysis with Nextclade
    //
    ch_nextclade_report = channel.empty()
    if (!params.skip_nextclade) {
        NEXTCLADE_RUN (
            consensus,
            nextclade_db
        )
        ch_nextclade_report = NEXTCLADE_RUN.out.csv
        ch_versions         = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())
    }

    //
    // Plot consensus base density
    //
    ch_bases_tsv = channel.empty()
    ch_bases_pdf = channel.empty()
    if (!params.skip_consensus_plots) {
        PLOT_BASE_DENSITY (
            consensus
        )
        ch_bases_tsv = PLOT_BASE_DENSITY.out.tsv
        ch_bases_pdf = PLOT_BASE_DENSITY.out.pdf
    }

    emit:
    quast_results    = ch_quast_results    // channel: [ val(meta), [ results ] ]
    quast_tsv        = ch_quast_tsv        // channel: [ val(meta), [ tsv ] ]

    pangolin_report  = ch_pangolin_report  // channel: [ val(meta), [ csv ] ]

    nextclade_report = ch_nextclade_report // channel: [ val(meta), [ csv ] ]

    bases_tsv        = ch_bases_tsv        // channel: [ val(meta), [ tsv ] ]
    bases_pdf        = ch_bases_pdf        // channel: [ val(meta), [ pdf ] ]

    versions         = ch_versions         // channel: [ versions.yml ]
}
