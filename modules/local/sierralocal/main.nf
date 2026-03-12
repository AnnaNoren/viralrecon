process SIERRALOCAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sierra-local:0.4.3--py310hdfd78af_0' :
        'biocontainers/sierra-local:0.4.3--py310hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path hivdb_xml
    path apobec_drm
    path apobec_csv
    path unusual_csv
    path sdrms_csv
    path mutation_csv
    output:
    tuple val(meta), path("*_resistance.json") , emit: json
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix              = task.ext.prefix     ?: "${meta.id}"
    def args                = task.ext.args       ?: ''
    def hivdb_xml           = params.hivdb_xml    ? "-xml  ${params.hivdb_xml}"            : ''
    def apobec_drm          = params.apobec_drm   ? "-json ${params.apobec_drm}"           : ''
    def apobec_csv          = params.apobec_csv   ? "-apobec_csv ${params.apobec_csv}"     : ''
    def unusual_csv         = params.unusual_csv  ? "-unusual_csv ${params.unusual_csv}"   : ''
    def sdrms_csv           = params.sdrms_csv    ? "-sdrms_csv ${params.sdrms_csv}"       : ''
    def mutation_csv        = params.mutation_csv ? "-mutation_csv ${params.mutation_csv}" : ''
    def sierralocal_version = "0.4.3"

    """
    sierralocal \\
        $args \\
        $hivdb_xml \\
        $apobec_drm \\
        $apobec_csv \\
        $unusual_csv \\
        $sdrms_csv \\
        $mutation_csv \\
        -o ${prefix}_resistance.json \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sierra-local: \$(echo "${sierralocal_version}")
    END_VERSIONS
    """
}
