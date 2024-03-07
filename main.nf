process run_singler {

    tag "Run SCTK SingleR cell ID calling"
    publishDir "${params.output_dir}/SctkSingler.cell_assignments", mode:"copy", pattern:"${output_prefix}*"
    container "ghcr.io/web-mev/mev-sctk-singler"
    cpus 4
    memory '10 GB'

    input:
        path exprs_file

    output:
        path "${output_prefix}*"

    script:
        output_prefix = "singler"
        """
        Rscript /opt/software/singler.R \
            -f ${exprs_file} \
            -o ${output_prefix} \
            --level ${params.annotation_level} \
            --reference "${params.reference}" \
            --featureType ${params.identifier_choice}
        """
}

workflow {
    run_singler(params.exprs_file)
}