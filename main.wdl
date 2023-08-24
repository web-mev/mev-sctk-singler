workflow SctkSingler {

    # The input matrix
    File exprs_file

    # 'main' or 'fine', which sets how fine-grained
    # the cell type calls are
    String annotation_level

    # which cell marker reference to use
    String reference

    # which identifier is used in the exp matrix? ensg or symbol
    String identifier_choice

    call runSingler {
        input:
            exprs_file = exprs_file,
            annotation_level = annotation_level,
            reference = reference,
            identifier_choice = identifier_choice
    }

    output {
        File cell_assignments = runSingler.cell_assignments
    }

}

task runSingler {
        
    File exprs_file

    # 'main' or 'fine', which sets how fine-grained
    # the cell type calls are
    String annotation_level

    # which cell marker reference to use
    String reference

    # which identifier is used in the exp matrix? ensg or symbol
    String identifier_choice

    String output_prefix = "singler"

    Int disk_size = 20

    command <<<
        Rscript /opt/software/singler.R \
            -f ${exprs_file} \
            -o ${output_prefix} \
            --level ${annotation_level} \
            --reference "${reference}" \
            --featureType ${identifier_choice}
    >>>

    output {
        File cell_assignments = "${output_prefix}.cell_types.${reference}.${annotation_level}.tsv"
    }

    runtime {
        docker: "ghcr.io/web-mev/mev-sctk-singler"
        cpu: 4
        memory: "10 G"
        disks: "local-disk " + disk_size + " HDD"
    }
}



