nextflow.enable.dsl=2

workflow {
    emit:
        filtered

    main:
        Channel.fromPath(params.x).set{blast_out}

        filtered = filter_blast(blast_out, params.min_blast_aln_len, params.min_blast_aln_frac)
}

process filter_blast_process {
    tag "${params.prefix}:filter_blast"

    publishDir "${params.outdir}/post_assembly/blast", mode: 'copy' 

    input:
        path blast_out
        val min_blast_aln_len
        val min_blast_aln_frac
    output:
        path "${blast_out.baseName}.filtered.txt" optional true
        
    """
    #!/usr/bin/env Rscript

    library(dplyr)

    tb <- read.table("${blast_out}", sep="\t", header=TRUE)
    filtered_tb <- tb %>%
            filter(
                .data[["ALN_LEN"]] > $min_blast_aln_len,
                .data[["ALN_LEN"]] / .data[["QUERY_LEN"]] > $min_blast_aln_frac
            )
    if (nrow(filtered_tb)>0) {
        sorted_filtered_tb <- tb[order(-tb\$"BITSCORE"),]
        write.table(sorted_filtered_tb, file ="${blast_out.baseName}.filtered.txt", row.names=FALSE, sep='\t', quote=FALSE)
    }
    """
}

workflow filter_blast_nucl {
    take:
        blast_out
    emit:
        filtered
    main:
        filtered = filter_blast_nucl_process(blast_out)
}

process filter_blast_nucl_process {
    tag "${params.prefix}:filter_blast_nucl"

    publishDir "${params.outdir}/post_assembly/blast", mode: 'copy' 

    input:
        path blast_out
    output:
        path "${blast_out.baseName}.filtered.txt" optional true
        
    """
    #!/usr/bin/env Rscript

    library(dplyr)

    tb <- read.table("${blast_out}", sep="\t", header=TRUE)
    filtered_tb <- tb %>%
            filter(
                .data[["ALN_LEN"]] > ${params.min_nucl_blast_aln_len},
                .data[["ALN_LEN"]] / .data[["QUERY_LEN"]] > ${params.min_nucl_blast_aln_frac}
            )
    if (nrow(filtered_tb)>0) {
        sorted_filtered_tb <- tb[order(-tb\$"BITSCORE"),]
        write.table(sorted_filtered_tb, file ="${blast_out.baseName}.filtered.txt", row.names=FALSE, sep='\t', quote=FALSE)
    }
    """
}

workflow filter_blast_prot {
    take:
        blast_out
    emit:
        filtered
    main:
        filtered = filter_blast_prot_process(blast_out)
}

process filter_blast_prot_process {
    tag "${params.prefix}:filter_blast_prot"

    publishDir "${params.outdir}/post_assembly/blast", mode: 'copy' 

    input:
        path blast_out
    output:
        path "${blast_out.baseName}.filtered.txt" optional true
        
    """
    #!/usr/bin/env Rscript

    library(dplyr)

    tb <- read.table("${blast_out}", sep="\t", header=TRUE)
    filtered_tb <- tb %>%
            filter(
                .data[["ALN_LEN"]] > ${params.min_prot_blast_aln_len},
                .data[["ALN_LEN"]] / .data[["QUERY_LEN"]] > ${params.min_prot_blast_aln_frac}
            )
    if (nrow(filtered_tb)>0) {
        sorted_filtered_tb <- tb[order(-tb\$"BITSCORE"),]
        write.table(sorted_filtered_tb, file ="${blast_out.baseName}.filtered.txt", row.names=FALSE, sep='\t', quote=FALSE)
    }
    """
}
