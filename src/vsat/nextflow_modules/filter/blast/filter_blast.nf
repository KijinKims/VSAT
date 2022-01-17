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

workflow filter_blast {
    take:
        blast_out
    emit:
        filtered
    main:
        filter_blast_by_algo(blast_out)
}

process filter_blast_by_algo {
    tag "${params.prefix}:filter_blast_by_algo"

    publishDir "${params.outdir}/post_assembly/blast", mode: 'copy' 

    input:
        path blast_out
    output:
        path "${blast_out.baseName}.filtered.txt" optional true
    script:
    if (blast_out.BaseName.Extension == "diamond")
        min_blast_aln_len=${params.min_prot_blast_aln_len}
        min_blast_aln_frac=${params.min_prot_blast_aln_frac}
    else
        min_blast_aln_len=${params.min_nucl_blast_aln_len}
        min_blast_aln_frac=${params.min_nucl_blast_aln_frac}
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