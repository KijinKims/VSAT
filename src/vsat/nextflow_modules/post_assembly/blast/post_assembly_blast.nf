nextflow.enable.dsl=2

workflow {
    Channel.fromPath(params.x).set{contigs}
    
    blast(contigs)
}

workflow blast {
    take:
        contigs

    emit:
        blast_outs

    main:

    blast_outs = Channel.empty()

    if (params.tool.contains("blastn")) {
        Channel.fromPath(params.blast_db_dir, type: 'dir').set{blast_db_dir}
        blastn(contigs, blast_db_dir)
        blast_outs.concat(blastn.out)
    }
    
    if (params.tool.contains("megablast")) {
        Channel.fromPath(params.blast_db_dir, type: 'dir').set{blast_db_dir}
        megablast(contigs, blast_db_dir)
        blast_outs.concat(megablast.out)
    }

    if (params.tool.contains("diamond")) {
        Channel.fromPath("${params.diamond_db}/rvdb-prot.dmnd", type: 'dir').set{diamond_db}
        diamond(contigs, diamond_db)
        blast_outs.concat(diamond.out)
    }
}

process blastn {
    tag "${params.prefix}:blastn"

    publishDir "${params.outdir}/post_assembly/blast", mode: 'copy'

    input:
        path contigs
        path blast_db_dir
    output:
        path "${params.prefix}.blastn.txt"
        val "blastn"
    """
    { echo 'QUERY_ID\tREF_ID\tTAX_ID\tREF_TITLE\tPER_IDENT\tQUERY_LEN\tALN_LEN\tMISMATCH\tGAPOPEN\tQUERY_START\tQUERY_END\tREF_START\tREF_END\tEVALUE\tBITSCORE'; \
    blastn -query $contigs -db "${blast_db_dir}/${params.blast_db_name}" -task blastn -evalue ${params.min_evalue} -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads ${params.threads}; \
    } > ${params.prefix}.blastn.txt
    """
}

process megablast {
    tag "${params.prefix}:megablast"

    publishDir "${params.outdir}/post_assembly/blast", mode: 'copy'

    input:
        path contigs
        path blast_db_dir
    output:
        path "${params.prefix}.megablast.txt"
        val "megablast"
    """
    { echo 'QUERY_ID\tREF_ID\tTAX_ID\tREF_TITLE\tPER_IDENT\tQUERY_LEN\tALN_LEN\tMISMATCH\tGAPOPEN\tQUERY_START\tQUERY_END\tREF_START\tREF_END\tEVALUE\tBITSCORE'; \
    blastn -query $contigs -db "${blast_db_dir}/${params.blast_db_name}" -task megablast -evalue ${params.min_evalue} -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads ${params.threads}; \
    } > ${params.prefix}.megablast.txt
    """
}

process diamond {
    tag "${params.prefix}:diamond"

    publishDir "${params.outdir}/post_assembly/blast", mode: 'copy'

    input:
        path contigs
        path diamond_db
    output:
        path "${params.prefix}.diamond.txt"
        val "diamond"
    """
    { echo 'QUERY_ID\tREF_ID\tTAX_ID\tREF_TITLE\tPER_IDENT\tQUERY_LEN\tALN_LEN\tMISMATCH\tGAPOPEN\tQUERY_START\tQUERY_END\tREF_START\tREF_END\tEVALUE\tBITSCORE'; \
    diamond blastx -d ${diamond_db} -q $contigs --max-target-seqs 1 --evalue ${params.min_evalue} --outfmt 6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore; \
    } > ${params.prefix}.diamond.txt
    """
}

