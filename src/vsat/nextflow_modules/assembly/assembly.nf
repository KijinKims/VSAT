nextflow.enable.dsl=2

workflow {
    emit:
        contigs

    main:

        if (params.platform == 'illumina') {
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
            contigs = assembly_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
            contigs = assembly_nanopore(fastx)
        } else { //hybrid
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
            Channel.fromPath(params.y).set{fastx}
            contigs = assembly_hybrid(fastq_pair, fastx)
        }
}

workflow assembly_illumina {
    take:
        fastq_pair
    emit:
        contigs
    main:
        
        contigs = spades(fastq_pair)
}

workflow assembly_nanopore {
    take:
        fastx
    emit:
        contigs
    main:
        
        if (params.tool == "canu") {
            contigs = canu(fastx)
        } else {
            contigs = megahit(fastx)
        }
}

workflow assembly_hybrid {
    take:
        fastq_pair
        fastx
    emit:
        contigs
    main:
        
        contigs = unicycler(fastq_pair, fastx)

}

process spades {
    tag "${params.prefix}:spades"

    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}
    
    input:
        tuple path(pe1), path(pe2)
    output:
        path "${params.prefix}/contigs.fasta"
    script:

    if (params.spades_meta)
        """
        spades.py --pe1-1 $pe1 --pe1-2 $pe2 -o ${params.prefix} --meta -t ${params.threads}
        """
    else
        """
        spades.py --pe1-1 $pe1 --pe1-2 $pe2 -o ${params.prefix} -t ${params.threads}
        """
}

process megahit {
    tag "${params.prefix}:megahit"
    
    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}

    input:
        path fastx
    output:
        path "${params.prefix}/${params.prefix}.contigs.fa"
    script:

    """
    megahit -r $fastx -o $params.prefix --out-prefix $params.prefix -t ${params.threads}
    """
}

process canu {
    tag "${params.prefix}:canu"
    
    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}

    input:
        path fastx
    output:
        path "${params.prefix}/${params.prefix}.contigs.fasta"
    script:

    """
    canu -p $params.prefix -d $params.prefix -nanopore $fastx $params.canu_options
    """
}

process unicycler {
    tag "${params.prefix}:unicycler"
    
    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}

    input:
        tuple path(pe1), path(pe2)
        path fastx
    output:
        path "${params.prefix}/${params.prefix}.contigs.fasta"
    script:
    //TO DO
    """
    unicycler 
    """
}