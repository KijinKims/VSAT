nextflow.enable.dsl=2
include { qc_illumina; qc_nanopore; qc_hybrid } from '../qc/qc'
include { filter_reads_illumina; filter_reads_nanopore; filter_reads_hybrid } from '../filter/reads/filter_reads'
include { filter_host_illumina; filter_host_nanopore; filter_host_hybrid } from '../filter/host/filter_host'
include { map_illumina; map_nanopore } from '../map/map'
include { filter_map } from '../filter/map/filter_map'
include { kmer_illumina; kmer_nanopore } from '../kmer/kmer' addParams(tool: 'kraken2 kaiju')
include { assembly } from '../assembly/assembly'
include { polish } from '../polish/polish' addParams(tool: 'racon medaka')
include { filter_contigs } from '../filter/contigs/filter_contigs'
include { blast } from '../post_assembly/blast/post_assembly_blast' addParams(tool: 'blastn megablast diamond')
include { filter_blast } from '../filter/blast/filter_blast'
include { match_taxonomy } from '../report/blast/report_blast'
include { zoonotic_rank } from '../post_assembly/zoonosis/post_assembly_zoonosis' addParams(tool: 'zoonotic_rank')

workflow {

    main:

    if (params.platform == 'illumina') {
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
        }

    if (params.platform == 'illumina') {
        qc_illumina(fastq_pair)
        filter_reads_illumina(fastq_pair)

        if (params.host_genome != null) {
            Channel.fromPath(params.host_genome).set{host_genome}
            filter_completed = filter_host_illumina(filter_reads_illumina.out, host_genome)
        } else {
            filter_completed = filter_reads_illumina.out
        }

        map_illumina(filter_completed)
        filter_map(map_illumina.out)

        kmer_illumina(filter_completed)

            } else if (params.platform == 'nanopore') {
        
        qc_nanopore(fastq_pair)
        filter_reads_nanopore(fastq_pair)

        if (params.host_genome != null) {
            Channel.fromPath(params.host_genome).set{host_genome}
            filter_completed = filter_host_nanopore(filter_host_nanopore.out, host_genome)
        } else {
            filter_completed = filter_host_nanopore.out
        }

        map_nanopore(filter_completed)
        filter_map(map_illumina.out)

        kmer_nanopore(filter_completed)
        }

    assembly(filter_completed)
    filter_contigs(assembly.out)

    blast(filter_contigs.out)
    filter_blast(blast.out, params.min_blast_aln_len, params.min_blast_aln_frac)
    match_taxonomy(filter_blast.out)

    zoonotic_rank(assembly.out)
}