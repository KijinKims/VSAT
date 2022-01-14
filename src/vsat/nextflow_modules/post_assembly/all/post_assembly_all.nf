nextflow.enable.dsl=2
include { blast } from '../blast/post_assembly_blast' addParams(tool: 'blastn megablast diamond')
include { filter_blast } from '../../filter/blast/filter_blast'
include { match_taxonomy } from '../../report/blast/report_blast'
include { zoonotic_rank } from '../zoonosis/post_assembly_zoonosis' addParams(tool: 'zoonotic_rank')

workflow {
    Channel.fromPath(params.x).set{contigs}

    blast(contigs)
    filter_blast(blast.out, params.min_blast_aln_len, params.min_blast_aln_frac)
    match_taxonomy(filter_blast.out)

    //zoonotic_rank(contigs)
}