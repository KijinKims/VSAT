nextflow.enable.dsl=2

workflow {
    main:

        if (params.platform == 'illumina') {
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
            kmer_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
            kmer_nanopore(fastx)
        } else { //hybrid
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
            Channel.fromPath(params.y).set{fastx}
            kmer_hybrid(fastq_pair, fastx)
        }
}

workflow kmer_illumina {
    take:
        fastq_pair

    main:
        tool = params.tool.tokenize()

        krona_inputs = Channel.empty()

        if ( tool.contains('kraken2') ) {
            Channel.fromPath(params.kraken2_db, type: 'dir').set{kraken2_db}
            kraken2_pair(fastq_pair, kraken2_db)
            kraken_krona_input = kraken2krona(kraken2_pair.out)
            krona_inputs = krona_inputs.concat(kraken_krona_input)

            kreport2list(kraken2_pair.out)
            metacomp(kreport2list.out)
        }
        
        if ( tool.contains('kaiju') ) {
            Channel.fromPath("${params.kaiju_fmi}").set{kaiju_fmi}
            Channel.fromPath("${params.kaiju_nodes}").set{kaiju_nodes}
            Channel.fromPath("${params.kaiju_names}").set{kaiju_names}
            kaiju_pair(fastq_pair, kaiju_fmi, kaiju_nodes)
            kaiju_krona_input = kaiju2krona(kaiju_pair.out, kaiju_nodes, kaiju_names)
            krona_inputs = krona_inputs.concat(kaiju_krona_input[0])
        }

        krona(krona_inputs)
}

workflow kmer_nanopore {
    take:
        fastx

    main:
        tool = params.tool.tokenize()
        
        krona_inputs = Channel.empty()

        if ( tool.contains('kraken2') ) {
            Channel.fromPath(params.kraken2_db, type: 'dir').set{kraken2_db}
            kraken2_single(fastx, kraken2_db)
            kraken_krona_input = kraken2krona(kraken2_single.out)
            krona_inputs = krona_inputs.concat(kraken_krona_input)

            kreport2list(kraken2_single.out)
            metacomp(kreport2list.out)
        }
        
        if ( tool.contains('kaiju') ) {
            Channel.fromPath("${params.kaiju_fmi}").set{kaiju_fmi}
            Channel.fromPath("${params.kaiju_nodes}").set{kaiju_nodes}
            Channel.fromPath("${params.kaiju_names}").set{kaiju_names}
            kaiju_single(fastx, kaiju_fmi, kaiju_nodes)
            kaiju_krona_input = kaiju2krona(kaiju_single.out, kaiju_nodes, kaiju_names)
            krona_inputs = krona_inputs.concat(kaiju_krona_input[0])
        }

        krona(krona_inputs)
}

process kraken2_pair {
    tag "${params.prefix}:kraken2_pair"

    input:
        tuple path(pe1), path(pe2)
        path kraken2_db
    output:
        path "${params.prefix}.kraken_report.csv"
    script:

    """
    kraken2 --db ${kraken2_db} \
        --report ${params.prefix}.kraken_report.csv \
        --paired --threads ${params.threads} --confidence ${params.kraken2_confidence_threshold}\
        $pe1 $pe2
    """
}

process kraken2_single {
    tag "${params.prefix}:kraken2_single"

    input:
        path single
        path kraken2_db
    output:
        path "${params.prefix}.kraken_report.csv"

    """
    kraken2 --db ${kraken2_db} \
        --report ${params.prefix}.kraken_report.csv \
        --threads ${params.threads} --confidence ${params.kraken2_confidence_threshold}\
        $single
    """
}

process kraken2krona {
    tag "${params.prefix}:kraken2krona"

    input:
        path kraken_report
    output:
        tuple path("${params.prefix}.kraken"), val("kraken")
    """
    kreport2krona.py -r $kraken_report -o ${params.prefix}.kraken
    """
}

process kreport2list {
    tag "${params.prefix}:kreport2list"

    input:
        path kraken_report
    output:
        path "${params.prefix}.list"

    """
    perl ~/convert_krakenRep2list.pl < $kraken_report > ${params.prefix}.list
    """
}

process metacomp {
    tag "${params.prefix}:metacomp"
    
    publishDir "${params.outdir}/kmer", mode: 'copy'

    input:
        path kraken_list
    output:
        path "${params.prefix}_family.svg"
        path "${params.prefix}_genus.svg"
        path "${params.prefix}_species.svg"
    """
    Rscript ~/generate_heatmap.R ${kraken_list} ${params.prefix}
    """
}

process kaiju_pair {
    tag "${params.prefix}:kaiju_pair"

    input:
        tuple path(pe1), path(pe2)
        path kaiju_fmi
        path kaiju_nodes
    output:
        path "${params.prefix}.kaiju.out"
    """
    kaiju -f ${kaiju_fmi} \
        -t ${kaiju_nodes} \
        -i $pe1 -j $pe2 \
        -o ${params.prefix}.kaiju.out \
        -a ${params.kaiju_mode} -E ${params.kaiju_confidence_threshold} -z ${params.threads} -v
    """
}

process kaiju_single {
    tag "${params.prefix}:kaiju_single"

    input:
        path single
        path kaiju_fmi
        path kaiju_nodes
    output:
        path "${params.prefix}.kaiju.out"

    """
    kaiju -f ${kaiju_fmi} \
        -t ${kaiju_nodes} \
        -i $single\
        -o ${params.prefix}.kaiju.out \
        -a ${params.kaiju_mode} -E ${params.kaiju_confidence_threshold} -z ${params.threads} -v
    """
}

process kaiju2krona {
    tag "${params.prefix}:kaiju2krona"
    publishDir "${params.outdir}/kmer", mode: 'copy', pattern: "${params.prefix}.kaiju_summary.tsv"

    input:
        path kraken_report
        path kaiju_nodes
        path kaiju_names
    output:
        tuple path("${params.prefix}.kaiju"), val("kaiju")
        path "${params.prefix}.kaiju_summary.tsv"
    """
    kaiju2table -t ${kaiju_nodes} -n ${kaiju_names} -r species -e ${params.prefix}.kaiju.out -o ${params.prefix}.kaiju_summary.tsv 
    kaiju2krona -t ${kaiju_nodes} -n ${kaiju_names} -i ${params.prefix}.kaiju.out -o ${params.prefix}.kaiju
    """

}

process krona {
    tag "${params.prefix}:krona"
    publishDir "${params.outdir}/kmer", mode: 'copy'

    input:
        tuple path(krona_input), val(tool)
    output:
        path "${params.prefix}_${tool}.html"
    """
    ktImportText $krona_input -o ${params.prefix}_${tool}.html
    """
}