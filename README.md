# VSAT

## QuickStart

Requirement: Docker, conda(or miniconda)

```bash
git clone https://github.com/KijinKims/VSAT.git
cd VSAT
conda env create -f environment.yml
conda activate vsat
pip install vsat-1.0.tar.gz
```

```bash
sh utils/docker_pull_list_of_images.sh docker_images.list
```

```bash
export VSAT_DB=/location/you/want/to/download/databases/for/vsat
```

Default path is `$HOME/VSAT_DB`

```bash
vsat database download --db kraken2-viral kaiju-viruses refseq-viral rvdb-prot taxdump accession2taxid
cat $VSAT_DB/refseq-viral/*fna > $VSAT_DB/refseq-viral/refseq_viral.fna
vsat database build blast --in $VSAT_DB/refseq-viral/refseq_viral.fna --outdir $VSAT_DB/blast/refseq-viral --dbname refseq-viral --title refseq-viral --parse_seqids --dbtype nucl
vsat database build diamond --in $VSAT_DB/rvdb-prot/U-RVDBv23.0-prot.fasta --dbname rvdb-prot --outdir $VSAT_DB/diamond/
vsat database build taxnomizr --outdir $VSAT_DB/taxonomizr
```

NOTE: The processes to build databases require comparable amount of memory. If your computer is not equipped with that amount of memory, the process will be killed. In that case, you need to look for pre-built database from another source.

```bash
mkdir tutorial && cd tutorial
fasterq-dump SRR13439799
esearch -db nucleotide -query "KJ942813.1" | efetch -format fasta > KJ942813_influenza_A.fasta 
```

```bash
#qc
vsat qc --platform nanopore -x SRR13439799_Sequecing_of_PR8_H1N1_culture_medium_1.fastq --prefix SRR13439799
#filter reads
vsat filter reads --platform nanopore -x SRR13439799_Sequecing_of_PR8_H1N1_culture_medium_1.fastq --prefix SRR13439799
#map
vsat map --platform nanopore -x SRR13439799/filter/SRR13439799.filtered.fastq --prefix SRR13439799 --ref KJ942813_influenza_A.fasta
#filter map
vsat filter map -x SRR13439799/map/SRR13439799.map.tsv --prefix SRR13439799
#kmer
vsat kmer --platform nanopore -x SRR13439799/filter/SRR13439799.filtered.fastq --prefix SRR13439799 
#assembly
vsat assembly --platform nanopore -x SRR13439799/filter/SRR13439799.filtered.fastq --prefix SRR13439799 --tool megahit
#filter contigs
vsat filter contigs -x SRR13439799/assembly/SRR13439799.contigs.fasta --prefix SRR13439799
#polish
vsat polish -x SRR13439799/assembly/SRR13439799.filtered_contigs.fasta --prefix SRR13439799 --reads SRR13439799_Sequecing_of_PR8_H1N1_culture_medium_1.fastq
#blast
vsat post_assembly blast -x SRR13439799/polish/SRR13439799.polished_contigs.fasta --prefix SRR13439799 --blast_db_dir $VSAT_DB/blast/refseq-viral --blast_db_name refseq-viral --diamond_db_dir $VSAT_DB/diamond --diamond_db_name rvdb-prot
#filter blast
vsat filter blast -x SRR13439799/post_assembly/blast/SRR13439799.megablast.txt --prefix SRR13439799
vsat filter blast -x SRR13439799/post_assembly/blast/SRR13439799.diamond.txt --prefix SRR13439799
#report blast
vsat report blast -x SRR13439799/post_assembly/blast/SRR13439799.megablast.filtered.txt --prefix SRR13439799 --taxonomizr_db $VSAT_DB/taxonomizr/accessionTaxa.sql
#zoonosis
vsat post_assembly zoonosis -x SRR13439799/polish/SRR13439799.polished_contigs.fasta --prefix SRR13439799
```



```bash
#post assembly all 
vsat post_assembly all -x SRR13439799/polish/SRR13439799.polished_contigs.fasta --prefix SRR13439799 --blast_db_dir $VSAT_DB/blast/refseq-viral --blast_db_name refseq-viral --diamond_db_dir $VSAT_DB/diamond --diamond_db_name rvdb-prot --taxonomizr_db $VSAT_DB/taxonomizr/accessionTaxa.sql
#end to end
vsat end_to_end --platform nanopore -x SRR13439799_Sequecing_of_PR8_H1N1_culture_medium_1.fastq --prefix SRR13439799
```

```bash
vsat consensus --platform nanopore -x SRR13439799_Sequecing_of_PR8_H1N1_culture_medium_1.fastq --prefix SRR13439799 --ref KJ942813_influenza_A.fasta
```