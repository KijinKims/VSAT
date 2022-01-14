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
vsat database download --db kraken2-viral kaiju-viruses rvdb rvdb-prot taxdump accession2taxid
vsat database build blast --in $VSAT_DB/rvdb/C-RVDBv23.0.fasta --outdir $VSAT_DB/blast/rvdb --dbname rvdb --title rvdb_blast --parse_seqids --dbtype nucl
vsat database build diamond --in $VSAT_DB/rvdb-prot/U-RVDBv23.0-prot.fasta
vsat database build taxnomizr
```

NOTE: The processes to build databases require comparable amount of memory. If your computer doesn't provide that amount of memory, the process will be killed. In that case, you need to look for pre-built database from another source.

```bash
mkdir tutorial && cd tutorial
SRX9853219
KJ942813.1 
```

```bash
#qc
vsat qc --platform nanopore -x SRR13439799_Sequecing_of_PR8_H1N1_culture_medium_1.fastq --prefix SRR13439799
#filter reads
vsat filter reads --platform nanopore -x SRR13439799_Sequecing_of_PR8_H1N1_culture_medium_1.fastq --prefix SRR13439799
#map
vsat map --platform nanopore -x SRR13439799/filter/SRR13439799.filtered.fastq --prefix SRR13439799 --ref KJ942813_influenza_A.fasta
#filter map
vsat filter map -x SRR13439799/map/SRR13439799.map.tsv
#kmer
vsat kmer --platform nanopore -x SRR13439799/filter/SRR13439799.filtered.fastq --prefix SRR13439799 
#assembly
vsat assembly --platform nanopore -x SRR13439799/filter/SRR13439799.filtered.fastq --prefix SRR13439799 --tool megahit
#filter contigs
vsat filter contigs -x SRR13439799/assembly/SRR13439799.contigs.fasta --prefix SRR13439799
#polish
vsat polish -x SRR13439799/assembly/SRR13439799.filtered_contigs.fasta --prefix SRR13439799 --reads SRR13439799_Sequecing_of_PR8_H1N1_culture_medium_1.fastq
#blast
vsat post_assembly blast -x SRR13439799/polish/SRR13439799.polished_contigs.fasta --prefix SRR13439799
#filter blast
vsat filter blast -x SRR13439799/post_assembly/blast/SRR13439799.megablast.txt --prefix SRR13439799
vsat filter blast -x SRR13439799/post_assembly/blast/SRR13439799.diamond.txt --prefix SRR13439799
#report blast
vsat report blast -x SRR13439799/post_assembly/blast/SRR13439799.megablast.filtered.txt --prefix SRR13439799
#zoonosis
vsat post_assembly zoonosis -x SRR13439799/assembly/polished_contigs.fasta --prefix SRR13439799
```

