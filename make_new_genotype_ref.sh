#!/usr/bin/bash
#SBATCH -p short --mem 24gb

module load gatk
module load samtools
module load bcftools
module load GAL
module unload miniconda2
module load miniconda3
source activate bcbio

# REF GENOME
RELEASE=47
PREFIX=AfumigatusA1163
DB=FungiDB-${RELEASE}_${PREFIX}_Genome.fasta
DBURL="https://fungidb.org/common/downloads/release-${RELEASE}/$PREFIX/fasta/data/$DB"

GFF=$(basename $DB _Genome.fasta)".gff"
GFFURL=https://fungidb.org/common/downloads/release-${RELEASE}/$PREFIX/gff/data/$GFF

# ISOLATE
STRAIN=W7
VCF=W7.gatk3.comb_selected.vcf.gz

if [ ! -f $DB ]; then
	curl -O -L $DBURL
fi
if [ ! -f $GFF ]; then
	curl -O -L $GFFURL
fi

if [ ! -s $(basename $DB .fasta)".dict" ]; then
	gatk CreateSequenceDictionary -R $DB
fi
if [ ! -f $DB.fai ]; then
	samtools faidx $DB
fi
if [ ! -s $VCF.tbi ]; then
	taxbix $VCF
fi
if [ ! -s $STRAIN.new_ref_from_${PREFIX}.fasta ]; then
	gatk FastaAlternateReferenceMaker -R $DB -O $STRAIN.new_ref_from_${PREFIX}.fasta -V $VCF
fi

./gff3_to_CDS.py FungiDB-47_AfumigatusA1163.gff W7.new_ref_from_AfumigatusA1163.fasta
