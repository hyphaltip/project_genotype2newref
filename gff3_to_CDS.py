#!/usr/bin/env python3
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
import sys
from pathlib import Path
from Bio.SeqRecord import SeqRecord

gff = sys.argv[1]
fasta = sys.argv[2]
print("%s %s : extracting."%(gff,fasta))

cdsfile = Path(gff).with_suffix('.cds.fasta')
pepfile = Path(gff).with_suffix('.aa.fasta')

in_seq_handle = open(fasta)
seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
in_seq_handle.close()


cds_seqs = []
pep_seqs = []
limit_info = {}

# Read the gff
for rec in GFF.parse(gff,limit_info=limit_info):
    faseq = seq_dict[rec.id]
    # only focus on the CDSs
    #print(rec.features[0])
    for feat in filter(lambda x: x.type == "gene", rec.features):
        # extract the locus tag
        gene_name = feat.id
        for mRNA in feat.sub_features:
            seq_CDS  = SeqRecord("",mRNA.id, "gene=%s"%(gene_name),"")
            seq_CDS.id = mRNA.id
            CDS = filter(lambda x: x.type == "CDS",mRNA.sub_features)
            for cds in sorted(CDS,key=lambda f: f.location.start * f.location.strand):
                #print(cds)
                seq_CDS.seq += cds.extract(faseq).seq
            print(mRNA.id,seq_CDS)
            cds_seqs.append(seq_CDS)
            #pep_seqs.append(cds_record.translate(cds=True))
        break
    break
        #dna_seq = str(feat.extract(faseq).seq)

        # simply print the sequence in fasta format
        #print('>%s\n%s' % (locus_tag, dna_seq))


SeqIO.write(cds_seqs, cdsfile, "fasta")
#SeqIO.write(pep_seqs, pepfile, "fasta")
