#!/usr/bin/env python3
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
import sys, re
from pathlib import Path
from Bio.SeqRecord import SeqRecord

gff = sys.argv[1]
fasta = sys.argv[2]
print("%s %s : extracting."%(gff,fasta))

cdsfile = Path(fasta).with_suffix('.cds.fasta')
pepfile = Path(fasta).with_suffix('.aa.fasta')

in_seq_handle = open(fasta)
seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
in_seq_handle.close()

new_seq_dict = {}
for s in seq_dict:
    desc = seq_dict[s].description
    m = re.search("^\d+\s+([^:]+)\:(\d+)\-(\d+)$",desc)
#    print("desc is '%s'"%(desc))
    if m:
        chrname = m.group(1)
        new_seq_dict[chrname] = seq_dict[s]
#        print("storing %s as chr"%(chrname))
    else:
#        print('storing old value %s'%(s))
        new_seq_dict[s] = seq_dict[s]

cds_seqs = []
pep_seqs = []
limit_info = {}

# Read the gff
for rec in GFF.parse(gff,limit_info=limit_info):
    faseq = new_seq_dict[rec.id]
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

            if len(seq_CDS) > 0:
                cds_seqs.append(seq_CDS)
                pep_seq = seq_CDS
                pep_seq.seq = seq_CDS.seq.translate()
                pep_seqs.append(pep_seq)

SeqIO.write(cds_seqs, cdsfile, "fasta")
SeqIO.write(pep_seqs, pepfile, "fasta")
