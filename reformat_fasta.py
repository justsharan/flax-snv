from Bio import SeqIO
import sys

# determine path of output file
ext = sys.argv[1].split(".")[-1]
outpath = f"{sys.argv[1][:-len(ext)]}out.{ext}"

# :)
with open(sys.argv[1], "r") as file, open(outpath, "w") as out:
    seqs = SeqIO.parse(file, "fasta")
    SeqIO.write(seqs, out, "fasta")
