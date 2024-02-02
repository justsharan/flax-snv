install.packages('vcfR')

library(vcfR)

vcf_file <- 'out/vcf/SRR1610817.vcf.gz'
dna_file <- 'Phytozome/PhytozomeV10/Lusitatissimum/assembly/Lusitatissimum_200_BGIv1.0.fa'
gff_file <- 'Phytozome/PhytozomeV10/Lusitatissimum/annotation/Lusitatissimum_200_v1.0.gene_exons.gff3'

vcf <- read.vcfR(vcf_file, verbose=FALSE)
dna <- ape::read.dna(dna_file, format='fasta')
gff <- read.table(gff_file, sep='\t', quote='')

unique(gff[,1])

chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)

plot(chrom)

chromoqc(chrom, dp.alpha=20)
