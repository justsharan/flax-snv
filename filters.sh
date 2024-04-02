#!/bin/bash

# Remove INDELs, bases around INDELs, multiallelic sites
bcftools filter -Ou -g 5:indel,other $1 | bcftools view -Ob -M 2 -v snps > remove_indels.bcf

# Hard filters
bcftools filter -Oz -e "INFO/DP>2.5*AVG(INFO/DP) | QUAL<30 | MQ<30 | SP>3" remove_indels.bcf > filtered.vcf.gz

# SnpEff
java -jar $snpEff/snpEff.jar flax filtered.vcf.gz -csvStats stats.csv | gzip > annotated.vcf.gz

# Filter for HIGH or MODERATE variants
java -jar $snpEff/SnpSift.jar filter "(ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE')" annotated.vcf.gz | gzip > high_moderate.vcf.gz
java -jar $snpEff/SnpSift.jar filter "ANN[*].IMPACT = 'HIGH'" high_moderate.vcf.gz | gzip > high.vcf.gz
java -jar $snpEff/SnpSift.jar filter "ANN[*].IMPACT = 'MODERATE'" high_moderate.vcf.gz | gzip > moderate.vcf.gz
