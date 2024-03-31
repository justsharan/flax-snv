library(tidyverse)
library(CMplot)
library(patchwork)

vcf.cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", 1:61)

# VCF file after hard filtering
hard.filtered.variants <- read_tsv("filtered.vcf.gz", comment = "#", col_names = vcf.cols)
# VCF file for high/moderate variants
high.moderate.variants <- read_tsv("high_moderate.vcf.gz", comment = "#", col_names = vcf.cols)

## PLOTTING

hard.filtered.snps <- hard.filtered.variants |>
    mutate(Type = paste(REF, ">", ALT, sep = "")) |>
    count(Type) |>
    mutate(Prop = n / sum(n)) |>
    ggplot(aes(x = reorder(Type, -Prop), y = Prop)) +
        geom_bar(stat = "identity") +
        labs(x = "SNP Type", y = "Relative Frequency") +
        ggtitle("All Filtered Variants") +
        scale_y_continuous(labels = scales::percent) +
        theme_classic()

high.moderate.snps <- high.moderate.variants |>
    mutate(Type = paste(REF, ">", ALT, sep = "")) |>
    count(Type) |>
    mutate(Prop = n / sum(n)) |>
    ggplot(aes(x = reorder(Type, -Prop), y = Prop)) +
        geom_bar(stat = "identity") +
        labs(x = "SNP Type", y = "Relative Frequency") +
        ggtitle("High and Moderate Impact Variants") +
        scale_y_continuous(labels = scales::percent) +
        theme_classic()

# Render bar plots
hard.filtered.snps + high.moderate.snps
ggsave("snp_type.jpg", hard.filtered.snps + high.moderate.snps, width = 9, height = 6)

# Generate SNP Density Plot
hard.filtered.variants |>
    select(CHROM, POS) |>
    mutate(SNP = paste("SNP", row_number(), sep = ""), CHROM = substring(CHROM, 3)) |>
    relocate(SNP, CHROM, POS) |>
    CMplot(plot.type = "d", chr.den.col = c("darkgreen", "yellow", "red"),
        file.name = "density_plot", file = "jpg")
