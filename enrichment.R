library(tidyverse)

samples <- c("E150016690_L01_122", "E150016690_L01_123", "E150016690_L01_124", "E150016690_L01_29", "E150016690_L01_30", "E150016690_L01_35", "E150016690_L01_39", "E150016690_L01_55", "E150016690_L01_94", "E150016690_L01_96", "E150016690_L01_97", "E150025109_L01_96", "E150025345_L01_104", "V350170890_L03_100", "V350170890_L03_101", "V350170890_L03_102", "V350170890_L03_103", "V350170890_L03_104", "V350170890_L03_121", "V350170890_L03_122", "V350170890_L03_123", "V350170890_L03_124", "V350170890_L03_125", "V350170890_L03_126", "V350170890_L03_127", "V350170890_L03_88", "V350170890_L03_89", "V350170890_L03_90", "V350170890_L03_91", "V350170890_L03_92", "V350170890_L03_93", "V350170890_L03_94", "V350170890_L03_95", "V350170890_L03_96", "V350170890_L03_97", "V350170890_L03_98", "V350170890_L03_99", "V350170890_L04_114", "V350170890_L04_115", "V350170890_L04_116", "V350170890_L04_117", "V350170890_L04_128", "V350170890_L04_25", "V350170890_L04_26", "V350170890_L04_28", "V350170890_L04_29", "V350170890_L04_30", "V350170890_L04_32", "V350170890_L04_33", "V350170890_L04_34", "V350170890_L04_35", "V350170890_L04_36", "V350170890_L04_37", "V350170890_L04_38", "V350170890_L04_39", "V350170890_L04_49", "V350170890_L04_50", "V350170890_L04_51", "V350170890_L04_52", "V350170890_L04_53", "V350170890_L04_55")

tsv.cols <- c("Chrom", "Pos", "Gene", "Transcript.ID", "Ref", "Alt", "Effect", "HGVS.C", "HGVS.P", "cDNA.Pos", "AA.Pos", samples)

ortholog.cols <- c("Transcript.ID", "Gene", "Transcript.Name", "Protein.Name", "PFAM", "Panther", "KOG", "KEGG.ec", "KEGG.Orthology", "Gene.Ontology", "TAIR10.Name", "TAIR10.Symbol", "TAIR10.Defline", "Rice.Name", "Rice.Symbol", "Rice.Defline")

# Ortholog table from Phytozome
orthologs <- read_tsv("Lusitatissimum_200_annotation_info.txt.gz", col_names = ortholog.cols) |>
    mutate(Transcript.ID = paste("PAC:", Transcript.ID, sep=""))

high.tsv <- read_tsv("out/high.tsv", skip = 1, col_names = tsv.cols) |>
    mutate(Ref = paste(Ref, ">", Alt, sep="")) |>
    rename(Change = Ref) |>
    select(-Alt)

high.tsv |>
    inner_join(orthologs, by=c("Gene", "Transcript.ID")) |>
    mutate(across(all_of(samples), function(s) str_match(s, "^(\\d/\\d):")[,2])) |>
    write_csv("annotated.csv")
