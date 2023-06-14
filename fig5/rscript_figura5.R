
library(ggplot2)
library(gggenes)
library(stringr)
library(tidyverse)

setwd("./")

df1 <- read.csv("all_incr.tsv", sep = "\t", header = T)
head(df1)

df1 <- df1 %>%
  mutate(id = case_when(
    endsWith(Sequence.Id, "contig_1") ~ "F15",
    endsWith(Sequence.Id, "contig_2") ~ "F19",
    endsWith(Sequence.Id, "contig_3") ~ "F24",
    endsWith(Sequence.Id, "contig_4") ~ "F28",
    endsWith(Sequence.Id, "contig_5") ~ "F29",
    endsWith(Sequence.Id, "contig_6") ~ "F32",
    endsWith(Sequence.Id, "contig_7") ~ "F34",
    endsWith(Sequence.Id, "contig_8") ~ "F35",
    endsWith(Sequence.Id, "contig_9") ~ "F41",
    endsWith(Sequence.Id, "contig_10") ~ "F43",
    endsWith(Sequence.Id, "contig_11") ~ "F47",
    endsWith(Sequence.Id, "contig_12") ~ "F53",
    endsWith(Sequence.Id, "contig_13") ~ "F54",
    endsWith(Sequence.Id, "contig_14") ~ "F87",
    endsWith(Sequence.Id, "contig_15") ~ "F88"))
head(df1)

df1[grepl("KPC", df1$Gene),]


df1$Strand <- gsub("[+]", "forward", df1$Strand)
df1$Strand <- gsub("[-]", "reverse", df1$Strand)
df1 <- df1 %>%
  mutate(orientation= case_when(
    endsWith(Strand, "reverse") ~ "1",
    endsWith(Strand, "forward") ~ "0"))

dff <- df1[-which(df1$id == "F32"),]

dff <- dff %>% filter(Start > 50000 & Stop < 75000)
head(dff)
dim(dff)
unique(dff$id)
dff$orientation <- dff$orientation %>% replace_na('0')
dff$orientation <- as.numeric(dff$orientation)

cl <- c("hypothetical protein" = "#ccff33", 
        "Integron gene cassette protein" = "#005f73", 
        "sulfonamide-resistant dihydropteroate synthase Sul1" = "#0a9396", 
        "quinolone resistance pentapeptide repeat protein QnrB2" = "#94d2bd", 
        "origin of replication" = "#e9d8a6", 
        "Signal peptide protein" = "#ee9b00", 
        "Dipeptide-binding protein" = "#ca6702", 
        "IS91 family transposase" = "#bb3e03", 
        "quaternary ammonium compound efflux SMR transporter QacE delta 1" = "#ae2012", 
        "trimethoprim-resistant dihydrofolate reductase DfrA25" = "#ffee32",
        "putative aminoglycoside riboswitch / attI site" = "#f94144", 
        "class 1 integron integrase IntI1" = "#f3722c", 
        "Resolvase/invertase-type recombinase catalytic domain-containing protein" = "#f8961e", 
        "IS1182 family ISKpn6 transposase" = "#f9844a", 
        "carbapenem-hydrolyzing class A beta-lactamase KPC-2" = "#f9c74f", 
        "IS21-like element ISKpn7 family helper ATPase IstB" = "#90be6d", 
        "IS21 family ISKpn7 transposase" = "#43aa8b", 
        "Tn3-like element Tn4401 family transposase" = "#4d908e", 
        "Tn3-like element Tn4401 family resolvase TnpR" = "#577590", 
        "Transposase" = "#277da1",
        "Nucleotidyl transferase AbiEii/AbiGii toxin family protein"  = "#8cb369", 
        "TEM family class A beta-lactamase" = "#f4e285", 
        "IS6 family IS15DIV transposase" = "#f4a259", 
        "IS6 family IS15DIV transposase" = "#5b8e7d", 
        "Integrase catalytic region" = "#bc4b51", 
        "broad-spectrum mercury transporter MerE" = "#8cb369", 
        "mercury(II) reductase" = "#f4e285", 
        "organomercurial transporter MerC" = "#f4a259", 
        "IS1 family ISKpn14 transposase ORF A" = "#5b8e7d")



fig <- ggplot(dff, aes(xmin = Start, xmax = Stop, y = id, fill = Product, label=Gene, forward = orientation)) +
  geom_gene_arrow(arrowhead_height = grid::unit(12, "mm"),
                  arrowhead_width = grid::unit(6, "mm"),
                  arrow_body_height = grid::unit(6, "mm")) +
  theme(legend.position="none") +
  geom_gene_label(align = "left") +
  guides(fill = guide_legend(nrow = 8)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.title = element_text(size=16),
        strip.text = element_text(size = 12),
        strip.background =element_rect(fill="white")) +
  facet_wrap(~id, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = cl) 
fig


dir.create("./figuras")
png("./figuras/gggenes_supp_all.png", res = 600, height = 30, width = 45, units = "cm")
fig
dev.off()


s <- df1[which(df1$id == "F19"),]
s <- s %>% filter(Start > 50000 & Stop < 75000)
head(s)
dim(s)
unique(s$id)
s$orientation <- s$orientation %>% replace_na('0')
s$orientation <- as.numeric(s$orientation)

fig2 <- ggplot(s, aes(xmin = Start, xmax = Stop, y = id, fill = Product, label=Gene, forward = orientation)) +
  geom_gene_arrow(arrowhead_height = grid::unit(12, "mm"),
                  arrowhead_width = grid::unit(6, "mm"),
                  arrow_body_height = grid::unit(6, "mm")) +
  geom_gene_label(height = grid::unit(5, "mm"), grow = TRUE, align = "centre") +
  theme(legend.position="bottom") +
  guides(fill = guide_legend(nrow = 8)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 14),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12),
        strip.background =element_rect(fill="white")) +
  facet_wrap(~id, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = cl)
fig2


dir.create("./figuras")
png("./figuras/incR.png", res = 600, height = 10, width = 45, units = "cm")
fig2
dev.off()
