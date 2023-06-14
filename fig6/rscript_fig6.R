library(ggplot2)
library(gggenes)
library(stringr)
library(tidyverse)

setwd("/mnt/raid2tb/bubble/tesis/capitulo1/fig6")

df1 <- read.csv("all_samp.tsv", sep = "\t", header = T)
head(df1)

df1 <- df1 %>%
  mutate(id = case_when(
    endsWith(Sequence.Id, "contig_1") ~ "F19",
    endsWith(Sequence.Id, "contig_2") ~ "F24",
    endsWith(Sequence.Id, "contig_3") ~ "F41"))
head(df1)

unique(df1$id)

df1[grepl("NDM", df1$Gene),]


df1$Strand <- gsub("[+]", "forward", df1$Strand)
df1$Strand <- gsub("[-]", "reverse", df1$Strand)
df1 <- df1 %>%
  mutate(orientation= case_when(
    endsWith(Strand, "reverse") ~ "1",
    endsWith(Strand, "forward") ~ "0"))

dff <- df1[which(df1$id == "F24"),]

dff <- dff %>% filter(Start > 110000 & Stop < 130000)
head(dff)
dim(dff)
unique(dff$id)
dff$orientation <- dff$orientation %>% replace_na('0')
dff$orientation <- as.numeric(dff$orientation)

unique(dff$Product)

cl <- c("Resolvase/invertase-type recombinase catalytic domain-containing protein" = "#001219",
        "class 1 integron integrase IntI1" = "#005f73",
        "trimethoprim-resistant dihydrofolate reductase DfrA12" = "#0a9396",
        "quaternary ammonium compound efflux SMR transporter QacE delta 1" = "#94d2bd",
        "IS91 family transposase" = "#e9d8a6",
        "N-(5'-phosphoribosyl)anthranilate isomerase" = "#ee9b00",
        "subclass B1 metallo-beta-lactamase NDM-1" = "#ca6702",
        "chloramphenicol efflux MFS transporter CmlA5" = "#bb3e03",
        "DUF3330 domain-containing protein" = "#ae2012",
        "putative aminoglycoside riboswitch / attI site" = "#9b2226",
        "ANT(3'')-Ia family aminoglycoside nucleotidyltransferase AadA2" = "#f94144",
        "sulfonamide-resistant dihydropteroate synthase Sul1" = "#f3722c",
        "hypothetical protein" = "#f8961e",
        "bleomycin binding protein Ble-MBL" = "#f9844a",
        "quinolone resistance pentapeptide repeat protein QnrA1" = "#f9c74f",
        "IS110 family IS5075 transposase" = "#90be6d")



fig <- ggplot(dff, aes(xmin = Start, xmax = Stop, y = id, fill = Product, label=Gene, forward = orientation)) +
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
fig


dir.create("./figuras")

dpp <- df1[which(df1$id == "F41"),]

dpp <- dpp %>% filter(Start > 5000 & Stop < 25000)
head(dpp)
dim(dpp)
unique(dpp$id)
dpp$orientation <- dpp$orientation %>% replace_na('0')
dpp$orientation <- as.numeric(dpp$orientation)

unique(dpp$Product)

fig2 <- ggplot(dpp, aes(xmin = Start, xmax = Stop, y = id, fill = Product, label=Gene, forward = orientation)) +
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

library(ggpubr)

a <- ggarrange(fig, fig2, nrow = 2, common.legend = T, legend = "bottom") 
a

dir.create("./figuras")
png("./figuras/incA.png", res = 600, height = 15, width = 45, units = "cm")
a
dev.off()
