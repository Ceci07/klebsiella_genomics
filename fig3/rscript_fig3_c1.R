library(ggplot2)
library(ggsignif)
library(dplyr)

setwd("./")

df <- read.table("tabla_resumen.tsv", sep = "\t", header = T)
head(df)
colnames(df) <- c("Muestra", "Largo_total", "Contigs", "contig_largo", "ST", "Plataforma")
df$Largo_total <- gsub("[.]", "", df$Largo_total)
df$Largo_total <- as.numeric(df$Largo_total)
df$contig_largo <- gsub("[.]", "", df$contig_largo)
df$contig_largo <- as.numeric(df$contig_largo)
df$Contigs <- as.numeric(df$Contigs)


cl <- c("#003049", "#d62828", "#f77f00", "#fcbf49", "#eae2b7", 
        "#d88c9a", "#f2d0a9", "#f1e3d3", "#99c1b9", "#8e7dbe",
        "#227c9d", "#17c3b2", "#ffcb77", "#fef9ef", "#fe6d73")

dt <- select(df, Muestra, Plataforma, Contigs)
dim(dt)
max(dt$Contigs)
colnames(dt)[1] <- c("Muestra")
head(dt)


p <- ggplot(dt, aes(x=Plataforma, y=Contigs)) + 
  geom_boxplot(alpha = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(size = 10)) +
  geom_signif(comparisons = list(c("ONT", "ONT + Illumina")), map_signif_level = T, y_position = c(120)) +
  geom_signif(comparisons = list(c("ONT", "Illumina")), map_signif_level = T, y_position = c(165)) +
  geom_jitter(size = 3, pch=21, aes(fill = Muestra), alpha = 1) +
  scale_fill_manual(values = cl) +ylab("Número de contigs") +
  scale_y_continuous(breaks=seq(0, 170, 20)) +
  guides(fill=guide_legend(nrow=2))
p



dt <- dt %>% group_by(Plataforma) %>%
  summarise(prom = mean(Contigs), ds = sd(Contigs))

shapiro.test(df$Contigs) #not normal distributed


l <- select(df, Muestra, Plataforma, Largo_total)
l <- l %>% group_by(Plataforma) %>%
  summarise(prom = mean(Largo_total), ds = sd(Largo_total))


q <- ggplot(df, aes(x=Plataforma, y=Largo_total)) + 
  geom_boxplot(alpha = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(size = 10)) +
  geom_signif(comparisons = list(c("ONT", "ONT + Illumina")), map_signif_level = T, y_position = c(5800000)) +
  geom_signif(comparisons = list(c("ONT", "Illumina")), map_signif_level = T, y_position = c(5750000)) +
  geom_jitter(size = 3, pch=21, aes(fill = Muestra), alpha = 1) +
  scale_fill_manual(values = cl) + ylab("Largo total (pb)")+
  scale_y_continuous(limits=c(5100000, 6100000))+
  guides(fill=guide_legend(nrow=2))
q



options(scipen = 999)
head(df)

df <- as.data.frame(df)
s <- ggplot(df, aes(x=Plataforma, y=contig_largo)) + 
  geom_boxplot(alpha = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(size = 10)) +
  geom_signif(comparisons = list(c("ONT", "ONT + Illumina")), map_signif_level = T, y_position = c(5600000)) +
  geom_signif(comparisons = list(c("ONT", "Illumina")), map_signif_level = T, y_position = c(5300000)) +
  geom_jitter(size = 3, pch=21, aes(fill = Muestra), alpha = 1) +
  scale_fill_manual(values = cl) +
  ylab("Contig mayor (pb)")+
  scale_y_continuous(limits=c(200000, 6100000)) +
  guides(fill=guide_legend(nrow=2))
s

cll <- select(df, Muestra, Plataforma, contig_largo)
cll <- cll %>% group_by(Plataforma) %>%
  summarise(prom = mean(contig_largo), ds = sd(contig_largo))

library(ggpubr)
fig <- ggarrange(p, q, s, ncol = 3, common.legend = T, legend = "bottom", labels = c("A", "B", "C"))
fig

png('fig3_cap1.png', res = 600, height = 10, width = 25, units = 'cm')
fig
dev.off()
