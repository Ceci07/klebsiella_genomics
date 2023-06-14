
library(dplyr)
library(ggplot2)

setwd("./")

df <- read.table("pavian_out_klebs_count.tsv", sep = "\t", header = T)
head(df)

colnames(df)[5] <- c("F15")
colnames(df)[6] <- c("F19")
colnames(df)[7] <- c("F24")
colnames(df)[8] <- c("F28")
colnames(df)[9] <- c("F29")
colnames(df)[10] <- c("F32")
colnames(df)[11] <- c("F34")
colnames(df)[12] <- c("F35")
colnames(df)[13] <- c("F41")
colnames(df)[14] <- c("F43")
colnames(df)[15] <- c("F47")
colnames(df)[16] <- c("F53")
colnames(df)[17] <- c("F54")
colnames(df)[18] <- c("F87")
colnames(df)[19] <- c("F88")

df$lineage <- NULL
df$taxRank <- NULL
df$taxID <- NULL
df$Max <- NULL

library(data.table)
long <- melt(setDT(df), id.vars = c("name"), variable.name = "sample")
head(long)

long <- as.data.frame(long)

target <- c("Klebsiella pneumoniae")
int <- filter(long, name %in% target)

others <- anti_join(long, int, by = "name") 
  
others$name <- NULL
others$name <- c("others")
others$value <- as.numeric(others$value)
others[is.na(others)] <- 0

others <- others %>% group_by(name, sample) %>%
  summarise(value2 = sum(value))
colnames(others)[3] <- c("value")

sam <- as.data.frame(rbind(int, others))

F15 <- sam[which(sam$sample == "F15"),]
F15$fraction = F15$value / sum(F15$value)
F15$ymax = cumsum(F15$fraction)
F15$ymin = c(0, head(F15$ymax, n=-1))

F19 <- sam[which(sam$sample == "F19"),]
F19$fraction = F19$value / sum(F19$value)
F19$ymax = cumsum(F19$fraction)
F19$ymin = c(0, head(F19$ymax, n=-1))

F24 <- sam[which(sam$sample == "F24"),]
F24$fraction = F24$value / sum(F24$value)
F24$ymax = cumsum(F24$fraction)
F24$ymin = c(0, head(F24$ymax, n=-1))

F28 <- sam[which(sam$sample == "F28"),]
F28$fraction = F28$value / sum(F28$value)
F28$ymax = cumsum(F28$fraction)
F28$ymin = c(0, head(F28$ymax, n=-1))

F29 <- sam[which(sam$sample == "F29"),]
F29$fraction = F29$value / sum(F29$value)
F29$ymax = cumsum(F29$fraction)
F29$ymin = c(0, head(F29$ymax, n=-1))


F32 <- sam[which(sam$sample == "F32"),]
F32$fraction = F32$value / sum(F32$value)
F32$ymax = cumsum(F32$fraction)
F32$ymin = c(0, head(F32$ymax, n=-1))


F34 <- sam[which(sam$sample == "F34"),]
F34$fraction = F34$value / sum(F34$value)
F34$ymax = cumsum(F34$fraction)
F34$ymin = c(0, head(F34$ymax, n=-1))


F35 <- sam[which(sam$sample == "F35"),]
F35$fraction = F35$value / sum(F35$value)
F35$ymax = cumsum(F35$fraction)
F35$ymin = c(0, head(F35$ymax, n=-1))

F41 <- sam[which(sam$sample == "F41"),]
F41$fraction = F41$value / sum(F41$value)
F41$ymax = cumsum(F41$fraction)
F41$ymin = c(0, head(F41$ymax, n=-1))

F43 <- sam[which(sam$sample == "F43"),]
F43$fraction = F43$value / sum(F43$value)
F43$ymax = cumsum(F43$fraction)
F43$ymin = c(0, head(F43$ymax, n=-1))

F47 <- sam[which(sam$sample == "F47"),]
F47$fraction = F47$value / sum(F47$value)
F47$ymax = cumsum(F47$fraction)
F47$ymin = c(0, head(F47$ymax, n=-1))

F53 <- sam[which(sam$sample == "F53"),]
F53$fraction = F53$value / sum(F53$value)
F53$ymax = cumsum(F53$fraction)
F53$ymin = c(0, head(F53$ymax, n=-1))

F54 <- sam[which(sam$sample == "F54"),]
F54$fraction = F54$value / sum(F54$value)
F54$ymax = cumsum(F54$fraction)
F54$ymin = c(0, head(F54$ymax, n=-1))

F87 <- sam[which(sam$sample == "F87"),]
F87$fraction = F87$value / sum(F87$value)
F87$ymax = cumsum(F87$fraction)
F87$ymin = c(0, head(F87$ymax, n=-1))

F88 <- sam[which(sam$sample == "F88"),]
F88$fraction = F88$value / sum(F88$value)
F88$ymax = cumsum(F88$fraction)
F88$ymin = c(0, head(F88$ymax, n=-1))

data <- as.data.frame(rbind(F15, F19, F24, F28, F29, F32, F34, F35, F41, F43, F47, F53, F54, F87, F88))
mean(data$value)


kra <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect() +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  scale_fill_manual(values = c("Klebsiella pneumoniae" = "#e63946", "others" = "#f4a261"), name = "Kraken2") +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        text = element_text(size = 12)) +
  theme_void() +
  facet_wrap(~sample, ncol = 3, nrow = 5) 
kra 

png('figura2.png', res = 600, height = 15, width = 15, units = 'cm')
kra
dev.off()

