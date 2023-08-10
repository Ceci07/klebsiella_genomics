library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)

setwd("./")

df <- read.csv("nano_vfdb.tsv", sep = "\t", header = T)
head(df)
dim(df)
df <- df %>% filter(X.COVERAGE > 95)
df <- df %>% filter(X.IDENTITY > 95)

grepl("f28",df$SEQUENCE)

sp <- as.data.frame(str_split_fixed(df$SEQUENCE, "_", 2))
head(sp)

df <- as.data.frame(cbind(df, sp$V1))
colnames(df)[16] <- c("sample")

list <- df %>%
  group_by(sample, GENE, PRODUCT) %>%
  tally() %>%
  ungroup()
list <- as.data.frame(list)
head(list)
list$tec <- c("ONT")

dff <- read.csv("hybrid_vfdb.csv", sep = "\t", header = T)
head(dff)
dim(dff)

dff <- dff %>% filter(X.COVERAGE > 95)
dff <- dff %>% filter(X.IDENTITY > 90)

spp <- as.data.frame(str_split_fixed(dff$SEQUENCE, "_", 2))
head(spp)

dff <- as.data.frame(cbind(dff, spp$V1))
colnames(dff)[16] <- c("sample")

list2 <- dff %>%
  group_by(sample, GENE, PRODUCT) %>%
  tally() %>%
  ungroup()
list2 <- as.data.frame(list2)
head(list2)
list2$tec <- c("ONT+Illumina")
list2$sample <-gsub("f28", "F28", list2$sample)

d <- as.data.frame(rbind(list, list2))
head(d)


fig <- ggplot(d, aes(sample, GENE, fill = n)) + 
  geom_tile(colour = "gray50") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size = 14),
        axis.text.y = element_text(size = 12)) +
  scale_fill_viridis(option = "F", direction = -1) +
  facet_wrap(~tec) +
  labs(fill= "Conteo" ) + xlab("Muestras") + ylab("") +
  ggtitle("VFDB")
fig

dir.create("figures")
png("./figures/vfdb.png", res = 600, height = 25, width = 30, units = "cm")
fig
dev.off()

