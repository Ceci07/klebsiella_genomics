library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)

setwd("./")

df <- read.csv("nano_card_out.tsv", sep = "\t", header = T)
head(df)
dim(df)
df <- df %>% filter(X.COVERAGE > 95)
df <- df %>% filter(X.IDENTITY > 90)

sp <- as.data.frame(str_split_fixed(df$SEQUENCE, "_", 2))
head(sp)

df <- as.data.frame(cbind(df, sp$V1))
colnames(df)[16] <- c("sample")

list <- df %>%
  group_by(sample, GENE, RESISTANCE) %>%
  tally() %>%
  ungroup()
list <- as.data.frame(list)
head(list)
list$tec <- c("ONT")

dff <- read.csv("hy_card_out.tsv", sep = "\t", header = T)
head(dff)
dim(dff)
dff <- dff %>% filter(X.COVERAGE > 95)
dff <- dff %>% filter(X.IDENTITY > 90)

spp <- as.data.frame(str_split_fixed(dff$SEQUENCE, "_", 2))
head(spp)

dff <- as.data.frame(cbind(dff, spp$V1))
colnames(dff)[16] <- c("sample")

list2 <- dff %>%
  group_by(sample, GENE, RESISTANCE) %>%
  tally() %>%
  ungroup()
list2 <- as.data.frame(list2)
head(list2)
list2$tec <- c("ONT+Illumina")

head(s)
head(ss)
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
        axis.text.y = element_text(size = 12),
  ) +
  scale_fill_viridis(option = "F", direction = -1) +
  facet_wrap(~tec) +
  labs(fill= "Conteo" ) + xlab("Muestras") + ylab("") +
  ggtitle("CARD")
fig

dir.create("figures")
png("./figures/ram.png", res = 600, height = 25, width = 30, units = "cm")
fig
dev.off()


################################################################################
df <- read.csv("nano_plasmid.tsv", sep = "\t", header = T)
head(df)
dim(df)

sp <- as.data.frame(str_split_fixed(df$SEQUENCE, "_", 2))
head(sp)

df <- as.data.frame(cbind(df, sp$V1))
colnames(df)[16] <- c("sample")

df[which(df$sample == "F15"),]
df <- df %>% filter(X.COVERAGE > 98)
df <- df %>% filter(X.IDENTITY > 90)


list <- df %>%
  group_by(sample, GENE, RESISTANCE) %>%
  tally() %>%
  ungroup()
list <- as.data.frame(list)
head(list)
list$tec <- c("ONT")

dff <- read.csv("hy_plasmid.tsv", sep = "\t", header = T)
head(dff)
dim(dff)
dff <- dff %>% filter(X.COVERAGE > 98)
dff <- dff %>% filter(X.IDENTITY > 90)

spp <- as.data.frame(str_split_fixed(dff$SEQUENCE, "_", 2))
head(spp)

dff <- as.data.frame(cbind(dff, spp$V1))
colnames(dff)[16] <- c("sample")

list2 <- dff %>%
  group_by(sample, GENE, RESISTANCE) %>%
  tally() %>%
  ungroup()
list2 <- as.data.frame(list2)
head(list2)
list2$tec <- c("ONT+Illumina")

head(s)
head(ss)
d <- as.data.frame(rbind(list, list2))
head(d)

figp <- ggplot(d, aes(sample, GENE, fill = n)) + 
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
  ggtitle("plasmidfinder")
figp

png("./figures/plasmid.png", res = 600, height = 25, width = 30, units = "cm")
figp
dev.off()



library(ggpubr)


a <- ggarrange(fig, figp, nrow = 2, heights = c(0.7,0.3), widths = c(0.7,0.3), align = "v", labels = c("A", "B"), font.label = list(size = 20), common.legend = T, legend = "right")
a
png("./figures/figura4.png", res = 600, height = 30, width = 25, units = "cm")
a
dev.off()




