library(ape)
library(treeio)
library(ggtree)
library(ggplot2)
library(phylotate)
library(tidytree)
library(ggnewscale)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(viridis)


setwd("/mnt/raid2tb/bubble/tesis/capitulo1/fig7/st11_patric3/P_2023_05_28_213813753074/")
tree <- read.newick(file = "parsnp.tree")
df <- read.csv("/home/csalazar/Descargas/genomes_kp/patric3/BVBRC_genome.csv", header = T)


tre <- read.newick(file = "parsnp.tree")
as_tibble(tre)
head(meta)
colnames(meta)[1] <- c("id")

t <- as.data.frame(as_tibble(tre))
head(t)
t$id <- t$label
t$id <- gsub("f_", "", t$id)
t$id <- gsub(".fna.fasta", "", t$id)

meta <- df
head(meta)
colnames(meta)[1] <- c("id")

metaa <- merge(t, meta, by = "id")
head(metaa)
dim(metaa)
dim(df)
meta[order(meta$Collection.Year),]


anti <- anti_join(t, metaa, by = "id")
anti <- anti[1:15,]
anti$Collection.Year <- c("2017")
anti$Isolation.Country <- c("Uruguay")
anti$Geographic.Group <- c("South America")
anti$Host.Name <- c("Human, Homo sapiens")
head(anti)
anti[15, 7] = "China"
anti[15, 8] = "Asia"

meta <- select(metaa, parent, node, branch.length, label, id, Collection.Year, Isolation.Country, Geographic.Group, Host.Name)
head(meta)


meta <- as.data.frame(rbind(meta, anti))
head(meta)
meta$Isolation.Country[meta$Isolation.Country==""]<-"SD"
meta$Geographic.Group[meta$Geographic.Group==""]<-"SD"



library(RColorBrewer)
coul <- brewer.pal(8, "Paired") 
coul <- colorRampPalette(coul)(15)

cls <- brewer.pal(8, "Set3") 
cls <- colorRampPalette(cls)(15)



p <- ggtree(tre, size =0.5, layout = "rectangular") %<+% meta +
  geom_tippoint(aes(color=Isolation.Country), size = 2) + 
  scale_colour_manual(values = coul) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25)) +
  geom_rootpoint(position = "identity") +
  labs(colour="País") 
p

p + geom_text(aes(label=node), hjust=-.3)


meta2 <- meta
rownames(meta2) <- meta[,4]
clade2 <- rep('Other', dim(meta2)[1])

library(data.table)
target2 <- meta2$Isolation.Country %like% 'Uruguay'
clade2[target2] <- 'target2'
meta2$clade2 <- clade2

grp <- list()
uph <- unique(meta2$clade2)
oth <- c()
nam <- c()
r   <- 1

for (u in 1:length(uph)) {
  
  w <- which(meta2$clade2 == uph[u])
  
  if (length(w) < 1) {
    
    oth <- c(oth, w)
    
  } else {
    
    s <- meta2[w,1]
    grp[[r]] <- s
    nam <- c(nam, uph[u])
    r <- r + 1
  }
}

names(grp) <- uph

target.clade <- MRCA(tre,  grp$target2)  

f <- p + geom_cladelab(node=322, label="*", offset = 0.005, barsize=1.5, fontsize=10, align=T)

f
head(meta2)
e <- collapse(f, node =373) %<+% meta2 +
  geom_tippoint(aes(color=Isolation.Country), size = 2) + 
  geom_tiplab(aes(label=Isolation.Country)) +
  geom_rootpoint(position = "identity") +
  geom_point2(aes(subset=(node==308)), shape=23, size=6, fill="#6c81d9", color = "#6c81d9") +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none") 
e

ggsave("supp_fig7.png", width = 40, height = 40, units = "cm", limitsize = FALSE)



count <- select(meta, label, Geographic.Group)
row.names(count) <- count$label
count$label <- NULL

d <- p + geom_cladelab(node=322, label="*", offset = 0.05, barsize=1.5, fontsize=10)


p1 <- gheatmap(d, count, offset=0.01, width=.1, colnames = F, color = FALSE) + 
    scale_fill_manual(values = cls, name = "Continente")
  
p1

ggsave("fig7.png", width = 40, height = 40, units = "cm", limitsize = FALSE)




library(ggpubr)

fig <- ggarrange(p1, e, ncol = 2, heights = c(0.7, 0.3), widths = c(0.7, 0.3), labels = c("A", "B"), font.label = 20)


e <- collapse(f, node =373) %<+% meta2 +
  geom_tippoint(aes(color=Isolation.Country), size = 2) + 
  geom_tiplab(aes(label=id)) +
  geom_rootpoint(position = "identity") +
  geom_point2(aes(subset=(node==308)), shape=23, size=6, fill="#6c81d9", color = "#6c81d9") +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none") 
e

head(df)
df[grepl("573.7245", df$Genome.ID),]
df[grepl("573.7359", df$Genome.ID),]
df[grepl("573.45119", df$Genome.ID),]
df[grepl("573.45852", df$Genome.ID),] "2017"

