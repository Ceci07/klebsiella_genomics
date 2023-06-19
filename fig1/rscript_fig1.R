library(ggplot2)
library(ggsignif)
library(dplyr)


setwd("./")

dt <- read.table("tabla_conteo_reads.tsv", sep = "\t", header = T)
head(dt)
colnames(dt) <- c("Muestra", "reads_ONT", "largo_ONT", "calidad_ONT", "N50_ONT", "reads_Illumina")
dt$reads_ONT <- gsub("[.]", "", dt$reads_ONT)
dt$reads_ONT <- as.numeric(dt$reads_ONT)
dt$largo_ONT <- gsub("[.]", "", dt$largo_ONT)
dt$largo_ONT <- gsub(",00", "", dt$largo_ONT)
dt$largo_ONT <- gsub(",30", "", dt$largo_ONT)
dt$largo_ONT <- as.numeric(dt$largo_ONT)
dt$calidad_ONT <- as.numeric(dt$calidad_ONT)
dt$N50_ONT <- gsub("[.]", "", dt$N50_ONT)
dt$N50_ONT <- as.numeric(dt$N50_ONT)
dt$reads_Illumina <- gsub("[.]", "", dt$reads_Illumina)
dt$reads_Illumina <- gsub("[,]", "", dt$reads_Illumina)
dt$reads_Illumina <- as.numeric(dt$reads_Illumina)

library(dplyr)
df1 <- select(dt, Muestra, reads_ONT)
colnames(df1) <- c("Muestra", "value")
df1$parameter <- c("Número de lecturas")
mean(df1$value)
sd(df1$value)

df2 <- select(dt, Muestra, largo_ONT)
colnames(df2) <- c("Muestra", "value")
df2$parameter <- c("Largo de lecturas")

df3 <- select(dt, Muestra, calidad_ONT)
colnames(df3) <- c("Muestra", "value")
df3$parameter <- c("Calidad de las lecturas")
mean(df3$value)
sd(df3$value)

df4 <- select(dt, Muestra, N50_ONT)
colnames(df4) <- c("Muestra", "value")
df4$parameter <- c("N50 de las lecturas")


df5 <- select(dt, Muestra, reads_Illumina)
colnames(df5) <- c("Muestra", "value")
df5$parameter <- c("Número de lecturas")
mean(df5$value)
sd(df5$value)

dtt <- as.data.frame(rbind(df1, df2, df3, df4))

mean(df1$value)

p1 <- ggplot(df1, aes(x = Muestra, y = value)) +
  geom_bar(stat = "identity", fill = "#735d78") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank()) +
  ylab("Cantidad") +
  geom_hline(yintercept = 106325.3, linetype = "dashed") +
  ggtitle("Lecturas de secuenciación ONT")

p1

p2 <- ggplot(df2, aes(x = Muestra, y = value)) +
  geom_bar(stat = "identity", fill = "#b392ac") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank()) +
  geom_hline(yintercept = 200, linetype = "dotted") +
  ylab("Largo promedio") 
p2

p3 <- ggplot(df3, aes(x = Muestra, y = value)) +
  geom_bar(stat = "identity", fill = "#d1b3c4") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank()) +
  ylab("Calidad promedio") +
  scale_y_continuous(breaks=seq(0, 13, 2)) +
  geom_hline(yintercept = 8, linetype = "dotted")
p3

p5 <- ggplot(df5, aes(x = Muestra, y = value)) +
  geom_bar(stat = "identity", fill = "#e8c2ca") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "bold")) +
  ylab("Canitidad") +
  geom_hline(yintercept = 1481288, linetype = "dashed") +
  ggtitle("Lecturas de secuenciación Illumina (Fwd/Rev)")
p5

fig1 <- ggarrange(p1, p2, p3, p5, nrow = 4, common.legend = T, legend = "right", heights = c(0.33, 0.33, 0.33), widths =c(0.33, 0.33, 0.33), align = "v", 
                  labels = c("A", "B", "C", "D"))

png('fig1_cap1.png', res = 600, height = 20, width = 20, units = 'cm')
fig1
dev.off()
