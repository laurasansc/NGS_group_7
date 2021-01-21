library(tidyverse)
library(ggthemes)

# read file
df <- read.csv("/STATS/basic_stats.csv")

# Run T-test

require(graphics)
raw_blood <- with(df, raw_reads[type=="blood"])
raw_saliva <- with(df, raw_reads[type=="saliva"])
t.test(x=raw_blood, y=raw_saliva)      # P = .0532

aln_blood <- with(df, reads_aligned[type=="blood"])
aln_saliva <- with(df, reads_aligned[type=="saliva"])
t.test(x=aln_blood, y=aln_saliva)      # P = .05828

in_blood <- with(df, median_insert_size[type=="blood"])
in_saliva <- with(df, median_insert_size[type=="saliva"])
t.test(x=in_blood, y=in_saliva)      # P = .0651


q_blood <- with(df, avg_quality[type=="blood"])
q_saliva <- with(df, avg_quality[type=="saliva"])
t.test(x=q_blood, y=q_saliva)      # P = .428


# STATS PLOTS
png("raw_reads.png", width = 800, height = 600, res = 150)
df %>% ggplot() +
  geom_boxplot(aes(y = raw_reads, x = type, fill = type, color = type), fill = c("#cb7f7f", "#7fc7ce"), alpha = 0.6) +
  geom_point(aes(y = raw_reads, x = type, fill = type)) +
  scale_y_continuous(name = "Number of  raw reads") +
  scale_x_discrete(name = "Sample type", labels = c("Blood", "Saliva")) +
  scale_fill_manual(
    labels = c("Blood", "Saliva"), values = c("#cb7f7f", "#7fc7ce"),
    breaks = c("blood", "saliva"),
    name = "Sample type"
  ) +
  scale_color_manual(
    labels = c("Blood", "Saliva"), values = c("#cb7f7f", "#7fc7ce"),
    breaks = c("blood", "saliva"),
    name = "Sample type"
  ) +
  theme_base()
dev.off()

png("aligned_reads.png", width = 800, height = 600, res = 150)
df %>% ggplot() +
  geom_boxplot(aes(y = reads_aligned, x = type, fill = type, color = type), fill = c("#cb7f7f", "#7fc7ce"), alpha = 0.6) +
  geom_point(aes(y = reads_aligned, x = type, fill = type)) +
  scale_y_continuous(name = "Number of aligned reads") +
  scale_x_discrete(name = "Sample type", labels = c("Blood", "Saliva")) +
  scale_fill_manual(
    labels = c("Blood", "Saliva"), values = c("#cb7f7f", "#7fc7ce"),
    breaks = c("blood", "saliva"),
    name = "Sample type"
  ) +
  scale_color_manual(
    labels = c("Blood", "Saliva"), values = c("#cb7f7f", "#7fc7ce"),
    breaks = c("blood", "saliva"),
    name = "Sample type"
  ) +
  theme_base()
dev.off()
png("paired_reads.png", width = 800, height = 600, res = 150)
df %>% ggplot() +
  geom_boxplot(aes(y = reads_aligned_paired, x = type, fill = type, color = type), fill = c("#cb7f7f", "#7fc7ce"), alpha = 0.6) +
  geom_point(aes(y = reads_aligned_paired, x = type, fill = type)) +
  scale_y_continuous(name = "Number of aligned and paired reads") +
  scale_x_discrete(name = "Sample type", labels = c("Blood", "Saliva")) +
  scale_fill_manual(
    labels = c("Blood", "Saliva"), values = c("#cb7f7f", "#7fc7ce"),
    breaks = c("blood", "saliva"),
    name = "Sample type"
  ) +
  scale_color_manual(
    labels = c("Blood", "Saliva"), values = c("#cb7f7f", "#7fc7ce"),
    breaks = c("blood", "saliva"),
    name = "Sample type"
  ) +
  theme_base()
dev.off()

png("insert.png", width = 800, height = 600, res = 150)
df %>% ggplot() +
  geom_boxplot(aes(y = median_insert_size, x = type, fill = type, color = type), fill = c("#cb7f7f", "#7fc7ce"), alpha = 0.6) +
  geom_point(aes(y = median_insert_size, x = type, fill = type)) +
  scale_y_continuous(name = "Median insert size") +
  scale_x_discrete(name = "Sample type", labels = c("Blood", "Saliva")) +
  scale_fill_manual(
    labels = c("Blood", "Saliva"), values = c("#cb7f7f", "#7fc7ce"),
    breaks = c("blood", "saliva"),
    name = "Sample type"
  ) +
  scale_color_manual(
    labels = c("Blood", "Saliva"), values = c("#cb7f7f", "#7fc7ce"),
    breaks = c("blood", "saliva"),
    name = "Sample type"
  ) +
  theme_base()
dev.off()

png("avg_quality.png", width = 800, height = 600, res = 150)
df %>% ggplot() +
  geom_boxplot(aes(y = avg_quality, x = type, fill = type, color = type), fill = c("#cb7f7f", "#7fc7ce"), alpha = 0.6) +
  geom_point(aes(y = avg_quality, x = type, fill = type)) +
  scale_y_continuous(name = "Average quality") +
  scale_x_discrete(name = "Sample Type", labels = c("Blood", "Saliva")) +
  scale_fill_manual(
    labels = c("Blood", "Saliva"), values = c("#cb7f7f", "#7fc7ce"),
    breaks = c("blood", "saliva"),
    name = "Sample type"
  ) +
  scale_color_manual(
    labels = c("Blood", "Saliva"), values = c("#cb7f7f", "#7fc7ce"),
    breaks = c("blood", "saliva"),
    name = "Sample type"
  ) +
  theme_base()
dev.off()

# Coverage plot
cov <- read.csv("STATS/coverage.txt", sep = "\t", header = FALSE)

cov <- cov %>%
  rename(
    group_steps = V1,
    steps = V2,
    SRR8595490 = V3,
    SRR8595491 = V4,
    SRR8595492 = V5,
    SRR8595493 = V6,
    SRR8595494 = V7,
    SRR8595495 = V8,
    SRR8595496 = V9,
    SRR8595497 = V10
  )

covNew <- cov %>% pivot_longer(cols = SRR8595490:SRR8595497, 
                                    names_to = "ID", 
                                    values_to = "Value")

covNew <- covNew %>% mutate(log_val = log(Value))
head(covNew)

COLORS <- c(SRR8595490 = "red", SRR8595491 ="orange",  
            SRR8595492 = "deeppink" ,SRR8595493 = "darkred",
            SRR8595494 = "cyan", SRR8595495="cornflowerblue", 
            SRR8595496="blue", SRR8595497="darkturquoise"
            )

library(scales)

png("coverage.png", width = 1800, height = 1200, res = 300)
ggplot(covNew, aes(x = steps, y = log(Value), group = ID, color = ID)) +
  geom_line(size = 0.5) +
  xlim(0, 70) +
  #ylim(12,20)+
  scale_color_manual(values = COLORS,name = "Sample ID")+
  theme(legend.position="bottom", 
        panel.background = element_rect(fill="white", colour="black", size=0.5), 
        panel.grid.major = element_line(colour="grey", size=0.3),
        panel.grid.minor = element_line(colour="grey", size=0.3)) +
  xlab("Coverage") + ylab("Number of mapped reads") +
 scale_y_continuous(breaks = 1:5,
               labels = c(100000,"1e+06","1e+07","1e+08","1e+09"),
               limits = c(12,20))  +
 annotation_logticks(sides="l") 
dev.off()

#### UPSET PLOTS #######

library(UpSetR)
input <- c(
  Blood = 121854,
  Saliva = 118659,
  "Blood&Saliva"=3710769
)
upset(fromExpression(input), 
      nintersects = 3, 
      nsets = 3, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
)

blood = c(rep(1,121854+3710769), rep(0,118659))
saliva = c(rep(0,121854), rep(1,118659+3710769))
dat <- data.frame(blood, saliva)
x<-c(rep("KnownSNPs", 3703817), rep("NovelSNPs", 6952), rep("KnownSNPs",112596), rep("NovelSNPs", 6063),rep("KnownSNPs", 115589), rep("NovelSNPs", 6265))
m<-cbind(dat, x)  
upset(m,
      queries = list(
        list(query = elements, 
             params = list("x", "KnownSNPs"), color = "#e69f00", active = T),
        list(query = elements, 
             params = list("x", "NovelSNPs"), color = "#cc79a7", active = T)))


## indels
blood_i = c(rep(1,98163+711990), rep(0,99676))
saliva_i = c(rep(0,98163), rep(1,99676+711990))
dat_i <- data.frame(blood_i, saliva_i)
x_i<-c(rep("KnownIndels",87437), rep("NovelIndels",10726), rep("KnownIndels", 701550), rep("NovelIndels", 10440), rep("KnownIndels", 87235), rep("NovelIndels", 12441))
m_i<-cbind(dat_i, x_i)  
upset(m_i,
      queries = list(
        list(query = elements, 
             params = list("x_i", c("KnownIndels","NovelIndels")), color = "grey", active = T),
        list(query = elements, 
             params = list("x_i", "NovelIndels"), color = "#e69f00", active = T)))


library(ComplexUpset)
groups <- c("blood_i","saliva_i")
upset(
  dat_i,
  groups,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      aes(fill=x_i)
    )
  ),
  width_ratio=0.5
)