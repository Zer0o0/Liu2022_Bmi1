library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggsci)
library(scico)
library(hrbrthemes)
library(scales)
library(ggsci)
library(pheatmap)
library(rtracklayer)

## Pie chart showing the distribution of peaks across the genome
gff <- readGFF("mm10.refGene.gtf")
id_map <- gff[, c("transcript_id", "gene_id", "gene_name")]
id_map <- unique(id_map)
f <- "path_to_peak_annotation_file"
data <- read_tsv(f)
data2 <- data %>%
  select(1:4, 8, 10, 12:17, 19) %>%
  separate(col = Annotation, into = c("feature", "annotation"), sep = " \\(") %>%
  separate(col = annotation, into = c("transcript_id", "xxx"), sep = "[,\\)]")
data2 <- merge(data2, id_map, by = "transcript_id", all.x = T)
data2 <- select(data2, 2:6, 1, 17, 8:15)
feature2 <- as.character(data2$feature)
feature2[feature2 == "non-coding"] <- "Intergenic"
feature2[feature2 == "intron"] <- "Intron"
feature2[feature2 == "TTS"] <- "Exon"
feature2[feature2 == "3' UTR"] <- "Exon"
feature2[feature2 == "5' UTR"] <- "Exon"
feature2[feature2 == "exon"] <- "Exon"
feature2[feature2 == "promoter-TSS"] <- "Promoter"
feature2[feature2 == "promoter"] <- "Promoter"
data2$feature2 <- feature2
data2 <- data2[!is.na(feature2), ]
feature2 <- as.character(data2$feature2)
feature3 <- c()
dist <- as.numeric(as.character(data2$`Distance to TSS`))
for (i in 1:length(dist)) {
  if (dist[i] > (-2500) & dist[i] < 500 & (feature2[i] != "Promoter")) {
    feature3 <- c(feature3, "Promoter")
  } else {
    feature3 <- c(feature3, feature2[i])
  }
}
data2$feature3 <- feature3
write.table(data2, paste0(f, "2"), row.names = F, sep = "\t")
mypal <- pal_npg("nrc", alpha = 0.6)(4)
show_col(mypal)
names(mypal) <- c("Intergenic", "Intron", "Exon", "Promoter")

data_pie <- data2 %>%
  group_by(feature3) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(prop = round(n / sum(n), 3) * 100) %>%
  mutate(ypos = cumsum(prop) - 0.5 * prop)
data_pie$feature3 <- factor(data_pie$feature3, levels = rev(data_pie$feature3))
data_pie$prop2 <- sprintf("%0.2f", data_pie$prop)
data_pie <- na.omit(data_pie)
mypalx <- mypal[as.character(data_pie$feature3)]
data_pie %>%
  ggplot(aes(x = 2, y = prop, fill = feature3)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  # scale_fill_npg(alpha=0.8,)+
  scale_fill_manual(values = mypalx) +
  # geom_text(aes(label = paste0(prop, "%")), position = position_stack(vjust=0.5),colour='white')+
  geom_text(aes(y = ypos, label = paste0(prop2, "%")), color = "black", size = 4, nudge_x = 0.8) +
  theme_void() +
  labs(fill = "Features") +
  theme(text = element_text(size = 14), legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"))

## Display gene enrichment analysis results from DAVID using bubble plots
# https://david.ncifcrf.gov/home.jsp
file_go <- "path_to_david_result_table_file"
go <- read_tsv(file_go)

go2 <- go %>%
  separate(Term, sep = "~", into = c("idx", "term")) %>%
  mutate(logp = log10(PValue) * (-1)) %>%
  select(term, logp, `Fold Enrichment`, `%`) %>%
  arrange(desc(logp)) %>%
  slice(1:10)
colnames(go2) <- c("Term", "log10pvalue", "Fold Enrichment", "Percentage")
scico(11, palette = "devon", direction = -1)
go2 %>%
  ggplot(aes(x = `Fold Enrichment`, y = reorder(Term, `Fold Enrichment`), size = Percentage, colour = log10FDR)) +
  geom_point() +
  scale_colour_gradient(low = "#BDB7F1", high = "#2E629F") +
  scale_size_continuous(range = c(2, 6)) +
  scale_x_continuous(breaks = c(1.5, 2.0, 2.5), labels = c("1.5", "2.0", "2.5")) +
  labs(y = "", colour = quote(-log[10] ~ pvalue), size = "Percentage(%)") +
  theme_light() +
  theme(
    panel.border = element_rect(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(size = 14),
    legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")
  )

#gene expression matrix and visualization
mat <- "path_to_expression_value_file"
anno_col <- data.frame("type" = rep(c("WT", "KO"), 3), row.names = colnames(mat))
cols <- pal_npg("nrc", alpha = 0.8)(2)
png("heatmap.png", width = 14, height = 8, units = "cm", res = 1200)
pheatmap(mat,
  scale = "row",
  color = colorRampPalette(scico(15, palette = "vikO", direction = 1)[2:14])(100),
  show_colnames = F, fontsize_row = 14,
  treeheight_row = 0, treeheight_col = 15,
  cellwidth = 34, cellheight = 34 * 0.618,
  # cutree_cols = 2
  annotation_col = anno_col, annotation_colors = list("Celltype" = c(WT = cols[1], KO = cols[2])),
  legend_breaks = seq(-1.6, 1.6, 0.8), legend_labels = seq(-1.6, 1.6, 0.8),
)
dev.off()

#bar chart with error line
file <- "path_to_file" #
data <- read_tsv(file)
data %>%
  pivot_longer(cols = rep1:rep3, names_to = "rep") %>%
  group_by(group, time) %>%
  summarise(mean = mean(value), sd = sd(value), .groups = "drop") %>%
  mutate(time2 = factor(time, levels = c("G0/G1", "S", "G2/M")), group2 = factor(group, levels = c("NC", "KD"))) %>%
  ggplot(aes(x = time2, y = mean, fill = group2)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.8, show.legend = F) +
  geom_signif(map_signif_level = TRUE, y_position = c(54, 84), xmin = c(0.8, 1.8), xmax = c(1.2, 2.2), tip_length = 0.02, annotations = "***", textsize = 3, size = 0.4, vjust = 0.2) +
  scale_fill_manual(values = mycol, name = "") +
  geom_errorbar(aes(ymin = mean, ymax = mean + sd), width = 0.3, position = position_dodge(0.8)) +
  ylim(0, 100) +
  theme_pubr() +
  labs(x = "", y = "Percent of Cell Number (%)", title = "") +
  theme(axis.title = element_text(size = 10))
