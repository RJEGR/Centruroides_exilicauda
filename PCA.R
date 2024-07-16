

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

dir <- "~/Documents/Centruroides/Results/03.Quantification/align-to-rnaspades-assembly/"
dir_out <- "~/Documents/Centruroides/Plots/"

f <- list.files(path = dir, pattern = "feature_counts.txt$", full.names = T) 

# cols <- read_tsv(f, skip = 1) %>% select(contains(".sorted.bam")) %>% names()

PCA <- function(f) {
  
  df <- read_tsv(f, skip = 1) 
  
  m <- df %>% select(contains(".sorted.bam")) %>% as("matrix")
  
  data <- DESeq2::vst(round(m))
  
  PCA <- prcomp(t(data), center = T, scale. = FALSE)
  
  PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Method = basename(f))
  
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  # percentVar <- round(PCA$sdev/sum(PCA$sdev)*100,1)
  
  PCAvar <- data.frame(
    Eigenvalues = PCA$sdev,
    percentVar = percentVar,
    Varcum = cumsum(percentVar),
    Method = basename(f))
  
  return(list(PCAdf, PCAvar))
}

PCAdf <- PCA(f)

percentVar <- PCAdf[[2]]$percentVar

PCAdf[[1]] %>%
  mutate(Method = "Rnaspades") %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  mutate(LIBRARY_ID = gsub("../|.sorted.bam", "", LIBRARY_ID)) %>%
  mutate(LIBRARY_ID = gsub("_CKDL240013335-1A_22JMN3LT3_L3", "", LIBRARY_ID)) %>%
  mutate(St = substr(LIBRARY_ID, 1,1)) %>%
  ggplot(., aes(PC1, PC2, color = St)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  ggforce::geom_mark_ellipse(aes(group = as.factor(St)),
    fill = 'grey70', color = NA) +
  scale_color_grey("") +
  facet_grid(~ Method) +
  ylim(-100,100) + xlim(-100,100) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_point(size = 7, alpha = 0.7) +
  geom_text( family = "GillSans", mapping = aes(label = LIBRARY_ID), size = 2.5, color = "black") +
  theme(legend.position = "top") -> p

ggsave(p, filename = '03.Quantification.PCA.png', path = dir_out, width = 4, height = 3, device = png, dpi = 300)


PCAdf[[2]] %>%
  group_by(Method) %>%
  mutate(Method = "Rnaspades") %>%
  mutate(Dim = row_number()) %>%
  ggplot(., aes(y = percentVar, x = as.factor(Dim), fill = Method, color = Method)) +
  geom_col(position = position_dodge2(), fill = "gray78", color = "gray78") +
  geom_line(aes(y = Varcum, group = Method)) +
  geom_point(aes(y = Varcum)) +
  # labs(x = "Component Number", y = "Eigenvalue") +
  labs(x = "Principal component", y = "Fraction variance explained (%)") +
  scale_fill_grey("") +
  scale_color_grey("") +
  guides(color=guide_legend(nrow = 1)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "none", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) -> p 

ggsave(p, filename = 'PCA.Eigenevalues.png', path = dir_out, width = 2.5, height = 2.5, device = png, dpi = 300)
