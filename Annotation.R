

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

dir <- "~/Documents/Centruroides/Results/04.Annotation/"

dir_out <- "~/Documents/Centruroides/Plots/"


f <- list.files(path = dir, pattern = "Trinotate.xls$", full.names = T) 

annot <- read_tsv(f, na = ".")

names(annot)[1] <- "gene_id"

# Using only ORFS

# annot %>% drop_na(SignalP) %>% view()

# annot <- annot %>% filter(!is.na(prot_id))

# OR, using only blastp

annot <- annot %>% drop_na(gene_ontology_BLASTP)

# annot %>% filter(grepl("Trpv1", sprot_Top_BLASTP_hit))

# annot %>% filter(grepl("Transient receptor potential cation", ))


# Plot NOGS

NOG.col <- c('J@Translation, ribosomal structure and biogenesis','A@RNA processing and modification','K@Transcription','L@Replication, recombination and repair','B@Chromatin structure and dynamics','D@Cell cycle control, cell division, chromosome partitioning','Y@Nuclear structure','V@Defense mechanisms','T@Signal transduction mechanisms','M@Cell wall/membrane/envelope biogenesis','N@Cell motility','Z@Cytoskeleton','W@Extracellular structures','U@Intracellular trafficking, secretion, and vesicular transport','O@Posttranslational modification, protein turnover, chaperones','X@Mobilome: prophages, transposons','C@Energy production and conversion','G@Carbohydrate transport and metabolism','E@Amino acid transport and metabolism','F@Nucleotide transport and metabolism','H@Coenzyme transport and metabolism','I@Lipid transport and metabolism','P@Inorganic ion transport and metabolism','Q@Secondary metabolites biosynthesis, transport and catabolism','R@General function prediction only','S@Function unknown')

into <- c("EggNM.COG_category", "COG_name")

NOG.col <- data.frame(NOG.col) %>% separate(NOG.col, sep = "@", into = into)



# 1) Split anchnoDB_BLASTP and peptides
# annot %>% drop_na(peptide)

anchnoDB_NOGS <- annot %>% drop_na(anchnoDB_BLASTP) %>%
  left_join(NOG.col, by = c("EggNM.COG_category")) %>%
  drop_na(COG_name) %>%
  count(COG_name, sort = T) %>%
  mutate(frac = n/sum(n), facet = "B) anchnoDB") 

signalP_NOGS <- annot %>% drop_na(SignalP) %>%
  left_join(NOG.col, by = c("EggNM.COG_category")) %>%
  drop_na(COG_name) %>%
  count(COG_name, sort = T) %>%
  mutate(frac = n/sum(n), facet = "c) SignalP") 



eggNOG_cols <- c("gene_id")

NOGS_DB <- annot  %>%  
  select_at(vars(contains(eggNOG_cols), starts_with("EggNM"))) %>%
  left_join(NOG.col, by = c("EggNM.COG_category"))


plotdf <- NOGS_DB %>%
  drop_na(COG_name) %>%
  count(COG_name, sort = T) %>%
  mutate(frac = n/sum(n)) %>%
  arrange(desc(frac)) %>%
  mutate(COG_name = factor(COG_name, levels = unique(COG_name))) %>%
  mutate(facet = "A) Transcriptome") %>% 
  rbind(anchnoDB_NOGS, signalP_NOGS)
  # filter(COG_name != "Function unknown")

data_text <- plotdf %>% 
  # group_by(COG_name, facet) %>% summarise(n = sum(n)) %>%
  mutate(n = paste0("(", n, ")"))

plotdf %>% 
  ggplot(aes(y = COG_name, x = frac, group = facet)) +
  labs(y = "Nested Orthologous Gene Group (NOGS)", x = "Enrichment ratio") +
  facet_grid(~ facet, scales = "free_x", space = "free_x", switch = "y") +
  geom_col(fill = "black") +
  theme_bw(base_size = 12, base_family = "GillSans") +
  xlim(0,0.4) +
  theme(legend.position = "top", 
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = 'grey95', color = 'white')) +
  geom_text(data = data_text, 
    aes(label = n), size = 2.5,
    hjust = -0.1, vjust = 0, 
    family = "GillSans", position = position_dodge(width = 1)) -> p

ggsave(p, filename = 'NOGS.png', path = dir_out, width = 8, height = 5, device = png, dpi = 300)

# Transcriptome blast results by sp -----


# df <- split_blast(annot %>% drop_na(anchnoDB_BLASTP))

# df %>% count(name, sort = T) %>% view()

hist(blastp$identity)

lineagedf <- blastp %>% 
  count(lineage,sort = T) %>%
  separate(lineage, sep = ";", into = paste0("Rank", rep(1:8))) %>%
  filter(Rank1 == "Eukaryota") %>%
  mutate(frac = n/sum(n)) %>% 
  mutate(x = Rank6, main_group = gsub(" ", "", Rank2))

lineagedf <- lineagedf %>% filter(frac > 1E-4)

order_r <- lineagedf %>% group_by(x) %>% tally(frac, sort = T) %>% pull(x)

lineagedf <- lineagedf %>% mutate(x = factor(x, levels = rev(order_r)))

labels <-  unique(lineagedf$main_group)

n_colors <- length(labels)

if(n_colors > 7) {
  library(ggsci)
  pal <- colorRampPalette(pal_igv()(7))(n_colors)
} else {
  pal <- pal_igv("default")(n_colors)
}

names(pal) <- labels


library(ggnested)

ggnested(lineagedf, 
  aes(x = x, 
    y = frac, 
    main_group = main_group, 
    sub_group = Rank4), main_palette = pal,
  legend_title = "Lineages groups") +
  geom_col() +
  labs(y = "Frac. of transcripts", x = "Lineage of annotated genes (Uniprot Blast)") +
  coord_flip() +
  # theme_bw(base_size = 12, base_family = "GillSans") +
  theme(
    text = element_text(size=12, family="GillSans"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid.minor.x = element_blank(),
    panel.grid = element_line(colour = "grey92"),
    strip.background = element_rect(fill = 'grey95', color = 'white')) -> p

# ggsave(p, filename = 'LINEAGE.png', path = dir_out, width = 6, height = 6, device = png, dpi = 300)

p <- p +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank()) 

p1 <- blastp %>% 
  select(lineage, identity, evalue) %>%
  separate(lineage, sep = ";", into = paste0("Rank", rep(1:8))) %>%
  filter(Rank1 == "Eukaryota") %>%
  mutate(Rank6 = factor(Rank6, levels = rev(order_r))) %>%
  drop_na(Rank6) %>% 
  ggplot(aes(x = identity, y = Rank6)) +
  labs(y = "Lineage of annotated genes (Uniprot Blast)", x = "Identity (%)")

p1 <- p1 + ggridges::geom_density_ridges() + xlim(0, 100)

library(patchwork)

p1 <- p1 +
  theme(
    text = element_text(size=12, family="GillSans"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid.minor.x = element_blank(),
    panel.grid = element_line(colour = "grey92"),
    strip.background = element_rect(fill = 'grey95', color = 'white'))


ps <- p1 + plot_spacer() + p + plot_layout(widths = c(0.3, -0.1, 0.5))

ps

ggsave(ps, filename = 'LINEAGE.png', path = dir_out, width = 8, height = 6, device = png, dpi = 300)

# focus on Protein toxins ------
# 0) Plot Pfam ontology of protein toxin by using 1, 2 and 3 subjects:

# 1) Related Centruroides blast
  
blastp <- annot %>% split_blast("BLASTP_hit")

q.genus <- c("Tityus", "Centruroides")

genusdf <- blastp %>% distinct(gene, transcript, genus) 

q.genusdf <- genusdf %>% filter(genus %in% q.genus)

q <- blastp %>% filter(genus %in% q.genus) %>% pull(transcript)

q.annot <- annot %>% filter(transcript_id %in% q)

Pfam.q <- q.annot %>% split_gene_ontology(hit = "Pfam") %>% 
  left_join(q.genusdf) %>% distinct() 

# 2) Related anchnoDB_BLASTP ----

Pfam.anchnoDB <- annot %>% 
  drop_na(anchnoDB_BLASTP) %>%
  split_gene_ontology(hit = "Pfam") %>%
  filter(!transcript %in% q) %>%
  mutate(genus = "anchnoDB") %>%
  distinct() 


# 3) Related signalP results ----

Pfam.SignalP <- annot %>% drop_na(SignalP) %>% 
  split_gene_ontology(hit = "Pfam") %>%
  filter(!transcript %in% q) %>%
  mutate(genus = "SignalP") 

grep_names <- "channel|inhibitor|protein-coupled|receptor"

# filter_names <- rbind(Pfam.anchnoDB, Pfam.q) %>% distinct(name) %>% pull()

Pfam.SignalP <- Pfam.SignalP %>%
  filter(grepl(grep_names, name)) 

Pfam.SignalP %>% count(name)

# 4) Grep other channel inhibitor activity in Pfam -----

Pfamdf <- annot %>% split_gene_ontology(hit = "Pfam") %>% 
  filter(!transcript %in% q) %>%
  mutate(genus = "Unknown")

Pfamdf <- Pfamdf %>% filter(grepl(grep_names, name))

Pfamdf %>% count(genus, ontology, name,sort = T) %>% view()
# Dataviz ----

DATAVIZ <- rbind(Pfam.q, Pfam.anchnoDB, Pfam.SignalP,Pfamdf) %>% 
  count(genus, ontology, name,sort = T) 

DATAVIZ %>%   group_by(name) %>% tally(n, sort = T) %>% pull(name) -> sort_names

unique(DATAVIZ$genus)
sort_genus <- c("Centruroides","Tityus","SignalP","anchnoDB","Unknown")


p <- DATAVIZ %>%
  mutate(genus = factor(genus, levels = sort_genus)) %>% 
  mutate(name = factor(name, levels = rev(sort_names))) %>% 
  ggplot() +
  geom_col(aes(y = name, x = n, fill = genus), position = position_stack(reverse = T)) +
  # scale_fill_grey("") +
  see::scale_fill_bluebrown() +
  labs(y = "Associated toxin Gene Ontology", x = "Number of transcrips (Pfam Domains)", fill = "") +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = 'grey95', color = 'white'))
  # see::scale_fill_metro(name = "", reverse = T) 
  # geom_tile(aes(fill = n, y = name, x = genus)) + 
  # facet_grid(~ genus, scales = "free", space = "free")


ggsave(p, filename = 'TOXINS.png', path = dir_out, width = 8.2, height = 10, device = png, dpi = 300)

# plot known anchnoDB toxins -----

DF <- annot %>% 
  drop_na(anchnoDB_BLASTP) %>%
  # drop_na(SignalP) %>%
  split_blast("BLASTP", prefix = "anchnoDB_") 



# DF <- DF %>% left_join(genusdf, by = "transcript")

sort_toxin <- DF %>% count(uniprot, sort = T) %>% pull(uniprot)

p1 <- DF %>% count(uniprot, sort = T) %>% 
  mutate(uniprot = factor(uniprot, levels = rev(sort_toxin))) %>%
  ggplot() + 
  geom_col(aes(y = uniprot, x = n), fill = "black") +
  labs(x = "Number of transcripts") +
  theme(
    text = element_text(size=10, family="GillSans"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid.minor.x = element_blank(),
    panel.grid = element_line(colour = "grey92"),
    strip.background = element_rect(fill = 'grey95', color = 'white'))

p1 <- p1 + theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank(),
  axis.line.y = element_blank()) 

p2 <- DF %>% 
  mutate(uniprot = factor(uniprot, levels = rev(sort_toxin))) %>%
  ggplot(aes(y = uniprot, x = identity,fill = stat(x))) +
  labs(y = "Protein toxin name (anchnoDB)", x = "Identity (%)") +
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 0.5, point_alpha = 1, alpha = 0.7) + 
  xlim(0, 100) +
  scale_fill_viridis_c(option = "C") +
  theme(legend.position = "none",
    text = element_text(size=10, family="GillSans"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid.minor.x = element_blank(),
    panel.grid = element_line(colour = "grey92"),
    strip.background = element_rect(fill = 'grey95', color = 'white'))

ps <- p2 + plot_spacer() + p1 + plot_layout(widths = c(0.3, -0.1, 0.5))

ggsave(ps, filename = 'anchnoDB.png', path = dir_out, width = 6, height = 6, device = png, dpi = 300)







