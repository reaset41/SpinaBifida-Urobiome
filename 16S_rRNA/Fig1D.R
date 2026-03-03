library(tidyr)
library(tidyverse)
library(patchwork)
library(phyloseq)

# subset by SampleType
ps_urine <- subset_samples(ps_merged, SampleType == "Urine")
ps_stool <- subset_samples(ps_merged, SampleType == "Stool")
ps_hand  <- subset_samples(ps_merged, SampleType == "Hand")
# get dataframes for each
get_family_abundance_df <- function(ps) {
  # Agglomerate at Family
  ps_family <- tax_glom(ps, taxrank = "Family")
  
  # Transform counts to relative abundance
  ps_relabund <- transform_sample_counts(ps_family, function(x) x / sum(x))
  
  # Melt to long format dataframe
  df <- psmelt(ps_relabund)
  
  # Keep only relevant columns and rename
  df %>%
    select(Sample, SampleType, Family, Abundance) %>%
    group_by(Sample, Family) %>%
    summarise(Abundance = sum(Abundance)) %>%
    ungroup()
}

df_urine <- get_family_abundance_df(ps_urine) %>% mutate(SampleType = "Urine")
df_stool <- get_family_abundance_df(ps_stool) %>% mutate(SampleType = "Stool")
df_hand  <- get_family_abundance_df(ps_hand) %>% mutate(SampleType = "Hand")

# combine for plotting
df_all <- bind_rows(df_urine, df_stool, df_hand)

top_families <- df_all %>%
  group_by(Family) %>%
  summarise(total_abund = sum(Abundance)) %>%
  arrange(desc(total_abund)) %>%
  slice_head(n = 9) %>%
  pull(Family)

df_all <- df_all %>%
  mutate(Family = ifelse(Family %in% top_families, Family, "Other")) %>%
  group_by(Sample, SampleType, Family) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup()

df_all_top10 <- df_all %>%
  group_by(Sample) %>%
  mutate(Family = fct_reorder(Family, Abundance, .fun = sum, .desc = TRUE)) %>%
  ungroup()

family_levels <- unique(df_all_top10$Family)

meta_df <- data.frame(sample_data(ps_urine))
meta_df$Sample <- rownames(meta_df)

target_family <- "Enterobacteriaceae"

sample_order <- df_all_top10 %>%
  filter(Family == target_family) %>%
  arrange(desc(Abundance)) %>%
  pull(Sample) %>%
  unique()

df_all_top10$Sample <- factor(df_all_top10$Sample, levels = sample_order)
meta_df$Sample <- factor(meta_df$Sample, levels = sample_order)


meta_long <- meta_df %>%
  select(Sample, Abx_30days, CIC) %>%
  pivot_longer(
    cols = c(Abx_30days, CIC),
    names_to = "Metadata",
    values_to = "Value"
  ) %>%
  mutate(fill_var = paste(Metadata, Value, sep = "_"))

meta_colors <- c(
  # Abx_30days palette
  "Abx_30days_Y" = "#1b9e77",
  "Abx_30days_N"  = "#D3D3D3",
  
  # CIC palette
  "CIC_Y" = "#d95f02",
  "CIC_N"  = "#D3D3D3"
)

meta_plot <- ggplot(meta_long,
                    aes(x = Sample, y = Metadata, fill = fill_var)) +
  geom_tile() +
  scale_fill_manual(values = meta_colors) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )


named_colors <- setNames(
  c(
    "#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#fdbf6f",
    "#ff7f00",
    "#6a3d9a",
    "gray80"),
  c(setdiff(family_levels, "Other"), "Other")
)
plot_family_bar <- function(df, sample_type) {
  ggplot(df %>% filter(SampleType == sample_type), aes(x = Sample, y = Abundance, fill = Family)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(sample_type, "Samples"), y = "Relative Abundance", x = "Sample") +
    #scale_fill_brewer(palette = "")+
    #scale_fill_manual(values = c(RColorBrewer::brewer.pal(9, "Set3"), "Other" = "gray80"))+
    scale_fill_manual(values = named_colors)+
    theme(
      axis.text.x = element_blank(),   # removes x-axis text labels
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=14, color="black"),
      axis.title.y = element_text(size=18, color="black") # removes x-axis tick marks
    )
}

p_urine <- plot_family_bar(df_all_top10, "Urine")


final_plot <- meta_plot / p_urine +
  plot_layout(heights = c(0.15, 1))

final_plot

sggsave(filename="Fig1D.pdf")