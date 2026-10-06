# Load packages -----------------------------------------------------------
library(tidyverse)
library(tidyplots)
library(patchwork)
library(nlme)
library(emmeans)


# Load data ---------------------------------------------------------------

# Read in named data (long)
long_named_data <- readr::read_csv("Data/SPP_S1 - long_named_data_asca.csv") %>% 
  mutate(variable = make.names(variable))


# Clean & process data ----------------------------------------------------
metabolites <- unique(long_named_data$variable)
timepoint_levels <- c("Pre", "Post", "8 h post", "22 h post")
session_levels <- c("Control", "Session 1", "Session 2")

# Pivot from long to wide
wide_named_data <- long_named_data %>% 
  mutate(session = factor(session, levels = session_levels), 
         sample_time = factor(sample_time, levels = timepoint_levels)) %>% 
  select(subject_id, session, sample_time, variable, value) %>% 
  pivot_wider(names_from = variable, values_from = value)


# Set up where results will be stored -------------------------------------
lmm_df <- data.frame(measure = character(), 
                     pval_time = numeric(), 
                     pval_session = numeric(), 
                     pval_int = numeric(), 
                     stringsAsFactors = FALSE)

emm_list <- list()


# Iterate over each metabolite --------------------------------------------
for (i in metabolites){
  print(i)
  
  ml <- as.formula(paste(i, "~ sample_time * session"))
  
  m1.nlme <- tryCatch({
    lme(ml, random = ~1|subject_id, 
        weights = varIdent(form = ~1|sample_time), 
        data = wide_named_data, 
        na.action = na.exclude)
  }, error = function(e){
    message("Error fitting model for variable ", i, ": ", e)
    return(NULL)
  })
  
  if(is.null(m1.nlme)) next
  
  # Omnibus ANOVA p-values 
  m1.nlme_anova <- anova(m1.nlme)
  pval_time    <- m1.nlme_anova$`p-value`[2]
  pval_session <- m1.nlme_anova$`p-value`[3]
  pval_int     <- m1.nlme_anova$`p-value`[4]
  
  lmm_df <- rbind(lmm_df, data.frame(measure = i,
                                     pval_time = pval_time,
                                     pval_session = pval_session,
                                     pval_int = pval_int,
                                     stringsAsFactors = FALSE))
  
  # Estimated marginal means at each timepoint x session
  emm <- emmeans(m1.nlme, specs = ~ sample_time * session)
  emm_list[[i]] <- as.data.frame(emm) %>% mutate(measure = i)
}


# Combine results across metabolites --------------------------------------

# Omnibus p-values with FDR correction
lmm_df$pval_time_fdr <- p.adjust(lmm_df$pval_time, method = "BH")
lmm_df$pval_session_fdr <- p.adjust(lmm_df$pval_session, method = "BH")
lmm_df$pval_int_fdr <- p.adjust(lmm_df$pval_int, method = "BH")

# Round to fixed number of decimal places
lmm_df <- lmm_df %>% 
  mutate(across(starts_with("pval"), ~round(.x, 4)))

# Combine all stored emmeans into one long data frame 
emm_df <- bind_rows(emm_list)

# Keep only metabolites with a significant interaction after FDR correction
sig_metabolites <- lmm_df %>% 
  filter(pval_int_fdr < 0.05)


# Metabolite families -----------------------------------------------------

acylcarnities <- c("carnitine", "l_acetylcarnitine_car_2_0", "propionylcarnitine_car_3_0",
                   "butyryl_isobutyrylcarnitine_car_4_0", "X3_hydroxybutyrylcarnitine_car_4_0_oh",
                   "tiglylcarnitine_car_5_1", "isovaleryl_valeryl_2_methylbutyryl_carnitine_car_5_0",
                   "glutarylcarnitine_car_5_0_dc", "hexanoylcarnitine_car_6_0", "X2_octenoylcarnitine_car_8_1",
                   "octanoylcarnitine_car_8_0", "decanoylcarnitine_car_10_0", "decenoylcarnitine_car_10_1",
                   "hydroxydecanoylcarnitine_car_10_0_oh", "lauroylcarnitine_car_12_0", "dodecenoylcarnitine_car_12_1",
                   "tetradecanoylcarnitine_car_14_0", "tetradecenoylcarnitine_car_14_1", "tetradecadienoylcarnitine_car_14_2", 
                   "hexadecenoylcarnitine_car_16_1", "palmitoylcarnitine_car_16_0", "octadecenoylcarnitine_car_18_1",
                   "octadecadienoylcarnitine_car_18_2", "stearoylcarnitine_car_18_0")

fatty_acids <- c("X12_hydroxydodecanoic_acid", "X5_hydroxydecanoate", "X2_hydroxyhexadecanoic_acid",
                 "myristic_acid", "heptadecanoic_acid", "linoleic_acid", "elaidic_acid",
                 "ricinoleic_acid", "nonadecanoic_acid", "phytanic_acid", "pinolenic_acid",
                 "X4z_7z_10z_13z_16z_19z_4_7_10_13_1_6_19_docosahexaenoic_acid", "behenic_acid", 
                 "tricosanoic_acid", "nervonic_acid", "tetracosanoic_acid", "pentacosanoic_acid", "hexacosanoic_acid")

amino_acid <- c("arginine", "citrulline", "homoarginine", "histidine", "l_histidine",
                "l_alanine", "l_glutamine", "l_methionine", "phenylalanine", "l_phenylalanine",
                "tryptophan", "l_tryptophan", "l_kynurenine", "l_homoserine", "l_norleucine",
                "d_serine", "d_cysteine", "l_cystine", "o_tyrosine", "n_acetyl_l_arginine",
                "symmetric_dimethylarginine", "glutamic_acid_1_3_nph_carbonyl_carboxy_phospho",
                "n_methyl_a_aminoisobutyric_acid")

tca_intermediates <- c("cis_aconitic_acid_1_3_nph_carbonyl_carboxy_phospho", "cis_aconitic_acid_3_3_nph_carbonyl_carboxy_phospho",
                       "succinic_acid_1_3_nph_carbonyl_carboxy_phospho", "succinic_acid_2_3_nph_carbonyl_carboxy_phospho",
                       "succinic_acid_2x_3_nph", "malic_acid_2_3_nph_carbonyl_carboxy_phospho", "pyruvic_acid_2x_3_nph",
                       "lactic_acid_1x_3_nph", "hydroxypropionic_acid", "citraconic_acid_1_3_nph_carbonyl_carboxy_phospho",
                       "citraconic_acid_2_3_nph_carbonyl_carboxy_phospho", "succinylacetone")


# Display names -----------------------------------------------------------

pretty_names <- c(
  # Acylcarnitines
  carnitine = "carnitine",
  l_acetylcarnitine_car_2_0 = "L-acetylcarnitine (C2:0)", 
  propionylcarnitine_car_3_0 = "propionylcarnitine (C3:0)", 
  butyryl_isobutyrylcarnitine_car_4_0 = "butyryl/isobutyrylcarnitine (C4:0)", 
  X3_hydroxybutyrylcarnitine_car_4_0_oh = "3-hydroxybutyrylcarnitine (C4:0-OH)", 
  tiglylcarnitine_car_5_1 = "tiglylcarnitine (C5:1)", 
  isovaleryl_valeryl_2_methylbutyryl_carnitine_car_5_0 = "isovaleryl/valeryl/2-methylbutyrylcarnitine (C5:0)", 
  glutarylcarnitine_car_5_0_dc = "glutarylcarnitine (C5:0-DC)", 
  hexanoylcarnitine_car_6_0 = "hexanoylcarnitine (C6:0)", 
  X2_octenoylcarnitine_car_8_1 = "trans-2-octenoylcarnitine (C8:1)", 
  octanoylcarnitine_car_8_0 = "octanoylcarnitine (C8:0)", 
  decanoylcarnitine_car_10_0 = "decanoylcarnitine (C10:0)", 
  decenoylcarnitine_car_10_1 = "decenoylcarnitine (C10:1)", 
  hydroxydecanoylcarnitine_car_10_0_oh = "hydroxydecanoylcarnitine (C10:0-OH)", 
  lauroylcarnitine_car_12_0 = "lauroylcarnitine (C12:0)", 
  dodecenoylcarnitine_car_12_1 = "dodecenoylcarnitine (C12:1)", 
  tetradecanoylcarnitine_car_14_0 = "tetradecanoylcarnitine (C14:0)", 
  tetradecenoylcarnitine_car_14_1 = "tetradecenoylcarnitine (C14:1)", 
  tetradecadienoylcarnitine_car_14_2 = "tetradecadienoylcarnitine (C14:2)", 
  palmitoylcarnitine_car_16_0 = "palmitoylcarnitine (C16:0)", 
  hexadecenoylcarnitine_car_16_1 = "hexadecenoylcarnitine (C16:1)", 
  stearoylcarnitine_car_18_0 = "stearoylcarnitine (C18:0)", 
  octadecenoylcarnitine_car_18_1 = "octadecenoylcarnitine (C18:1)", 
  octadecadienoylcarnitine_car_18_2 = "octadecadienoylcarnitine (C18:2)",
  # Amino acids
  arginine = "arginine", 
  citrulline = "citrulline", 
  homoarginine = "homoarginine", 
  histidine = "histidine", 
  l_histidine = "L-histidine", 
  l_alanine = "L-alanine", 
  l_glutamine = "L-glutamine", 
  l_methionine = "L-methionine", 
  phenylalanine = "phenylalanine", 
  l_phenylalanine = "L-phenylalanine",
  tryptophan = "tryptophan", 
  l_tryptophan = "L-tryptophan", 
  l_kynurenine = "L-kynurenine", 
  l_homoserine = "L-homoserine", 
  l_norleucine = "L-norleucine",
  d_serine = "D-serine", 
  d_cysteine = "D-cysteine", 
  l_cystine = "L-cystine", 
  o_tyrosine = "O-tyrosine", 
  n_acetyl_l_arginine = "n-acetyl-l-arginine (NA-Arg)",
  symmetric_dimethylarginine = "symmetric dimethylarginine (SDMA)", 
  glutamic_acid_1_3_nph_carbonyl_carboxy_phospho = "glutamic acid (3-NPH-cpp 1)",
  n_methyl_a_aminoisobutyric_acid = "n-methyl-alpha-aminoisobutyric acid (MeAIB)",
  # TCA intermediates
  cis_aconitic_acid_1_3_nph_carbonyl_carboxy_phospho = "cis-aconitic acid (3-NPH-ccp 1)", 
  cis_aconitic_acid_3_3_nph_carbonyl_carboxy_phospho = "cis-aconitic acid (3-NPH-ccp 3)",
  succinic_acid_1_3_nph_carbonyl_carboxy_phospho = "succinic acid (3-NPH-ccp 1)", 
  succinic_acid_2_3_nph_carbonyl_carboxy_phospho = "succinic acid (3-NPH-ccp 2)",
  succinic_acid_2x_3_nph = "succinic acid (di-3-NPH)", 
  malic_acid_2_3_nph_carbonyl_carboxy_phospho = "malic acid (3-NPH-ccp 2)", 
  pyruvic_acid_2x_3_nph = "pyruvic acid (di-3-NPH)",
  lactic_acid_1x_3_nph = "lactic acid (3-NPH)", 
  hydroxypropionic_acid = "hydroxypropionic acid", 
  citraconic_acid_1_3_nph_carbonyl_carboxy_phospho = "citraconic acid (3-NPH-ccp 1)",
  citraconic_acid_2_3_nph_carbonyl_carboxy_phospho = "citraconic acid (3-NPH-ccp 2)", 
  succinylacetone = "succinylacetone",
  # Fatty acids
  X12_hydroxydodecanoic_acid = "12-hydroxydodecanoic acid (C12:0-OH)", 
  X5_hydroxydecanoate = "5-hydroxydecanoate (C10:0-OH)", 
  X2_hydroxyhexadecanoic_acid = "2-hydroxyhexadecanoic acid (C16:0-OH)",
  myristic_acid = "myristic acid (C14:0)", 
  heptadecanoic_acid = "heptadecanoic acid (C17:0)", 
  linoleic_acid = "linoleic acid (C18:2)", 
  elaidic_acid = "elaidic acid (trans-C18:1)",
  ricinoleic_acid = "ricinoleic acid (C18:1-OH)", 
  nonadecanoic_acid = "nonadecanoic acid (C19:0)", 
  phytanic_acid = "phytanic acid (C20:0)", 
  pinolenic_acid = "pinolenic acid (C18:3)",
  X4z_7z_10z_13z_16z_19z_4_7_10_13_1_6_19_docosahexaenoic_acid = "docosahexaenoic acid (DHA, C22:6)", 
  behenic_acid = "behenic acid (C22:0)", 
  tricosanoic_acid = "tricosanoic acid (C23:0)", 
  nervonic_acid = "nervonic acid (C24:1)", 
  tetracosanoic_acid = "tetracosanoic acid (C24:0)", 
  pentacosanoic_acid = "pentacosanoic acid (C25:0)", 
  hexacosanoic_acid = "hexacosanoic acid (C26:0)"
)


# Heatmaps ----------------------------------------------------------------

# Family lookup
family_lookup <- bind_rows(
  data.frame(measure = acylcarnities, family = "Acylcarnitines"), 
  data.frame(measure = fatty_acids, family = "Fatty acids"), 
  data.frame(measure = amino_acid, family = "Amino acids"), 
  data.frame(measure = tca_intermediates, family = "TCA intermediates")
)

# Join family + FDR p-value onto all emmeans results (not just significant ones)
emm_df_all <- emm_df %>% 
  left_join(family_lookup, by = "measure") %>% 
  filter(!is.na(family)) %>% 
  left_join(lmm_df %>% select(measure, pval_int_fdr), by = "measure") %>% 
  group_by(measure) %>% 
  mutate(emmean_z = as.numeric(scale(emmean))) %>% 
  ungroup() %>% 
  mutate(measure_pretty = coalesce(unname(pretty_names[measure]), measure),
         measure_label  = measure_pretty)

# Check which metabolites still need a display name
setdiff(unique(emm_df_all$measure), names(pretty_names))

# Order metabolites within a family by clustering similarity of their pattern
order_within_family <- function(df) {
  wide_mat <- df %>%
    distinct(measure, sample_time, session, emmean_z) %>%
    unite(cell, sample_time, session) %>%
    pivot_wider(names_from = cell, values_from = emmean_z) %>%
    column_to_rownames("measure") %>%
    as.matrix()
  hc <- hclust(dist(wide_mat))
  rownames(wide_mat)[hc$order]
}

# Shared colour limits
z_lim <- max(abs(emm_df_all$emmean_z), na.rm = TRUE)

# Heatmap for a single family, using ALL its metabolites
heatmap_plot <- function(fam_name) {
  df <- emm_df_all %>% filter(family == fam_name)
  ord <- order_within_family(df)
  df <- df %>% mutate(measure_label = factor(measure_label, 
                                             levels = measure_label[match(ord, measure)]))
  
  # Append a star to labels of significant metabolites
  sig_lookup <- df %>% 
    distinct(measure_label, pval_int_fdr) %>% 
    mutate(star_label = if_else(pval_int_fdr < 0.05, 
                                paste0(measure_label, " \u2731"),   # ✱ heavy asterisk
                                as.character(measure_label)))
  
  df %>% 
    ggplot(aes(x = sample_time, y = measure_label, fill = emmean_z)) + 
    geom_tile(color = "white", linewidth = 0.3) +
    facet_grid(~session) +
    scale_fill_gradient2(low = "#2e5597", mid = "white", high = "#fb495c", 
                         midpoint = 0, 
                         limits = c(-z_lim, z_lim), 
                         name = "Z-scored\nEMM") + 
    scale_y_discrete(labels = function(x) sig_lookup$star_label[match(x, sig_lookup$measure_label)]) + 
    labs(title = fam_name, x = NULL, y = NULL) +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 24, face = "italic"), 
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = 14), 
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), 
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14), 
          strip.background = element_rect(fill = "#cdc9c9"), 
          strip.text = element_text(size = 16, face = "bold"))
}

# Generate one family's heatmap at a time
acylcarnitine_heatmap <- heatmap_plot("Acylcarnitines")
ggsave("Images/acylcarnitine_heatmap.png", acylcarnitine_heatmap, 
       width = 10, height = 6.5, units = "in", dpi = 600)

fatty_acid_heatmap <- heatmap_plot("Fatty acids")
ggsave("Images/ffa_heatmap.png", fatty_acid_heatmap, 
       width = 10, height = 5, units = "in", dpi = 600)

amino_acid_heatmap <- heatmap_plot("Amino acids")
ggsave("Images/aa_heatmap.png", amino_acid_heatmap, 
       width = 10, height = 6, units = "in", dpi = 600)

tca_intermediates_heatmap <- heatmap_plot("TCA intermediates")
ggsave("Images/tca_heatmap.png", tca_intermediates_heatmap, 
       width = 10, height = 4, units = "in", dpi = 600)

# Put plots into 2 figures

# Figure 1: Acylcarnitines + Amino acids (side by side)
fig1 <- (acylcarnitine_heatmap | amino_acid_heatmap) +
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 24, face = 'bold'))

# Figure 2: Fatty acids + TCA intermediates (side by side)
fig2 <- (fatty_acid_heatmap | tca_intermediates_heatmap) +
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 24, face = "bold"))

# 9 x 5.6 in leaves room for a 3-line caption on a landscape letter page
# with 1 in margins (usable area is 9 x 6.5 in)
ggsave("Images/fig_acyl_aa.png", fig1, width = 18, height = 8, dpi = 600)
ggsave("Images/fig_fa_tca.png",  fig2, width = 18, height = 8, dpi = 600)



