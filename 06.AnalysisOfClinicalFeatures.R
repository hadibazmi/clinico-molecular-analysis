# ============================================
# ANALYSIS OF CLINICAL FEATURES
# ============================================

# ---- Main function ----
# 1. Setup
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, writexl, ggplot2, dplyr)

# 2. Worker Function
analyze_single_gene_tissue = function(
    dataset, gene_col, clinical_cols, gene_name, tissue_label,              
    cutoff_list = list(), output_dir_base = "Results", color_palette = NULL
) {
  
  gene_dir = file.path(output_dir_base, gene_name)
  if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)
  
  if (is.null(color_palette)) {
    color_palette = c("#512DA8", "#E91E63", "#1976D2", "#FF6F00", "#FFB300")
  }
  
  results_accumulator = list()
  full_y_label = paste0("Relative Expression of ", gene_name, "\n(", tissue_label, ")")
  
  # RESTORING ORIGINAL THEME
  custom_theme = theme_minimal() +
    theme(
      axis.text = element_text(size = 12, face = "bold", color = "black"),
      axis.title.y = element_text(size = 12, face = "bold", color = "black", margin = margin(r = 10)),
      axis.title.x = element_text(size = 12, face = "bold", color = "black", margin = margin(t = 10)),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      axis.line = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  
  for (clin_var in clinical_cols) {
    res = tryCatch({
      clean_var = trimws(clin_var)
      if (!clean_var %in% names(dataset)) next
      
      df_sub = dataset %>%
        select(Expr = all_of(gene_col), Clin = all_of(clean_var)) %>%
        filter(!is.na(Expr), !is.na(Clin))
      
      if (nrow(df_sub) < 2) next
      
      is_num = is.numeric(df_sub$Clin)
      is_cutoff = clean_var %in% names(cutoff_list)
      
      if (is_cutoff || !is_num) {
        if (is_cutoff) {
          val = cutoff_list[[clean_var]]
          df_sub$Group = factor(ifelse(df_sub$Clin < val, paste0("< ", val), paste0(">= ", val)))
          tname = "Mann-Whitney"
        } else {
          df_sub$Group = as.factor(df_sub$Clin)
          tname = ifelse(length(unique(df_sub$Group)) == 2, "Mann-Whitney", "Kruskal-Wallis")
        }
        
        n_grps = length(unique(df_sub$Group))
        if (n_grps < 2) next
        
        # Statistics
        pval = if(n_grps == 2) wilcox.test(Expr ~ Group, data = df_sub)$p.value else kruskal.test(Expr ~ Group, data = df_sub)$p.value
        sig = ifelse(pval < 0.05, "*", "ns")
        
        # PLOT - RESTORING ORIGINAL STYLE
        internal_label = paste0(clean_var, "\n(", tname, " | P: ", format.pval(pval, digits=3), " ", sig, ")")
        center_x = (n_grps + 1) / 2
        
        p = ggplot(df_sub, aes(x = Group, y = Expr, fill = Group)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black", lwd = 0.6) +
          geom_jitter(width = 0.1, alpha = 0.6, size = 1.8) +
          scale_fill_manual(values = color_palette) +
          labs(y = full_y_label, x = NULL) +
          custom_theme +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
          annotate("text", x = center_x, y = Inf, label = internal_label, vjust = 1.5, size = 3.5, color = "gray40", fontface = "bold")
        
        ggsave(file.path(gene_dir, paste0(tissue_label, "_", clean_var, ".tiff")), p, width = 6, height = 6, dpi = 300, compression = "lzw")
        
        # EXCEL DATA - SEPARATED COLUMNS
        df_sub %>%
          group_by(SubGroup = Group) %>%
          summarise(N = n(), Mean = round(mean(Expr), 4), SD = round(sd(Expr), 4), .groups = "drop") %>%
          mutate(Gene = gene_name, Tissue = tissue_label, Variable = clean_var, Test = tname, P_Value = pval, Sig = sig)
        
      } else {
        # Correlation Case
        ctest = cor.test(df_sub$Clin, df_sub$Expr, method = "spearman")
        
        # PLOT CORRELATION STYLE
        center_x = mean(range(df_sub$Clin, na.rm = TRUE))
        internal_label = paste0(clean_var, "\n(Spearman r: ", round(ctest$estimate, 3), " | P: ", format.pval(ctest$p.value, digits=3), ")")
        
        p = ggplot(df_sub, aes(x = Clin, y = Expr)) +
          geom_point(color = color_palette[1], alpha = 0.6, size = 2) +
          geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
          labs(y = full_y_label, x = clean_var) +
          custom_theme +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
          annotate("text", x = center_x, y = Inf, label = internal_label, vjust = 1.5, size = 3.5, color = "gray40", fontface = "bold")
        
        ggsave(file.path(gene_dir, paste0(tissue_label, "_", clean_var, "_Correlation.tiff")), p, width = 6, height = 6, dpi = 300, compression = "lzw")
        
        data.frame(Gene = gene_name, Tissue = tissue_label, Variable = clean_var, SubGroup = "Continuous",
                   N = nrow(df_sub), Mean = NA, SD = NA, Test = "Spearman", P_Value = ctest$p.value, Sig = ifelse(ctest$p.value<0.05,"Sig","ns"))
      }
    }, error = function(e) return(NULL))
    results_accumulator[[clin_var]] = res
  }
  return(bind_rows(results_accumulator))
}

# 3. Manager Function
run_multi_gene_analysis = function(dataset, gene_config_list, clinical_cols, cutoff_list, 
                                   label_case = "Tumor", label_control = "Margin", color_palette = NULL) {
  master_list = list()
  for (gene_info in gene_config_list) {
    message("Analyzing: ", gene_info$name)
    t_res = analyze_single_gene_tissue(dataset, gene_info$case, clinical_cols, gene_info$name, label_case, cutoff_list, color_palette = color_palette)
    m_res = analyze_single_gene_tissue(dataset, gene_info$control, clinical_cols, gene_info$name, label_control, cutoff_list, color_palette = color_palette)
    master_list[[gene_info$name]] = bind_rows(t_res, m_res)
  }
  final_df = bind_rows(master_list)
  write_xlsx(final_df, "Master_Clinical_Results.xlsx")
  message(">> DONE! Plots restored to original style & Excel updated with separate Mean/SD columns.")
}





# ============================================
# EXAMPLE USAGE
# ============================================

my_data = read_excel("filename.xlsx")
class(my_data$`FAB-Subtype`)
my_data$`FAB-Subtype`=as.factor(my_data$`FAB-Subtype`)
# NOTE: Since we reverted to the original style, remember to convert 
# your categorical columns (like FAB-Subtype) to factors MANUALLY here if needed.
# Example: my_data$`FAB-Subtype` = as.factor(my_data$`FAB-Subtype`)

genes_config = list(
  list(
    name = "P53",             
    case = "TumoralRelativeExp",    
    control = "MarginalRelativeExp" 
  )
)

clinical_vars = c(
  "Age", 
  "Histology", 
  "Differentiation", 
  "Tumor-size", 
  "Lymph-node", 
  "Depth-of-cervical-invasion", 
  "Stage", 
  "HPV"
)

my_cutoffs = list("Age" = 50)


run_multi_gene_analysis(
  dataset = my_data,
  gene_config_list = genes_config,
  clinical_cols = clinical_vars,
  cutoff_list = my_cutoffs,
  color_palette = c("#512DA8", "#E91E63", "#1976D2", "#FF6F00", "#FFB300"),
  label_case = "Tumor",
  label_control = "Margin"
)
