# ============================================
# RainCloudPlot for tumor and marginal exp
# ============================================

# ---- Main function ----
generate_raincloud <- function(
    file_path,               # [REQUIRED] Path to Excel file
    gene_name,               # [REQUIRED] Name of the gene
    col_case,                # [REQUIRED] Column name for Case/Tumor
    col_ctrl,                # [REQUIRED] Column name for Control/Margin
    test_method,             # "mann_whitney" or "wilcoxon"
    # --- Optional Arguments ---
    label_case = "Tumor",    
    label_control = "Marginal", 
    color_case = "#FF6347",  
    color_control = "#4682B4"  
    
) {
  
  # A. Validation Checks
  if (!file.exists(file_path)) stop("Error: Excel file not found.")
  
  raw_data <- read_excel(file_path)
  
  if (!all(c(col_case, col_ctrl) %in% names(raw_data))) {
    stop("Error: Column names not found in the Excel file.")
  }
  
  # B. Setup Output Directory
  out_dir <- file.path("Results_RainCloud", gene_name)
  if (!dir.exists(out_dir)) dir_create(out_dir)
  
  # C. Perform Statistical Test
  case_vals <- raw_data[[col_case]]
  ctrl_vals <- raw_data[[col_ctrl]]
  test_method <- tolower(test_method)
  
  if (test_method == "mann_whitney") {
    test_res <- wilcox.test(case_vals, ctrl_vals, paired = FALSE)
    test_name_str <- "Mann-Whitney U"
  } else if (test_method == "wilcoxon") {
    test_res <- wilcox.test(case_vals, ctrl_vals, paired = TRUE)
    test_name_str <- "Wilcoxon Signed-Rank"
  } else {
    stop("Error: Invalid test_method.")
  }
  
  p_val <- test_res$p.value
  
  # Determine Significance Stars
  sig_stars <- case_when(
    p_val < 0.0001 ~ "****",
    p_val < 0.001  ~ "***",
    p_val < 0.01   ~ "**",
    p_val < 0.05   ~ "*",
    TRUE           ~ "ns"
  )
  
  # Stats Label String (Gray, inside frame)
  stats_label_text <- paste0(
    test_name_str, " | P-value: ", format.pval(p_val, digits = 3), " (", sig_stars, ")"
  )
  
  # D. Prepare Data
  plot_df <- data.frame(
    Group = factor(rep(c(label_case, label_control), each = nrow(raw_data)), 
                   levels = c(label_control, label_case)), 
    Expression = c(case_vals, ctrl_vals)
  )
  
  color_map <- setNames(c(color_case, color_control), c(label_case, label_control))
  
  # E. Create RainCloud Plot
  p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group, color = Group)) +
    
    # 1. Cloud
    stat_halfeye(
      adjust = 0.5, width = 0.6, justification = -0.2, 
      .width = 0, point_colour = NA, alpha = 0.6
    ) +
    
    # 2. Boxplot
    geom_boxplot(
      width = 0.15, outlier.shape = NA, alpha = 0.5, 
      color = "black", lwd = 0.6
    ) +
    
    # 3. Rain
    geom_jitter(
      width = 0.1, alpha = 0.6, size = 1.8
    ) +
    
    # 4. Styling
    scale_fill_manual(values = color_map) +
    scale_color_manual(values = color_map) +
    labs(
      title = NULL,      
      subtitle = NULL,
      # Y-axis: Uses bquote to render Greek Beta symbol
      y = bquote("Relative Expression of " ~ .(gene_name) ~ "/" ~ beta ~ "-actin"),
      x = NULL
    ) +
    
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12, face = "bold", color = "black"),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
      legend.position = "none",
      panel.grid = element_blank(), 
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      axis.line = element_blank() 
    ) +
    
    # 5. Internal Annotation (Smaller font size = 3)
    annotate("text", x = 1.5, y = Inf, label = stats_label_text, 
             vjust = 2, size = 3, color = "gray40", fontface = "bold")
  
  # =========================================================
  # F. Save Results (CHANGED HERE FOR 300 DPI TIFF)
  # =========================================================
  
  # Changed extension to .tiff and added explicit units and compression
  ggsave(
    filename = file.path(out_dir, paste0(gene_name, "_RainCloud.tiff")), 
    plot = p, 
    width = 6, 
    height = 6, 
    units = "in",       # Fixes dimension issues
    dpi = 300,          # High resolution
    compression = "lzw" # Standard compression
  )
  
  # Save Stats Excel
  stats_out <- data.frame(
    Gene = gene_name,
    Test = test_name_str,
    P_Value = p_val,
    Significance = sig_stars,
    Mean_Case = mean(case_vals),
    Mean_Control = mean(ctrl_vals)
  )
  write_xlsx(stats_out, file.path(out_dir, paste0(gene_name, "_Stats.xlsx")))
  
  message(sprintf(">> Success! Plot saved in: %s", normalizePath(out_dir)))
  return(p)
}

# ============================================
# EXAMPLE USAGE
# ============================================

generate_raincloud(
  file_path = "filename.xlsx",
  gene_name = "P53",
  col_case  = "TumoralRelativeExp",
  col_ctrl  = "MarginalRelativeExp",
  label_case = "Tumor",
  label_control = "Marginal",
  test_method = "wilcoxon"
)

