# ============================================
# Normality Tests & Plots
# ============================================

# ---- Main function ----
normality <- function(excel_path, column, preferred_name = NULL) {
  
  # 1) Read Excel file
  if (!file.exists(excel_path)) stop(paste("File not found:", excel_path))
  dat <- read_excel(excel_path)
  
  # 2) Resolve the column
  if (is.character(column)) {
    if (!column %in% names(dat)) stop(sprintf("Column '%s' not found.", column))
    x <- dat[[column]]
    col_label <- column
  } else if (is.numeric(column) && length(column) == 1) {
    x <- dat[[column]]
    col_label <- names(dat)[column]
  } else {
    stop("Column must be name or index.")
  }
  
  # 3) Output label
  out_label <- if (is.null(preferred_name) || preferred_name == "") col_label else preferred_name
  
  # 4) Create folder
  if (!dir.exists(out_label)) dir.create(out_label)
  
  # 5) Validate Data
  if (!is.numeric(x)) stop("Selected column is not numeric.")
  x <- x[!is.na(x)]
  if (length(x) < 3) stop("Not enough observations (n < 3).")
  if (sd(x) == 0) stop("Standard deviation is zero.")
  
  # 6) Run Statistical Tests
  shapiro_res <- tryCatch(shapiro.test(x), error = function(e) NULL)
  shapiro_p   <- if (!is.null(shapiro_res)) shapiro_res$p.value else NA
  
  ad_res      <- tryCatch(ad.test(x), error = function(e) NULL)
  ad_stat     <- if (!is.null(ad_res)) unname(ad_res$statistic) else NA
  ad_p        <- if (!is.null(ad_res)) ad_res$p.value else NA
  
  lillie_res  <- tryCatch(lillie.test(x), error = function(e) NULL)
  lillie_p    <- if (!is.null(lillie_res)) lillie_res$p.value else NA
  
  # Descriptive stats
  n_val   <- length(x)
  mean_x  <- mean(x)
  sd_x    <- sd(x)
  skew_x  <- mean((x - mean_x)^3) / sd_x^3
  kurt_x  <- mean((x - mean_x)^4) / sd_x^4
  
  # 7) Prepare Results Dataframe
  results <- data.frame(
    Label          = out_label,
    Source         = col_label,
    N              = n_val,
    Mean           = mean_x,
    SD             = sd_x,
    Skewness       = skew_x,
    Kurtosis       = kurt_x,
    Shapiro_P      = shapiro_p,
    Anderson_P     = ad_p,
    Lilliefors_P   = lillie_p,
    stringsAsFactors = FALSE
  )
  
  # 8) Save Excel
  excel_out <- file.path(out_label, paste0(out_label, "_Results.xlsx"))
  write_xlsx(list(Stats = results), path = excel_out)
  
  # 9) Create Professional Plots (Combined)
  df_plot <- data.frame(val = x)
  
  # A. Histogram with Density
  p_hist <- ggplot(df_plot, aes(x = val)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, 
                   fill = "#5D8AA8", color = "white", alpha = 0.8) + # Teal color
    stat_function(fun = dnorm, args = list(mean = mean_x, sd = sd_x), 
                  color = "#B22222", linewidth = 1.2) + # Firebrick red curve
    labs(title = "Histogram & Normal Curve",
         subtitle = sprintf("Mean: %.2f | SD: %.2f", mean_x, sd_x),
         x = "Value", y = "Density") +
    theme_light() +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, color = "gray30"))
  
  # B. QQ-Plot
  p_qq <- ggplot(df_plot, aes(sample = val)) +
    stat_qq(color = "#2F4F4F", alpha = 0.7, size = 2) + # Dark Slate Gray dots
    stat_qq_line(color = "#B22222", linewidth = 1.2, linetype = "dashed") +
    labs(title = "Q-Q Plot",
         subtitle = sprintf("Shapiro P: %.3f | Lillie P: %.3f", shapiro_p, lillie_p),
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_light() +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, color = "gray30"))
  
  # Combine plots side-by-side using arrangeGrob to create an object
  combined_plot <- arrangeGrob(p_hist, p_qq, ncol = 2, 
                               top = textGrob(paste("Normality Analysis:", out_label), 
                                              gp = gpar(fontsize = 14, font = 2)))
  
  # ---------------------------------------------------------
  # SAVE AS TIFF (300 DPI)
  # ---------------------------------------------------------
  plot_file <- file.path(out_label, paste0("Plot_", out_label, ".tiff"))
  
  ggsave(
    filename = plot_file, 
    plot = combined_plot, 
    width = 10, 
    height = 5, 
    units = "in",       # Ensure units are inches
    dpi = 300,          # High resolution
    compression = "lzw" # Standard lossless compression
  )
  
  message(sprintf(">> Processed: %s | Saved in folder: %s", out_label, out_label))
  return(results)
}


# ============================================
# EXAMPLE USAGE
# ============================================

VariableName = normality("filename.xlsx", "Column1", "Column2")
