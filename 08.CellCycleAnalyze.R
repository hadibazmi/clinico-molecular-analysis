# ==============================================================================
# Cell Cycle Analyze
# ==============================================================================

# ---- Main function ----
analyze_cell_cycle <- function(
    input_data, 
    # phase_list: A named list where each element is a vector of c(control_col, case_col)
    # Example: list("G1" = c("G1_Ctrl", "G1_Treat"), "S" = c("S_Ctrl", "S_Treat"))
    phase_list,
    control_label = "Control",
    case_label    = "Treatment",
    x_axis_title  = "Cell Cycle Phase",
    y_axis_title  = "% Population",
    palette       = NULL, 
    save_dir      = "CellCycle-Results"
) {
  
  # --- INIT ---
  dataset <- NULL
  if (is.character(input_data)) {
    if (file.exists(input_data)) { dataset <- read_excel(input_data) } 
    else { stop("Error: File not found.") }
  } else { dataset <- input_data }
  
  if (!dir.exists(save_dir)) dir.create(save_dir)
  if (is.null(palette)) palette <- c("#999999", "#E69F00")
  group_levels <- c(control_label, case_label)
  
  all_stats_excel <- data.frame()
  combined_data_long <- data.frame()
  bracket_lines_data <- data.frame()
  bracket_text_data <- data.frame()
  
  # Themes
  theme_combined <- theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      legend.position = "right", 
      legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5),
      axis.line = element_blank()
    )
  
  theme_individual <- theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(face = "bold", color = "black"),
      axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5)
    )
  
  # --- DYNAMIC ANALYSIS LOOP ---
  phase_names <- names(phase_list)
  phase_counter <- 1 
  
  for (p_name in phase_names) {
    cols <- phase_list[[p_name]]
    ctrl_col <- cols[1]
    case_col <- cols[2]
    
    if (!ctrl_col %in% names(dataset) || !case_col %in% names(dataset)) {
      warning(paste("Columns missing for phase:", p_name, "- Skipping."))
      next
    }
    
    c_vals <- as.numeric(na.omit(dataset[[ctrl_col]]))
    t_vals <- as.numeric(na.omit(dataset[[case_col]]))
    
    if (length(c_vals) < 2 || length(t_vals) < 2) next
    
    # Stats
    t_res <- t.test(t_vals, c_vals)
    p_val <- t_res$p.value
    
    # ANOVA (Optional check)
    long_tmp <- data.frame(V=c(c_vals, t_vals), G=factor(c(rep("C",length(c_vals)), rep("T",length(t_vals)))))
    aov_p <- summary(aov(V ~ G, data = long_tmp))[[1]][["Pr(>F)"]][1]
    
    sig <- "ns"
    if (p_val < 0.001) { sig <- "***" } else if (p_val < 0.01) { sig <- "**" } else if (p_val < 0.05) { sig <- "*" }
    
    p_plot_text <- format(round(p_val, 6), scientific = FALSE)
    
    # Mean and SD
    mean_c <- mean(c_vals); sd_c <- sd(c_vals)
    mean_t <- mean(t_vals); sd_t <- sd(t_vals)
    c_stat_str <- paste0(round(mean_c, 2), " \u00B1 ", round(sd_c, 2))
    t_stat_str <- paste0(round(mean_t, 2), " \u00B1 ", round(sd_t, 2))
    
    # Excel Report
    all_stats_excel <- rbind(all_stats_excel, data.frame(
      Phase = p_name,
      Control_Stats = c_stat_str,
      Case_Stats    = t_stat_str,
      T_Test_Pvalue = p_val,
      ANOVA_Pvalue  = aov_p,
      Significance  = sig
    ))
    
    # Individual Plot logic remains same but uses p_name
    indiv_df <- data.frame(
      Group = factor(c(rep(control_label, length(c_vals)), rep(case_label, length(t_vals))),
                     levels = group_levels),
      Value = c(c_vals, t_vals)
    )
    stat_txt <- paste0("P = ", p_plot_text, " ", sig)
    
    p_ind <- ggplot(indiv_df, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
      geom_jitter(aes(fill = Group), shape = 21, color = "black", stroke = 0.8,
                  width = 0.1, size = 3, alpha = 0.9, show.legend = FALSE) +
      scale_fill_manual(values = palette) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
      labs(title = p_name, y = y_axis_title) + 
      theme_individual +
      annotate("text", x = 1.5, y = Inf, label = stat_txt, 
               vjust = 1.5, size = 4, color = "gray40", fontface = "bold.italic")
    
    ggsave(file.path(save_dir, paste0("Individual_", p_name, ".tiff")), p_ind, width=4, height=5, dpi=300)
    
    # Combined Data
    combined_data_long <- rbind(combined_data_long, data.frame(
      Phase = p_name, 
      Group = factor(c(rep(control_label, length(c_vals)), rep(case_label, length(t_vals))),
                     levels = group_levels),
      Value = c(c_vals, t_vals)
    ))
    
    # Bracket Calculation
    max_val <- max(c(c_vals, t_vals))
    y_bar <- max_val * 1.05 
    y_leg_c <- max(c_vals) + (max_val * 0.02)
    y_leg_t <- max(t_vals) + (max_val * 0.02)
    
    path_coords <- data.frame(
      x = c(phase_counter - 0.2, phase_counter - 0.2, phase_counter + 0.2, phase_counter + 0.2),
      y = c(y_leg_c, y_bar, y_bar, y_leg_t),
      Bracket_ID = paste0("Bracket_", phase_counter)
    )
    bracket_lines_data <- rbind(bracket_lines_data, path_coords)
    
    bracket_text_data <- rbind(bracket_text_data, data.frame(
      x = phase_counter,
      y = y_bar,
      label = paste0(sig, "\n(P=", p_plot_text, ")") 
    ))
    
    phase_counter <- phase_counter + 1
  }
  
  if (nrow(all_stats_excel) > 0) write_xlsx(all_stats_excel, file.path(save_dir, "Full_Stats_Report.xlsx"))
  
  # --- PLOTTING COMBINED ---
  combined_data_long$Phase <- factor(combined_data_long$Phase, levels = phase_names)
  
  p_combined <- ggplot(combined_data_long, aes(x = Phase, y = Value, fill = Group)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_point(aes(fill = Group), shape = 21, color = "black", stroke = 0.6,
               position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
               size = 2.5, alpha = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = palette) +
    labs(x = x_axis_title, y = y_axis_title) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) + 
    theme_combined +
    geom_path(data = bracket_lines_data, 
              aes(x = x, y = y, group = Bracket_ID), 
              inherit.aes = FALSE, 
              color = "black", linewidth = 0.6, lineend = "round", linejoin = "round") +
    geom_text(data = bracket_text_data, 
              aes(x = x, y = y, label = label), 
              inherit.aes = FALSE, vjust = -0.2, 
              size = 3.5, color = "gray30", fontface = "bold.italic")
  
  ggsave(file.path(save_dir, "Combined_Plot_Final.tiff"), p_combined, width=2 + (2 * length(phase_names)), height=6, dpi=600, compression = "lzw")
  
  message(paste("Done! Results saved in:", save_dir))
  print(p_combined)
}




# ============================================
# EXAMPLE USAGE
# ============================================

P53 = read_excel("filename.xlsx")

# Define the phase mapping in a list
# Structure: "Display Label" = c("Control_Column", "Case_Column")
P53_phases <- list(
  "S"    = c("S-Control", "S-Si"),
  "G2/M" = c("G2-Control", "G2-Si")
)
P53_phases2 <- list(
  "G1" = c("G1-Control", "G1-Si"),
  "S"    = c("S-Control", "S-Si"),
  "G2/M" = c("G2-Control", "G2-Si")
)

analyze_cell_cycle(
  input_data    = P53,
  phase_list    = P53_phases2,
  control_label = "Control",
  case_label    = "Si-Treated",
  x_axis_title  = "",
  palette       = c("#512DA8", "#E91E63"),
  save_dir      = "path"
)



