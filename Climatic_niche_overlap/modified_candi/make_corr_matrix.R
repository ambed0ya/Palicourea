#' make corr matrix
#'
#' A functiont that creates a correlation matrix between the environmental data variables of correlation within the area of the maximum extents of occurrence points
#'
#' SpatRaster inputs (e.g., WorldClim from {geodata}).
#'
#' @param occurrences A list of data frames with occurrence data for each species.
#'        Must contain columns named 'latitude' and 'longitude'.
#' @param environment_data A terra SpatRaster (preferred) or a raster Raster* object.
#' @return a pdf of the correlation matrix is saved to the local disk
make_corr_matrix <- function(
  occurrences,
  environment_data,
  out_dir = getwd(),
  abs_highlight = 0.8
) {

  # ---- deps (fail early with helpful message) ----
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required. Install it with install.packages('terra').")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
  }

  # ---- normalize occurrences ----
  if (is.data.frame(occurrences)) {
    occ_df <- occurrences
  } else {
    occ_df <- do.call(rbind, occurrences)
  }
  if (!all(c("latitude", "longitude") %in% names(occ_df))) {
    stop("occurrences must contain columns named 'latitude' and 'longitude'.")
  }

  # ---- extent from occurrences ----
  max_lat <- max(occ_df$latitude, na.rm = TRUE)
  min_lat <- min(occ_df$latitude, na.rm = TRUE)
  max_lon <- max(occ_df$longitude, na.rm = TRUE)
  min_lon <- min(occ_df$longitude, na.rm = TRUE)

  # ---- coerce rasters to terra ----
  if (inherits(environment_data, "Raster")) {
    environment_data <- terra::rast(environment_data)
  }
  if (!inherits(environment_data, "SpatRaster")) {
    stop("environment_data must be a terra SpatRaster (or a raster Raster* object that can be converted).")
  }

  # ---- crop ----
  crop_ext <- terra::ext(min_lon, max_lon, min_lat, max_lat)
  cropped <- terra::crop(environment_data, crop_ext)

  # ---- correlation matrix across layers ----
  message("Calculating correlation coefficients of climatic variables")
  # terra::layerCor expects argument name `fun` and (for Pearson correlation)
  # the built-in option is "cor" (not "pearson").
  # Using an unsupported string can lead to confusing downstream errors.
  cor_res <- terra::layerCor(cropped, fun = "cor", use = "complete.obs")

  # layerCor returns either a matrix or a list depending on 'fun'
  # (for standard stats it can return a list with $cor)
  if (is.list(cor_res) && !is.null(cor_res$cor)) {
    c_matrix <- cor_res$cor
  } else {
    c_matrix <- cor_res
  }

  # ---- output dir ----
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- save matrix ----
  csv_name <- paste0("correlationBioClim_", Sys.Date(), ".csv")
  csv_path <- file.path(out_dir, csv_name)
  message("Saving correlation matrix to: ", csv_path)
  utils::write.csv(c_matrix, file = csv_path)

  # ---- bubble-style plot (lower triangle) ----
  # Goal: match your reference style:
  # - show LOWER triangle only (incl. diagonal)
  # - print r values in every cell
  # - draw GREY bubbles ONLY where |r| >= abs_highlight (e.g., 0.8)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" is required. Install it with install.packages(\"ggplot2\").")
  }

  vars <- colnames(c_matrix)
  if (is.null(vars)) {
    vars <- paste0("layer", seq_len(ncol(c_matrix)))
    colnames(c_matrix) <- vars
    rownames(c_matrix) <- vars
  }

  # Long-format table of correlations
  cor_df <- as.data.frame(as.table(c_matrix), stringsAsFactors = FALSE)
  names(cor_df) <- c("Var1", "Var2", "cor")
  cor_df$cor <- as.numeric(cor_df$cor)

  # Factor ordering for axes
  cor_df$Var1 <- factor(cor_df$Var1, levels = vars)
  cor_df$Var2 <- factor(cor_df$Var2, levels = vars)

  # Keep lower triangle (including diagonal): Var1 index >= Var2 index
  cor_df <- cor_df[as.integer(cor_df$Var1) >= as.integer(cor_df$Var2), , drop = FALSE]

  # Reverse y-axis order so the kept half *looks* like the lower triangle
  cor_df$Var2 <- factor(cor_df$Var2, levels = rev(levels(cor_df$Var2)))

  # Threshold used for grey bubbles (explicit, so it's easy to see/change)
  threshold <- abs_highlight

  # Labels and highlight mask
  cor_df$lab <- sprintf("%.2f", cor_df$cor)
  cor_df$high_corr <- abs(cor_df$cor) >= threshold

  plot_obj <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Var1, y = Var2)) +
    # Grey bubbles ONLY for |r| >= threshold
    ggplot2::geom_point(
      data = cor_df[cor_df$high_corr, , drop = FALSE],
      ggplot2::aes(size = abs(cor)),
      color = "grey40"
    ) +
    # Correlation values for all kept cells
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 3) +
    ggplot2::scale_size(range = c(4, 10), guide = "none") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0),
      axis.text.y = ggplot2::element_text(),
      panel.grid = ggplot2::element_blank()
    )

  pdf_name <- paste0("correlationsBioClim_", Sys.Date(), ".pdf")
  pdf_path <- file.path(out_dir, pdf_name)
  message("Saving PDF of correlation plot to: ", pdf_path)
  message("(For reference: |correlation| >= ", threshold, " is often treated as 'high'.)")

  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot_obj,
    device = "pdf",
    width = 11,
    height = 8.5
  )

  return(c_matrix)
}
