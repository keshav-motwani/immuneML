pdf(NULL)

plot_heatmap = function(matrix,
                        row_annotations = NULL,
                        column_annotations = NULL,
                        one_hot_row_annotations = NULL,
                        one_hot_column_annotations = NULL,
                        split_rows = NULL,
                        split_columns = NULL,
                        palette = list(),
                        row_names = NULL,
                        column_names = NULL,
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        show_row_dend = TRUE,
                        show_column_dend = TRUE,
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        show_legend_row,
                        show_legend_column,
                        legend_position = "side",
                        lower_quantile = 0.01,
                        upper_quantile = 0.99,
                        text_size = 10,
                        row_names_size = 4,
                        column_names_size = 4,
                        value_name = "value",
                        title = character(0),
                        row_standardization = "scale",
                        heatmap_color = "BWR",
                        color_by_quantile = TRUE,
                        height = 10,
                        width = 10,
                        result_name = "test_heatmap",
                        result_path = getwd()) {

  params = as.list(match.call())
  params[[1]] = NULL
  saveRDS(params, file.path(result_path, paste0(result_name, ".rds")))

  if (is.character(matrix)) {
    matrix = read.table(matrix, header = FALSE, sep = ",")
  }
  if (!is.null(column_names))
    colnames(matrix) = as.character(column_names)
  if (!is.null(row_names))
    rownames(matrix) = as.character(row_names)
  if (is.character(row_annotations)) {
    row_annotations = as.data.frame(readr::read_csv(row_annotations))
    row_annotations$X1 = NULL
  }
  if (is.character(column_annotations)) {
    column_annotations = as.data.frame(readr::read_csv(column_annotations))
    column_annotations$X1 = NULL
  }
  if (is.character(one_hot_row_annotations)) {
    one_hot_row_annotations = as.data.frame(readr::read_csv(one_hot_row_annotations))
    one_hot_row_annotations$X1 = NULL
  }
  if (is.character(one_hot_column_annotations)) {
    one_hot_column_annotations = as.data.frame(readr::read_csv(one_hot_column_annotations))
    one_hot_column_annotations$X1 = NULL
  }
  if (!is.list(palette))
    palette = rjson::fromJSON(palette)
  show_legend_row = as.character(show_legend_row)
  show_legend_column = as.character(show_legend_column)

  if (row_standardization == "scale") {
    matrix = t(scale(t(matrix)))
  } else if (row_standardization == "min_max") {
    indices = rowSums(matrix) > 0
    matrix = matrix[indices, ]
    row_annotations = row_annotations[indices, , drop = FALSE]
    one_hot_row_annotations = one_hot_row_annotations[indices, , drop = FALSE]
    matrix = t(apply(
      t(matrix),
      MARGIN = 2,
      FUN = function(X)
        (X - min(X)) / diff(range(X))
    ))
  }

  heatmap = ggexp::plot_heatmap(
    matrix = matrix,
    row_annotations = row_annotations,
    column_annotations = column_annotations,
    one_hot_row_annotations = one_hot_row_annotations,
    one_hot_column_annotations = one_hot_column_annotations,
    split_rows = split_rows,
    split_columns = split_columns,
    palette = palette,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_dend = show_row_dend,
    show_column_dend = show_column_dend,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    show_legend_row = show_legend_row,
    show_legend_column = show_legend_column,
    legend_position = legend_position,
    lower_quantile = lower_quantile,
    upper_quantile = upper_quantile,
    text_size = text_size,
    row_names_size = row_names_size,
    column_names_size = column_names_size,
    value_name = value_name,
    title = title,
    heatmap_color = heatmap_color,
    color_by_quantile = color_by_quantile
  )

  png(
    file = file.path(result_path, paste0(result_name, ".png")),
    height = height,
    width = width,
    units = "in",
    res = 1200
  )

  print(heatmap)

  dev.off()

  return(heatmap)
}
