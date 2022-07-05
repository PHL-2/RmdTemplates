## ===========================================
##   Order the factors of one column based on the other of another column
## ===========================================

##this function is the same as:
##mutate(col2order = factor(col2order, levels = unique(col2order[order(sortingcol, decreasing = TRUE)])))

order_on_other_col <- function(df, col2order, sortingcol, decreasing = TRUE) {
  
  col2order <- enquo(col2order)
  sorting_name <- gsub("!!", "", deparse(substitute(!!sortingcol)))
  ordering <- order(unlist(df[, sorting_name]), decreasing = decreasing)
  uniq <- unique(unlist(df[,gsub("!!~", "", deparse(substitute(!!col2order)))])[ordering])
  
  df %>%
    mutate(!!quo_name(col2order) := factor(!!col2order,  levels = uniq))
  
}

## ===========================================
##   Read in SampleSheet.csv file with loop. Returns a list
## ===========================================

load_sample_sheet <- function(fp) {
  
  read_delim(fp, delim = "no_delim", col_names = FALSE) %>%
    mutate(col_names = case_when(grepl("^\\[", X1) & grepl("\\]$", X1) ~ X1,
                                 TRUE ~ NA_character_)) %>%
    fill(col_names, .direction = "down") %>%
    filter(col_names != X1) %>%
    mutate(col_names = gsub("\\[|\\]", "", col_names)) %>%
    pivot_wider(names_from = "col_names", values_from = "X1", values_fn = list)
  
}

## ========================
##   End line with new line
## ========================

start_w_newline <- function(string) {
  gsub("^", "\n", string)
}
