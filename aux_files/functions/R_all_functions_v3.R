## ======================================================================
##   Order the factors of one column based on the other of another column
## ======================================================================

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

## ========================================================
##   Read in SampleSheet.csv file with loop. Returns a list
## ========================================================

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

## =========================================
##   Function for turning tables into kables
## =========================================

kable_style <- function(data, font_sizing = 10) {

  row_num <- nrow(data)

  ##substitute underscore with escaped underscores and remove na in p.value columns
  data_return <- data %>%
    #replace underscores with escaped underscores in column names
    select_all(~gsub("_", "\\\\_", .)) %>% ##need to escape the escape
    #change underscores to escaped underscores
    mutate_if(function(x) is.character(x) | is.factor(x), ~gsub("_", "\\\\_", .)) %>%
    #escape percent signs
    mutate_if(function(x) is.character(x) | is.factor(x), ~gsub("%", "\\\\%", .))

  # ... should be column number
  if (row_num > 40) {
    data_return <- data_return %>%
      kable("latex", longtable = T, digits=2, booktabs=T, escape=F) %>%
      kable_styling(latex_options = c("repeat_header", "HOLD_position"), font_size = font_sizing) %>%
      row_spec(0, bold = T, color="#7C0A02")

  }
  else {
    data_return <- data_return %>%
      kable("latex", longtable = F, digits=2, booktabs=T, escape=F) %>%
      kable_styling(latex_options = c("scale_down", "repeat_header", "HOLD_position")) %>%
      row_spec(0, bold = T, color="#7C0A02")

    if(row_num > 1) { ##always collapse row unless there is only 1 row
      data_return <- data_return %>%
        collapse_rows(columns = 1, valign = "top")
    }
  }

  return(data_return)

}

## ==============================
##   Filter samples for reporting
## ==============================

filter4report <- function(data) {

  data %>%

    #filter out unassigned reads sample
    filter(!grepl("None|Undetermined", sample_id)) %>%

    #filter out samples with any pangolin data; this should be the same as filtering out is.na(lineage)
    filter(!grepl(" - No Variant Calls$", sample_type)) %>%

    #filter out controls
    filter(!isControl) %>%

    #filter out bad results
    filter(nc_qc_status != "bad") %>%

    #filter if no coverage
    filter(!is.na(median_coverage)) %>%

    #filter out mediocre results with lower coverages
    filter(!(nc_qc_status == "mediocre" & pct_genome_coverage_over_30x < .8)) %>%

    #filter out samples that didn't have a pangolin report
    filter(pango_qc_status != "fail")

}

## ============================================================
##   Resubmit the shell script if proxy authentication required
## ============================================================

cli_submit <- function(exe_path, cli_command, sh_arguments, shQuote_type = "sh") {

  for(i in 1:5) {
    Sys.sleep(30)
    message(paste0("Trying ", basename(cli_command), " ", i, " time"))
    shell_return <- system2(exe_path,                               #the git shell executable on windows
                            args = c(                               #some arguments need to be in quotes to be passed to shell
                              shQuote(cli_command, type = shQuote_type), #the shell script to run
                              sh_arguments                          #arguments
                            ), stdout = TRUE)
    message(shell_return)

    if(!any(grepl("Proxy Authentication Required|502 Bad Gateway: Reason unknown", shell_return))) {
      break
    }

  }


  if(i == 5) {
    stop(simpleError("Retried 5 times"))
  } else{
    return(shell_return)
  }

}

