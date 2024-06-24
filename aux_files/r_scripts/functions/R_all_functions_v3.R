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

  read.delim(fp, header = FALSE) %>%
    mutate(V1 = as.character(V1)) %>%
    mutate(col_names = case_when(grepl("^\\[", V1) & grepl("\\]", V1) ~ V1,
                                 TRUE ~ NA_character_)) %>%
    fill(col_names, .direction = "down") %>%
    filter(col_names != V1) %>%
    mutate(col_names = gsub(",*$", "", col_names)) %>%
    mutate(col_names = gsub("\\[|\\]", "", col_names)) %>%
    pivot_wider(names_from = "col_names", values_from = "V1", values_fn = list)

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
  if (row_num > 14) {
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
    filter(!grepl(" - No Variant Calls$| - Poor$", sample_type)) %>%

    #filter out controls
    filter(!isControl) %>%

    #filter if no coverage
    filter(!is.na(median_coverage)) %>%

    #filter out bad results
    filter(if_any(contains(c("nc_qc_status", "nextclade_qc_overallstatus")), ~ .x != "bad")) %>%

    #filter out samples with less than 80% coverage
    filter(pct_genome_coverage_over_30x >= .80) %>%

    #filter out samples that didn't have a pangolin report
    filter(if_any(contains(c("pango_qc_status", "pangolin_qc_status")), ~ .x != "fail"))

}

## ============================================================
##   Resubmit the shell script if proxy authentication required
## ============================================================

cli_submit <- function(exe_path, bs_cli_command, sh_arguments, shQuote_type = "sh") {

  n_tries <- 5
  time2sleep <- 10
  success <- FALSE

  for(i in 1:n_tries) {
    Sys.sleep(time2sleep)
    message(paste0("\nTrying ", basename(bs_cli_command), " ", i, " time"))
    shell_return <- system2(exe_path,                               #the git shell executable on windows
                            args = c(                               #some arguments need to be in quotes to be passed to shell
                              shQuote(bs_cli_command, type = shQuote_type), #the shell script to run
                              sh_arguments                          #arguments
                            ), stdout = TRUE, stderr = TRUE)
    message(shell_return)

    if(!any(grepl("Proxy Authentication Required|502 Bad Gateway: Reason unknown", shell_return))) {
      success <- TRUE # set i to 0 so it doesn't trigger the error if it retried 5 times and was successful on the 5th time
      break
    }

  }


  if(success) {
    return(shell_return)
  } else{
    stop(simpleError(paste0("Retried ", n_tries, " time")))
  }

}

## =============================================================
##   Read in Excel file with sheet name. Return NULL if no sheet
## =============================================================

read_excel_safely <- function(file, sheet, skip_row = 0) {

  try_read_excel <- purrr::safely(readxl::read_excel, otherwise = NULL)

  try_read_excel(file, sheet = sheet, skip = skip_row)$result

}

## ===================================================================================
##   Run the following commands through the RStudio terminal (mainly for ssh commands)
## ===================================================================================

run_in_terminal <- function(command2run = "", command2print = "") {

  init_terminal <- rstudioapi::terminalExecute(command = command2run)

  # sleep if command hasn't finished running
  while (is.null(rstudioapi::terminalExitCode(init_terminal))) {
    Sys.sleep(1)
  }
  # throw error for non-zero exit codes
  if(rstudioapi::terminalExitCode(init_terminal) != 0) {
    rstudioapi::executeCommand('activateConsole')
    stop(simpleError(paste0("\nThere was an issue running the ssh command through the terminal!\n",
                            "Run the following command through the EC2 instance on AWS or follow the instructions:\n\n",
                            command2print)))
  }

  Sys.sleep(5)
  rstudioapi::terminalKill(init_terminal)
}

## =============================================================================
##   Use the run_in_terminal function to submit jobs to EC2 instance through ssh
## =============================================================================

submit_screen_job <- function(message2display = "Running function to submit screen job", ec2_login = "", screen_session_name = "", command2run = "") {

  run_in_terminal(paste("echo", paste0("'", message2display, " through a screen command on EC2 instance [", ec2_login, "]';"),
                        "ssh -tt", ec2_login,
                        shQuote(paste("sleep 5;",
                                      "if screen -ls | grep", paste0("'", screen_session_name, "'"), "-q;",
                                        "then echo 'Detached", screen_session_name, "session detected. Skipping ahead to check status';",
                                        "else echo 'Submitting job now';",
                                        "mkdir -p ~/.tmp_screen;",
                                        "rm -f", paste0("~/.tmp_screen/", screen_session_name, ".screenlog;"),
                                        "sleep 5;",
                                        "screen -dm -S", screen_session_name, "-L -Logfile", paste0("~/.tmp_screen/", screen_session_name, ".screenlog"), "bash -c",
                                        paste0("\"", command2run, "\";"),
                                        "sleep 5;",
                                      "fi"), type = "sh")),
                  command2run)
}

## =============================================================
##   Use the run_in_terminal function to check on submitted jobs
## =============================================================

check_screen_job <- function(message2display = "Running function to check screen job", ec2_login = "", screen_session_name = "") {

  run_in_terminal(paste("echo 'Checking for dead jobs';",
                        "ssh -tt", ec2_login,
                        shQuote(paste("while (screen -ls | grep 'Dead' -q);",
                                      "do echo 'ERROR DETECTED!!! The previous submitted job is dead';",
                                      "echo 'Close this terminal and rerun the whole script or just the previous job';",
                                      "echo 'This message will remain here indefinitely until the RStudio terminal is closed';",
                                      "screen -wipe;",
                                      "sleep infinity;",
                                      "done;"), type = "sh")))

  # monitor the screen log file
  run_in_terminal(paste("echo", paste0("'", message2display, "';"),
                        "ssh -tt", ec2_login,
                        shQuote(paste("sleep 5;",
                                      "SCREEN_PID=`screen -ls | grep", screen_session_name, "| cut -f1 -d'.' | sed 's/\\W//g'`;",
                                      "if test -z ${SCREEN_PID};",
                                      "then echo 'Job finished already';",
                                      "echo -n 'Log last modified: ';",
                                      "TZ='US/Eastern'",
                                      "date '+%F %r' -r", paste0("~/.tmp_screen/", screen_session_name, ".screenlog;"),
                                      "echo '\n';",
                                      "echo", paste0(c(rep("-", 100), ";"), collapse = ""),
                                      "echo '\n';",
                                      "tail -n 100", paste0("~/.tmp_screen/", screen_session_name, ".screenlog;"),
                                      "sleep 15;",
                                      "else echo Monitoring screen session: $SCREEN_PID;",
                                      "echo '\n';",
                                      "echo", paste0(c(rep("-", 100), ";"), collapse = ""),
                                      "echo '\n';",
                                      "tail --pid=$SCREEN_PID -f", paste0("~/.tmp_screen/", screen_session_name, ".screenlog;"),
                                      "sleep 15;",
                                      "fi"), type = "sh")))
}

## =============================================
##  Convert sample size to font size (inversely)
## =============================================

convert_sample_size_2_font_size <- function(sample_size = x,
                                            max_sample_size = 96,
                                            min_font = 0.5,
                                            max_font = 8) {

  if(sample_size > max_sample_size) {
    stop(simpleError(paste("Sample size input of", sample_size, "is greater than the max sample size of", max_sample_size,
                           "\nPlease adjust these inputs for this function accordingly")))
  }

  font_range = min_font - max_font

  ceiling((((sample_size*font_range)/max_sample_size) + max_font)*2)/2

}

## ==============================
##  Reverse complement DNA string
## ==============================

reverse_complement <- function(index) {

  non_atcg <- gsub("[ATCG]", "", index, ignore.case = TRUE)

  if(all(nchar(non_atcg) != 0)) {
    stop(simpleError("There are non-canonical bases in your indices"))
  }

  reversed_index <- stringi::stri_reverse(index)
  chartr(old = "atcgATCG", new = "tagcTAGC", reversed_index)
}

## =================================================================
##  Filter the dataframe and pivot wider with specified column names
## =================================================================

filter_n_pivot <- function(df, db_field, db_name, col_name, col_value) {

  db_field <- enquo(db_field)

  df %>%
    filter(!!db_field == db_name) %>%
    select(col_name, col_value) %>%
    pivot_wider(names_from = col_name, values_from = col_value)
}

## ======================================================================================
##  Calculate final concentration after accounting for dilution and concentration factors
## ======================================================================================

stock_concentration_calc = function(concentration, rxn_vol, input_vol, extraction_vol, elution_vol, end_concentration_vol, initial_concentration_vol, efficiency) {
  result = (concentration * rxn_vol / input_vol) * (elution_vol / extraction_vol) * (end_concentration_vol / initial_concentration_vol) / efficiency
  return(result)
}
