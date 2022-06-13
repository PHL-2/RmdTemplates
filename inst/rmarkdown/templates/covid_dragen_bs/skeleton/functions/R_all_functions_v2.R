## ================
## Get alpha diversity measures
## ================

get_alpha <- function(toTest, counts_matrix, colMerge = "SampleID") {
  toTest %>%  
    merge(diversity(t(counts_matrix)), by.x=colMerge, by.y="row.names", all.x=T) %>%
    dplyr::rename(shannon = y) %>%
    merge(rarefy(t(counts_matrix), 1000), by.x=colMerge, by.y="row.names", all.x=T) %>%
    dplyr::rename(richness = y) %>%
    gather("metric", "alpha", c("richness", "shannon")) %>%
    mutate(metric = fct_recode(metric, Richness="richness", Shannon="shannon"))
}

## ================
##   Run linear model
## ================

run_lm <- function(props_toTest, s_toTest, form1, p_cutoff) {
  set.seed(2)
  props_toTest[,s_toTest$SampleID] %>%
    as.data.frame() %>%
    rownames_to_column(var = "Taxa") %>%
    gather(key = "SampleID", value = "Prop", s_toPlot$SampleID) %>%
    fix_zeros_v2(Prop) %>%
    merge(s_toTest, by="SampleID") %>%
    group_by(Taxa) %>%
    do(tidy_lm(lm(as.formula(form1), data=., na.action=na.omit))) %>%
    setNames(c("Taxa","term","Estimate","Std.Error","t.value","p.value")) %>%
    ungroup() #%>%
  #filter(term != '(Intercept)') #%>%
  #group_by(term) %>%
  #mutate(fdr = p.adjust(p.value, method="BH")) %>%
  #ungroup() %>%
  #filter(p.value < p_cutoff)
}

## ================
##   Tidy lm results
## ================

tidy_lm <- function(lm_test) {
  set.seed(2)
  mod <- summary(lm_test)
  data.frame(term  = rownames(mod$coefficients), mod$coefficients, row.names=NULL)
}

## ================
##   Add a small value to zero's
## ================

fix_zeros_v2 <- function(props_long, props_label = props) {
  props_label <- enquo(props_label)
  props_long %>%
    mutate(!!props_label := !!props_label + min(filter(props_long, !!props_label > 0) %>% pull(!!props_label)) / 10) %>%
    mutate(props_log := log(!!props_label))
}

## ================
##   Rearranging date format
## ================

change_date_format <- function(d) { #change date format to MM-DD-YY
  if(grepl("-", d)) {
    paste(substr(d,6,7), substr(d,9,10), substr(d,1,4), sep="-")
  }
  else if (grepl("/", d)) {
    gsub("/", "-", d)
  } 
  else if (nchar(unique(d)) == 8) {
    paste(substr(d,5,6), substr(d,7,8), substr(d,1,4), sep="-")
  }
  else {
    stop (simpleError(paste0("Your date ", d, " is not in YYYY-MM-DD, MM/DD/YY, or MMDDYYYY format.")))
  }
}

## ================
##   Summarizing model functions
## ================

tidy_lmer <- function(lmer_test) { ##get summary of table from lmer
  set.seed(2)
  mod <- summary(lmer_test)
  data.frame(term  = rownames(mod$tTable), mod$tTable, row.names=NULL) #%>%
  #mutate(AIC = AIC(lmer_test))
}
tidy_lmer2 <- function(lmer_test, term_string) { ##get summary of table from lmer when using emmeans
  set.seed(2)
  mod <- anova(lmer_test)
  form_pairwise <- as.formula(paste0("pairwise~", term_string))
  bind_rows(data.frame(contrast = rownames(mod), mod, row.names=NULL),
            data.frame(emmeans(lmer_test, list(form_pairwise), adjust="tukey")[[2]])) %>%
    #mutate(AIC = AIC(lmer_test)) %>%
    select(contrast, p.value, everything()) %>%
    mutate(estimate = -estimate)
}

## ================
##   Summarize just the p-value 
## ================

tidy_anova_lme <- function (lme_model) {
  set.seed(2)
  anova(lme_model) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    select(-denDF) %>%
    rename(df = numDF) %>%
    filter(!(term %in% "(Intercept)")) %>%
    rename(`p.value`='p-value')
}

## ================
##   PERMANOVA functions
## ================

tidy_permanova <- function(anov){ ##get table of permanova results
  set.seed(2)
  data.frame(Term = rownames(anov$aov.tab), anov$aov.tab, row.names = NULL) %>%
    rename(p.value = Pr..F.)
}

permanova_with_shuffle_1_group_posthoc <- function(dist_matrix, s_toTest, group_label, rep_mes_label, covariates, perm, is_within, p_cutoff=0.05){
  set.seed(2)
  s_toTest <- data.frame(s_toTest)
  a_ixn <- permanova_with_shuffle_1_group(dist_matrix, s_toTest, group_label, rep_mes_label, covariates, perm=perm, is_within)
  #print(class(s_toTest))
  combs <- combn(unique(s_toTest[[group_label]]), 2)
  num_tests <- dim(combs)[2]
  # do post hoc tests
  if (filter(a_ixn, Term==group_label)$p.value < p_cutoff) {
    post_hocs <- lapply(1:num_tests,
                        function(x) data.frame(comparison = paste(combs[,x], collapse=' - '),
                                               permanova_with_shuffle_1_group(dist_matrix, s_toTest[s_toTest[[group_label]] %in% combs[,x],], group_label, rep_mes_label, covariates, perm=perm, is_within)))
    a_ixn <- rbind(data.frame(comparison="all", a_ixn), do.call(rbind, post_hocs))
  }
  a_ixn
}

##if perm is 100, the lowest p-value you can get per test is 1/100
##if perm is 20, the lowest p-value you can get per test is 0.05
##the higher the perm, the more significance you can get from the test
permanova_with_shuffle_1_group <- function(dist_matrix, s_toTest, group_label, rep_mes_label, covariates, perm, is_within=F){
  set.seed(2)
  s_toTest <- data.frame(s_toTest)
  dist_toTest <- usedist::dist_subset(dist_matrix, s_toTest$SampleID)
  form1 <- paste("dist_toTest", "~", group_label)
  if (!is.na(covariates)) {
    form1 <- paste(form1, " + ", covariates)
  }
  a_ixn <- adonis(as.formula(form1), data=s_toTest, permutations=perm)
  f_ixn <- a_ixn$aov.tab[1, 4]
  
  fs_permuted <- replicate(perm, {
    s_permuted <- s_toTest
    if (is_within){
      s_permuted[,group_label] <- shuffle_within_groups(s_permuted[,group_label], s_permuted[,rep_mes_label])
    } else {
      s_permuted[,group_label] <- shuffle_between_groups(s_permuted[,group_label], s_permuted[,rep_mes_label])
    }
    a_permuted <- adonis(as.formula(form1), s_permuted, permutations = 4)
    a_permuted$aov.tab[1, 4]
  })
  p_ixn <- sum(c(f_ixn, fs_permuted) >= f_ixn) / (length(fs_permuted) + 1)
  a_ixn$aov.tab[1,6] <- p_ixn
  tidy_permanova(a_ixn)
}

shuffle_within_groups <- function(x,g) { ##shuffle x within each group in g
  ave(x, g, FUN = function(a) if(length(a)>1) sample(a) else a)
}

permanova_with_shuffle_2_groups <- function(dist_matrix, s_toTest, group_label1, group_label2, rep_mes_label, covariates, perm, first_within=F, second_within=F){
  set.seed(2)
  s_toTest <- as.data.frame(s_toTest)
  dist_toTest <- usedist::dist_subset(dist_matrix, s_toTest$SampleID)
  form1 <- paste("dist_toTest", "~", group_label1, " * ", group_label2)
  if (!is.na(covariates)) {
    form1 <- paste(form1, " + ", covariates)
  }
  a_ixn_orj <- adonis(as.formula(form1), data=s_toTest, permutations=perm)
  
  terms_perm <- c(group_label1, group_label2, paste0(group_label1, ":", group_label2))
  
  tidy_output <- tidy_permanova(a_ixn_orj)
  f_ixn_all <- tidy_output[match(terms_perm, tidy_output$Term),"F.Model"]
  #select(Term, F.Model)
  
  fs_permuted <- replicate(perm, {
    s_permuted <- s_toTest
    
    if (first_within) {
      s_permuted[,group_label1] <- shuffle_within_groups(s_permuted[,group_label1], s_permuted[,rep_mes_label])
    } else {
      s_permuted[,group_label1] <- shuffle_between_groups(s_permuted[,group_label1], s_permuted[,rep_mes_label])
    }
    
    if (second_within) {
      s_permuted[,group_label2] <- shuffle_within_groups(s_permuted[,group_label2], s_permuted[,rep_mes_label])
    } else {
      s_permuted[,group_label2] <- shuffle_between_groups(s_permuted[,group_label2], s_permuted[,rep_mes_label])
    }
    
    a_permuted <- adonis(as.formula(form1), s_permuted, permutations = 4)
    
    temp_output <- tidy_permanova(a_permuted)
    temp_output[match(terms_perm, temp_output$Term),"F.Model"]
    #c(a_permuted_g1$aov.tab[1, 4], a_permuted_g2$aov.tab[1, 4], a_permuted$aov.tab[3, 4])
  })
  
  p_ixn <- rowSums(cbind(f_ixn_all, fs_permuted) >= f_ixn_all, na.rm = T) / (dim(fs_permuted)[2] + 1)
  
  tidy_output[match(terms_perm, tidy_output$Term),"p.value"] <- p_ixn
  tidy_output  
}

## ================
##   Filter low coverage
## ================
filter_low_coverage <- function(props, frac_cutoff=0.6, min_ab=0){
  frac_nonzero <- function (x) sum(x > min_ab) / length(x)
  apply(props, 1, frac_nonzero) >= frac_cutoff
}
###=====
###  make_pcoa_plot <- function(uu, s, shape_by, color_by, title)
###  uu: distance, s: mapping file, shape_by: variable used for shape, color_by: variable used for color
###=====
make_pcoa_plot <- function(dm, s, shape_by, color_by) {
  dm <- usedist::dist_subset(dm, s$SampleID)
  pc <- pcoa(dm)
  pc_df <- merge(s, pc$vectors[, 1:3], by.x="SampleID", by.y="row.names")
  pc_pct <- round(pc$values$Relative_eig * 100)
  
  pcoa_plot = ggplot(pc_df, aes(x=Axis.1, y=Axis.2)) +
    theme_bw() +
    scale_shape_discrete(name=sub("_", " ", shape_by)) +
    scale_colour_discrete(name=sub("_", " ", color_by)) +
    labs(
      x=paste0("PCoA axis 1 (", pc_pct[1], "%)"),
      y=paste0("PCoA axis 2 (", pc_pct[2], "%)")
    )
  
  if (is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by))))
  } else if (!is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by)), shape=factor(get(shape_by))))
  } else {
    pcoa_plot <- pcoa_plot + geom_point()
  }
  return(pcoa_plot)
}
###=====
###  heatmap_grouped <- function(genus_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1)
###  option=1: rows_to_keep <- filter_low_coverage(heatmap_props, perc_cutoff=thre) ## taxa found in at least 80% of samples
###  option=2: rows_to_keep <- apply(heatmap_props,1,max) >= 0.01 ## taxa with abundance in any sample exceeding 1%
###=====
heatmap_grouped <- function(summed_props, heatmap_s, group_colors = NULL, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1, prop_cut=0.01, satu_limit=0.4, gaps_col = NULL, max_breaks = 1){
  
  #color = saturated_rainbow(101)
  color = saturated_rainbow(101, saturation_limit=satu_limit)
  breaks = c(0, 1e-10, seq(0.001, max_breaks, length.out = 100))
  
  heatmap_props <- summed_props[,heatmap_s$SampleID]
  
  if (option == 1) {
    rows_to_keep <- filter_low_coverage(heatmap_props, frac_cutoff=thre)
  } else if (option == 2) {
    rows_to_keep <- apply(heatmap_props,1,max) >= prop_cut
  }
  heatmap_props <- heatmap_props[rows_to_keep,]
  
  ## group the SampleIDs
  #heatmap_s %<>% arrange_(.dots=grps)
  heatmap_props <- heatmap_props[, heatmap_s$SampleID]
  
  ## update the annotation
  annc <- heatmap_s[,grps] %>% as.data.frame()
  rownames(annc) <- heatmap_s$SampleID
  ##groups must be factors; the factors must also include any NAs (for controls)
  colnames(annc) <- grps
  
  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = annc, annotation_colors = group_colors, color = color, breaks = breaks, filename = fname,
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8, gaps_col = gaps_col)
  else
    pheatmap(heatmap_props, annotation = annc, annotation_colors = group_colors, color = color, breaks = breaks,
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8, gaps_col = gaps_col)
}

## ================
##   Get FDR values
## ================

get_fdr <- function(data, ...) {
  
  grouping_col <- enquos(...)
  
  data_return <- data
  
  if(sum(grepl("Pr...t..", colnames(data_return))) > 0)  {
    data_return <- data_return %>%
      rename(`p.value` = "Pr...t..")
  }
  
  data_return <- data_return %>%
    filter(!is.na(p.value)) %>%
    group_by(!!!grouping_col) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>% 
    ungroup()
  
}

## ================
##   Get most significant values for kable
## ================

###always have to pass a column to sort by based on lowest fdr (usually 1st column)
kable_sort <- function(data, comparison_col, sort_col = 1) {
  
  data_return <- data %>%
    filter(!is.na(p.value))
  
  ##if a column is passed to calculate the fdr
  if(!missing(comparison_col)) {
    ##get the column with the intercepts
    intercept_col <- enquo(comparison_col)
    
    ##sort df by Variable based on lowest fdr first (usually the first column; this can change)
    ##have to turn it to symbol to use in dplyr
    Variable <- sym(colnames(data_return)[sort_col])
    dep_var <- ""
    ##if dependent variable (column name) is Value (depends on the type of table)
    if(sum(grepl("Value", colnames(data_return))) > 0)  {
      dep_var <- sym("Value")
    }
    else if (sum(grepl("Estimate", colnames(data_return))) > 0){
      dep_var <- sym("Estimate")
    }
    else if (sum(grepl("estimate", colnames(data_return))) > 0){
      dep_var <- sym("estimate")
    }
    else if (sum(grepl("F-value", colnames(data_return))) > 0){
      dep_var <- sym("F-value")
    }
    else {
      dep_var <- sym("R2")
    }
    
    ##initiate another df and sort the first column by the lowest FDR and then the value
    sort_df <- data_return %>%
      ##filter out the intercept in intercept_col/comparison_col
      filter(!grepl("Intercept", !!intercept_col)) %>%
      filter(!is.na(!!dep_var)) %>%
      arrange(fdr, -abs(!!dep_var))
    
    ##if the only signif fdr values are from intercepts
    if(only_intercept(data_return, intercept_col)) {
      
      sort_df <- sort_df %>%
        ##get only the top five hits with lowest fdr
        head(5) %>%
        ##change df into vector
        pull(!!Variable)
      
      ##return dataframe with the 5 lowest FDR
      data_return <- data_return %>%
        filter(!!Variable %in% sort_df)
    }
    else {
      sort_df <- sort_df %>%
        filter(fdr < 0.05) %>%
        pull(!!Variable)
      
      data_return <- data_return %>%
        ##get only the groups with significant fdr values
        filter(!!Variable %in% sort_df) %>%
        #filter(fdr < 0.05) %>%
        #filter(!grepl("Intercept", !!intercept_col)) %>%
        mutate(fdr = cell_spec(signif(fdr, 2), "latex",
                               bold = ifelse((!grepl("\\(", !!intercept_col)&(fdr<0.05)), T, F)))
    }
  }
  ##bold the p-value
  else {
    data_return <- data_return %>%
      mutate(p.value = cell_spec(signif(p.value, 2), "latex", bold = p.value<0.05))
  }
  
  return(data_return)
  
}

## ===========================================
##   Check if table only contains the intercept
## ===========================================

only_intercept <- function(data, grouping_col) {
  
  ##get filtered dataframe when fdr < 0.05
  data <- data %>%
    filter(fdr < 0.05)
  ##get dataframe when it's only the intercept
  intercept_df <- data %>%
    filter(grepl("Intercept", !!grouping_col))
  
  if(nrow(data)==nrow(intercept_df)) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

## ===========================================
##   Function for turning tables to kables
## ===========================================

kable_style <- function(data) {
  
  row_num <- nrow(data)
  
  ##substitute underscore with escaped underscores and remove na in p.value columns
  data_return <- data %>%
    select_all(~gsub("_", "\\\\_", .)) %>% ##need to escape the escape
    select_all(~gsub("#", "\\\\#", .)) %>% ##need to escape the escape
    mutate_if(function(x) is.character(x) | is.factor(x), ~gsub("_", "\\\\_", .)) %>%
    mutate_if(function(x) is.character(x) | is.factor(x), ~gsub("%", "\\\\%", .))
  
  ##if Taxa is a column in the dataframe
  if(sum(grepl("Taxa", colnames(data_return))) > 0)  {
    data_return <- data_return %>%
      mutate(Taxa = gsub("[pcofgs]  ", "", Taxa))
  }
  
  # ... should be column number
  if (row_num > 40) {
    data_return <- data_return %>%
      kable("latex", longtable = T, digits=2, booktabs=T, escape=F) %>%
      kable_styling(latex_options = c("repeat_header", "HOLD_position"), font_size = 7) %>%
      row_spec(0, bold = T, color="#7C0A02") #%>%
      #collapse_rows can help add lines in-between rows
      #collapse_rows(columns = 1, valign = "top") 
    
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

## ===========================================
##   Function for calculating distances in parallel (stolen and modified from usedist package)
## ===========================================

##the maximum number of parallel processes running will be decided by availableCores()
##plan(multiprocess) will typically use multicores and then multiple sessions to parallelize processes
  #but because we are using RStudio environment, multicores is disabled and the process is parallelized by opening multiple background R sessions


para_dist_make <- function (x, distance_fcn, tree = NULL, ...) {
  set.seed(2)
  distance_from_idxs <- function (idxs) {
    i1 <- idxs[1]
    i2 <- idxs[2]
    if(is.null(tree)) {
      distance_fcn(x[i1,], x[i2,], ...)
    }
    else {
      distance_fcn(x[i1,], x[i2,], tree, ...)
    }
  }
  
  size <- nrow(x)
  
  plan(multiprocess) #parallelize on local computer
  
  d <- future_apply(utils::combn(size, 2), 2, distance_from_idxs)
  
  attr(d, "Size") <- size
  xnames <- rownames(x)
  if (!is.null(xnames)) {
    attr(d, "Labels") <- xnames
  }
  attr(d, "Diag") <- FALSE
  attr(d, "Upper") <- FALSE
  class(d) <- "dist"
  d
}

## ===========================================
##   Make directory of desired file's filepath if not exist. Returns file's fp
## ===========================================

mk_dir <- function(fp = "") {
  if (dir.exists(fp)) {
    fp
  } else {
    dir.create(dirname(fp), recursive = TRUE)
    fp
  }
}

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
##   Gaussian mixture models: Each group/cluster should have at least 15 points!!
## ===========================================

#Loop through grouping factors and calculate the parameters for GMM and then assign probabilities for a range of inputs
GMM_probs <- function(df, variable, group1, group2, k = 2, samp = 100, tails = 0.1) {
  
  df <- df %>% as.data.frame()
  variable <- enquo(variable)
  group1 <- enquo(group1)
  group2 <- enquo(group2)
  
  GMM <- data.frame(group1 = NULL, group2 = NULL, lambda = NULL, mu = NULL, sigma = NULL, comp = NULL)
  
  for(t in unique(df %>% select(!!group1) %>% pull())) {
    for(f in unique(df %>% select(!!group2) %>% pull())) {
      df_temp <- df %>%
        filter(grepl(t, !!group1) & grepl(f, !!group2))
      
      set.seed(2)
      em <- df_temp %>%
        select(!!variable) %>%
        pull() %>%
        normalmixEM(., k = k, verb = FALSE)
      
      GMM <- rbind(GMM, data.frame(group1 = t, group2 = f, lambda = em$lambda, mu = em$mu, sigma = em$sigma, comp = 1:k))
    }
  }
  
  GMM %<>%
    rename(!!group1 := "group1", !!group2 := "group2")
  
  prob_df <- df %>%
    group_by(!!group1, !!group2) %>%
    summarize(max_range = max(!!variable), min_range = min(!!variable)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(ranges = list(seq(from = min_range*tails,
                             to = max_range*(1 + tails),
                             length.out = samp))) %>%
    unnest(ranges) %>%
    merge(GMM, by = c(quo_name(group1), quo_name(group2)), all.y = TRUE) %>%
    mutate(probability = lambda*dnorm(ranges, mean = mu, sigma)) %>%
    group_by(!!group1, !!group2) %>%
    mutate(prob_density = n()*probability) %>%
    ungroup() %>%
    rename(!!variable := "ranges")
  
}

## =====
##   BIC
## =====

#BIC is Bayesian information criterion which is one way to estimate how good the GMM is in predicting the observed data.
#Lowest BIC or derivative of BIC scores that doesnt change is best
gmm_bic <- function (m) {
  num_pts <- length(m$x)
  num_params <- length(m$sigma) + length(m$mu)
  (log(num_pts) * num_params) - (2 * m$loglik)
}

## =========
##   get BIC
## =========

mixEM_stats_bic <- function(df, variable, group1, group2, k = 2, max_k = 5) {
  
  df <- df %>% as.data.frame()
  variable <- enquo(variable)
  group1 <- enquo(group1)
  group2 <- enquo(group2)
  
  return_df <- data.frame(group1 = NULL, group2 = NULL, k = NULL, bic = NULL)
  
  for(k_i in k:max_k) {
    for(t in unique(df %>% select(!!group1) %>% pull())) {
      for(f in unique(df %>% select(!!group2) %>% pull())) {
        df_temp <- df %>%
          filter(grepl(t, !!group1) & grepl(f, !!group2))
        
        set.seed(2)
        bic <- df_temp %>%
          select(!!variable) %>%
          pull() %>%
          normalmixEM(., k = k_i, verb = FALSE) %>%
          gmm_bic()
        
        return_df <- rbind(return_df, data.frame(group1 = t, group2 = f, k = k_i, bic = bic))
      }
    }
  }
  return_df %>%
    #group_by(group1, group2) %>%
    #arrange(k) %>%
    #mutate(dBIC = c(diff(bic), 0)) %>%
    #ungroup() %>%
    pivot_wider(id_cols = c("group1", "group2"), names_from = "k", names_prefix = "k ", values_from = "bic") %>%
    rename(!!group1 := "group1", !!group2 := "group2")
  
}

## =============================================
##   Remove underscore and capitalize first char
## =============================================

cap_first <- function(x) {
  stringr::str_to_sentence(gsub("_", " ", x))
}

## ========================
##   Make taxon names nicer
## ========================

make_taxon_label <- function (taxa, other = "Other", show_phylum = FALSE, phylum_sep = " - ") {
  reformat_phylum <- if (show_phylum) {
    function (x) str_replace(x, " ", phylum_sep)
  } else {
    function (x) str_remove(x, "^\\w+ ")
  }
  
  taxa %>%
    if_else(str_detect(., "s__") | (. %in% other), ., paste(., "spp.")) %>%
    reformat_phylum() %>%
    str_remove_all("[kpcofgs]__")
}

## ========================
##   End line with new line
## ========================

start_w_newline <- function(string) {
  gsub("^", "\n", string)
}
