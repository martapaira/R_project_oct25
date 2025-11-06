#' @title Task 1.1: Filter and Summarize Bulk Counts
#'
#' @description Filters bulk RNA-Seq counts for specific conditions and gene patterns,
#' then computes the mean and median expression for the remaining genes.
#'
#' @param dt_counts A data.table containing the columns 'gene', 'count', and 'sample_id'.
#' @param dt_metadata A data.table containing the columns 'sample_id' and 'condition'.
#' @param condition_filter String specifying the condition to keep (e.g., "treated").
#' @param gene_pattern Regex string to filter genes (e.g., "^GENE_00").
#'
#' @return A data.table with columns 'gene', 'media_counts', and 'median_counts'.
#' @import data.table
task1_filter_summarize <- function(dt_counts, dt_metadata, condition_filter, gene_pattern) {
  
  # Pre-join data to get 'condition' column
  dt_joined <- dt_counts[dt_metadata, on = "sample_id"]
  
  # 1.1.2: Data.table Implementation
  dt_result <- dt_joined[
    grepl(gene_pattern, gene) & condition == condition_filter,
    list(media_counts = mean(count), 
         median_counts = median(count)),
    by = gene
  ]
  return(dt_result)
}

#' @title Task 1.2: Join Metadata and Summarize Mean Counts by Condition
#'
#' @description Performs an equi-join between count data and metadata, then calculates
#' the mean expression (media_counts) for every gene within each experimental condition.
#'
#' @param dt_counts A data.table containing the columns 'gene', 'count', and 'sample_id'.
#' @param dt_metadata A data.table containing the columns 'sample_id' and 'condition'.
#'
#' @return A data.table with columns 'gene', 'condition', and 'media_counts'.
#' @import data.table
task1_join_and_summarize <- function(dt_counts, dt_metadata) {
  
  # 1.2.2: Data.table Implementation
  dt_result <- dt_counts[dt_metadata, on = "sample_id"][,
                                                        list(media_counts = mean(count)),
                                                        by = list(gene, condition)
  ]
  return(dt_result)
}

#' @title Task 2.1: Add Log2 Counts and Fixed Threshold Flag (In-Place)
#'
#' @description Adds a log2(count + 1) column and a binary flag 'high' based on a fixed
#' count threshold (> 100). The operation is performed in-place using `data.table:::=`.
#'
#' @param dt_counts_t2 A data.table (modified in-place) containing the 'count' column.
#' @param fixed_threshold Numeric value for the fixed threshold (default: 100).
#'
#' @return The modified data.table, invisibly. Returns NULL to emphasize the in-place operation.
#' @import data.table
task2_add_fixed_qc_columns <- function(dt_counts_t2, fixed_threshold = 100) {
  
  # Check if 'count' exists before proceeding
  if (!"count" %in% names(dt_counts_t2)) {
    stop("Input data.table must contain a column named 'count'.")
  }
  
  # Add log2_counts column in-place
  dt_counts_t2[, log2_counts := log2(count + 1)]
  
  # Add initial 'high' flag based on fixed threshold (Task 2.1)
  dt_counts_t2[, high := ifelse(count > fixed_threshold, 1, 0)]
  
  # In-place operations traditionally return the object invisibly or NULL
  invisible(NULL)
}

#' @title Task 2.2: Overwrite Flag using Gene-wise Median Threshold (In-Place)
#'
#' @description Overwrites the 'high' flag in-place using a group-wise threshold:
#' a count is 'high' (1) if it exceeds the median count for its specific gene.
#'
#' @param dt_counts_t2 A data.table (modified in-place) containing 'count', 'gene', and the 'high' flag.
#'
#' @return The modified data.table, invisibly. Returns NULL to emphasize the in-place operation.
#' @import data.table
task2_overwrite_gene_threshold <- function(dt_counts_t2) {
  
  # Check if necessary columns exist
  if (!all(c("count", "gene") %in% names(dt_counts_t2))) {
    stop("Input data.table must contain columns 'count' and 'gene'.")
  }
  
  # Overwrite 'high' flag using grouped calculation (Task 2.2)
  dt_counts_t2[, high := ifelse(count > median(count, na.rm = TRUE), 1, 0), by = gene]
  
  invisible(NULL)
}

#' @title Task 3.1: Optimized Equi-Join using Setkey
#'
#' @description Performs an equi-join between count data and metadata. It optimizes the
#' join by setting the key on the metadata table first, enabling a fast binary search.
#'
#' @param dt_counts The larger data.table (counts, the 'i' table).
#' @param dt_metadata The smaller data.table (metadata, the 'x' table). This table will be keyed in-place.
#'
#' @return A data.table containing the joined result (counts + metadata).
#' @import data.table
task3_optimized_join <- function(dt_counts, dt_metadata) {
  
  # 1. Define the key on the lookup table (metadata)
  setkey(dt_metadata, sample_id)
  
  # 2. Execute the join: dt_metadata[dt_counts] on="sample_id"
  dt_result <- dt_metadata[dt_counts, on = "sample_id", nomatch=0]
  
  return(dt_result)
}

#' @title Task 3.2: Add Secondary Index for Fast Lookups (In-Place)
#'
#' @description Adds a secondary index on the specified columns of the data.table
#' to accelerate future subsetting and filtering operations. This operation is performed in-place.
#'
#' @param dt_table The data.table to modify in-place.
#' @param index_cols Character vector of column names to index (e.g., c("gene", "sample_id")).
#'
#' @return The modified data.table, invisibly. Returns NULL to emphasize the in-place operation.
#' @import data.table
task3_add_secondary_index <- function(dt_table, index_cols = c("gene", "sample_id")) {
  
  # Add secondary index
  setindexv(dt_table, index_cols)
  
  invisible(NULL)
}

#' @title Task 4.1: Compute Total Counts per Patient
#'
#' @description Joins counts with metadata to link samples to patients, then computes
#' the total aggregated count (sum(count)) for each patient across all genes/samples.
#'
#' @param dt_counts A data.table containing 'count' and 'sample_id'.
#' @param dt_metadata A data.table containing 'sample_id' and 'patient_id'.
#'
#' @return A data.table with 'patient_id' and 'total_counts'.
#' @import data.table
task4_total_counts_per_patient <- function(dt_counts, dt_metadata) {
  
  # 4.1.2: Data.table Implementation
  dt_result <- dt_counts[dt_metadata, on = "sample_id"][,
                                                        list(total_counts = sum(count)),
                                                        by = patient_id 
  ]
  return(dt_result)
}

#' @title Task 4.2: Find Top N Genes by Mean Count per Condition
#'
#' @description Performs a "Top N per Group" analysis. It calculates the mean count 
#' for every gene within each condition and returns the top N genes with the highest 
#' mean expression per condition (default N=10).
#'
#' @param dt_counts A data.table containing 'gene', 'count', and 'sample_id'.
#' @param dt_metadata A data.table containing 'sample_id' and 'condition'.
#' @param N Numeric value specifying the number of top genes to return per condition (default 10).
#'
#' @return A data.table containing the top N genes per condition, including 'gene', 'condition', 'mean_count', and 'rank'.
#' @import data.table
task4_top_n_genes_by_condition <- function(dt_counts, dt_metadata, N = 10) {
  
  # 1. Calculate mean per (gene, condition)
  dt_ranked <- dt_counts[dt_metadata, on = "sample_id"][,
                                                        list(mean_count = mean(count)),
                                                        by = list(gene, condition)
  ]
  
  # 2. Assign rank and filter
  dt_result <- dt_ranked[
    # Assign a rank for each condition group
    , rank := rank(-mean_count, ties.method = "min"), 
    by = condition
  ][
    # Filter only those with rank <= N
    rank <= N
  ]
  
  return(dt_result)
}

#' @title Task 5.1: Classify Lab Values using Non-Equi Join (Normal vs. Out-of-Range)
#'
#' @description Classifies clinical lab values (dt_labs) as "Normal" or "Out_of_Range"
#' by performing a non-equi join against reference intervals (dt_ranges).
#'
#' @param dt_labs A data.table containing 'patient_id', 'lab_name', and 'value'.
#' @param dt_ranges A data.table containing 'lab_name', 'lower', and 'upper' reference bounds.
#'
#' @return A data.table containing all lab measurements with an added 'status' column 
#' (either 'Normal' or 'Out_of_Range').
#' @import data.table
task5_classify_lab_status <- function(dt_labs, dt_ranges) {
  
  # This implementation uses the Join + Update method (X[I, j:=V]) 
  # to assign 'Normal' status only where intervals overlap.
  
  # 1. Preparation: Copy and standardize column names
  dt_labs <- copy(dt_labs)
  dt_ranges <- copy(dt_ranges)
  setnames(dt_labs, old = "lab", new = "lab_name", skip_absent = TRUE)
  setnames(dt_ranges, old = "lab", new = "lab_name", skip_absent = TRUE)
  
  # 2. Assign default status to 'Out_of_Range'
  dt_labs[, status := "Out_of_Range"]
  
  # 3. Define keys on the range limits for the join
  setkey(dt_labs, lab_name)
  setkey(dt_ranges, lab_name)
  
  # 4. Execute the Non-Equi Join and UPDATE status to 'Normal' for matching rows.
  # dt_labs_local[dt_ranges_local, on=...]: Joins dt_labs_local using dt_ranges_local criteria.
  # j=status:="Normal": Updates the 'status' column in dt_labs_local IN-PLACE for corresponding rows.
  dt_labs[dt_ranges, 
          on = .(lab_name, value >= lower, value <= upper),
          j = status := "Normal"]
  
  # 5. Clean up temporary columns that might be created by the join
  dt_labs[, c("lower", "upper") := NULL] 
  
  return(dt_labs)
}

#' @title Task 5.2: Summarize Abnormal Rates by Patient and Lab
#'
#' @description Calculates the total number of tests and the rate of abnormal results
#' (AbnormalRate) for each combination of patient and lab test.
#'
#' @param dt_classified_data A data.table resulting from Task 5.1 (dt_final_t5_base) 
#' containing 'patient_id', 'lab_name', and 'status'.
#'
#' @return A data.table with 'patient_id', 'lab_name', 'TotalTests', 'AbnormalTests', and 'AbnormalRate'.
#' @import data.table
task5_summarize_abnormal_rates <- function(dt_classified_data) {
  
  # 5.2.2: Data.table Implementation
  dt_result <- dt_classified_data[, 
                                  list(
                                    TotalTests = .N,
                                    AbnormalTests = sum(status == "Out_of_Range"),
                                    AbnormalRate = sum(status == "Out_of_Range") / .N
                                  ),
                                  by = list(patient_id, lab_name)
  ]
  
  return(dt_result)
}

#' @title Task 6.1: Optimized Rolling Join for Nearest Vitals
#'
#' @description Performs time-series preparation and an optimized rolling join (`roll=TRUE`)
#' to attach the nearest prior HR and SBP vital readings to each lab measurement.
#'
#' @param dt_labs A data.table of lab measurements with 'patient_id' and 'time_iso'.
#' @param dt_vitals A data.table of vital sign time-series data with 'patient_id', 'time_iso', 'vital', and 'value'.
#'
#' @return A data.table containing merged lab and vital data, including the calculated 'lag_minutes'.
#' @import data.table
#' @importFrom lubridate ymd_hms
task6_rolling_join_vitals <- function(dt_labs, dt_vitals) {
  
  # TEMPORAL DATA PREPARATION
  # 1. Convert text to POSIXct datetime object (Crucial for time calculation)
  if ("time_iso" %in% names(dt_labs)) { dt_labs[, time_iso := lubridate::ymd_hms(as.character(time_iso), tz = "UTC")] }
  if ("time_iso" %in% names(dt_vitals)) { dt_vitals[, time_iso := lubridate::ymd_hms(as.character(time_iso), tz = "UTC")] }
  
  # 2. PIVOT: Transform vitals from 'long' to 'wide' (HR and SBP as columns)
  dt_vitals_wide <- dcast(dt_vitals, patient_id + time_iso ~ vital, value.var = "value", fun.aggregate = mean, na.rm = TRUE)
  
  # Order the lookup data (vitals) for the rolling join 
  setorder(dt_vitals_wide, patient_id, time_iso)
  
  # Rename time column in 'x' (vitals) BEFORE rolling join
  setnames(dt_vitals_wide, "time_iso", "vitals_time", skip_absent = TRUE)
  
  # Preserve 'vitals_time' in a copy, because the join will overwrite the key column.
  dt_vitals_wide[, vitals_time_matched := vitals_time]
  
  # --- EXECUTE ROLLING JOIN ---
  # 3. Execute Join: dt_vitals_wide[dt_labs]
  # The join matches on patient_id (equi-join) and time (rolling join)
  dt_join_result <- dt_vitals_wide[dt_labs, 
                                   on = .(patient_id, vitals_time = time_iso), 
                                   roll = TRUE, 
                                   mult = "last"]
  
  # --- POST-PROCESSING ---
  
  # 4. Rename columns after the join
  setnames(dt_join_result, 
           old = c("vitals_time", "vitals_time_matched", "HR", "SBP"), 
           new = c("lab_time", "vitals_time", "nearest_hr", "nearest_sbp"), 
           skip_absent = TRUE)
  
  # 5. Calculate lag (lab_time - vitals_time)
  dt_join_result[, lag_minutes := as.numeric(difftime(lab_time, vitals_time, units = "mins"))]
  
  # 6. Rename the lab value column
  setnames(dt_join_result, old = "value", new = "lab_value", skip_absent = TRUE)
  setnames(dt_join_result, old = "lab", new = "lab_name", skip_absent = TRUE) # Assicurati che sia lab_name
  
  # Clean up temporary columns (if any are created during the join)
  cols_to_remove <- intersect(names(dt_join_result), c("i.patient_id", "i.time_iso", "i.vital", "i.value"))
  if (length(cols_to_remove) > 0) dt_join_result[, (cols_to_remove) := NULL]
  
  return(dt_join_result)
}

#' @title Task 6.2: Summarize Correlation (CRP vs. Vitals) per Patient
#'
#' @description Filters the joined dataset for CRP measurements, ensures sufficient data 
#' points (n >= 2), and computes the Pearson correlation between CRP and the nearest 
#' measured vital signs (HR and SBP) for each patient.
#'
#' @param dt_joined_data A data.table resulting from Task 6.1, containing 'patient_id', 
#' 'lab_name', 'lab_value', 'nearest_hr', and 'nearest_sbp'.
#'
#' @return A data.table with 'patient_id', 'cor_CRP_HR', and 'cor_CRP_SBP'.
#' @import data.table
task6_summarize_correlations <- function(dt_joined_data) {
  
  # 1. Filter for CRP and complete rows (no NA in vitals or lab_value)
  dt_cor_base <- dt_joined_data[lab_name == "CRP" & !is.na(nearest_hr) & !is.na(nearest_sbp) & !is.na(lab_value)]
  
  # 2. Identify patients with sufficient data points (n >= 2) to calculate correlation
  dt_valid_patients <- dt_cor_base[, .N, by = patient_id][N >= 2]
  
  if (nrow(dt_valid_patients) == 0) {
    warning("Nessun paziente trovato con almeno 2 misurazioni CRP/Vitals complete per la correlazione.")
    return(data.table(patient_id=character(), cor_CRP_HR=numeric(), cor_CRP_SBP=numeric()))
  }
  
  # 3. Join back to the base data and calculate the grouped correlation
  dt_result <- dt_cor_base[dt_valid_patients, on = "patient_id"][, 
                                                                 list(
                                                                   cor_CRP_HR = cor(lab_value, nearest_hr, use = "pairwise.complete.obs"),
                                                                   cor_CRP_SBP = cor(lab_value, nearest_sbp, use = "pairwise.complete.obs")
                                                                 ),
                                                                 by = patient_id
  ]
  
  return(dt_result)
}

#' @title Task 7.1: Filter ATAC Peaks by Genomic Region
#'
#' @description Filters ATAC peaks data based on chromosome and a start position range, 
#' returning only the peaks that fall within the specified genomic window.
#'
#' @param dt_peaks A data.table of ATAC peaks, containing 'chr', 'start', 'end', and 'score'.
#' @param chr_target String specifying the target chromosome (e.g., "chr2").
#' @param start_min Numeric minimum boundary for the start position (e.g., 2000000).
#' @param start_max Numeric maximum boundary for the start position (e.g., 4000000).
#'
#' @return A data.table containing the subset of peaks within the window.
#' @import data.table
task7_filter_genomic_region <- function(dt_peaks, chr_target, start_min, start_max) {
  
  # 7.1.2: Data.table Implementation (Filter)
  dt_result <- dt_peaks[chr == chr_target & start >= start_min & start <= start_max]
  
  return(dt_result)
}

#' @title Task 7.2: Find Top N Peaks by Score
#'
#' @description Takes a filtered data.table of peaks, orders them by the 'score' column 
#' in descending order, and returns the top N peaks.
#'
#' @param dt_filtered_peaks A data.table of peaks, assumed to contain the 'score' column.
#' @param N Numeric value specifying the number of top peaks to return (default 50).
#'
#' @return A data.table containing the top N peaks by score.
#' @import data.table
task7_find_top_n_by_score <- function(dt_filtered_peaks, N = 50) {
  
  # 7.2.2: Data.table Implementation (Sorting and Head)
  
  # 1. Copy date to not modify external table with setorder
  dt_peaks_copy <- copy(dt_filtered_peaks)
  
  # 2. Order in-place to score  (-score)
  setorder(dt_peaks_copy, -score)
  
  # 3. Take first N
  dt_result <- head(dt_peaks_copy, N)
  
  return(dt_result)
}

#' @title Task 8.1: Calculate Robust Summary Statistics
#'
#' @description Joins counts with metadata and calculates robust summary statistics 
#' (mean, median, Q1/Q3) for each gene within each experimental condition.
#'
#' @param dt_counts A data.table containing 'count' and 'sample_id'.
#' @param dt_metadata A data.table containing 'sample_id' and 'condition'.
#'
#' @return A data.table with 'gene', 'condition', and robust summary metrics.
#' @import data.table
task8_calculate_robust_stats <- function(dt_counts, dt_metadata) {
  
  # 1. Join data
  dt_joined <- dt_counts[dt_metadata, on = "sample_id"]
  
  # 2. Calculate statistics by group (gene, condition)
  dt_result <- dt_joined[,
                         list(
                           mean_count = mean(count, na.rm = TRUE),
                           median_count = median(count, na.rm = TRUE),
                           Q1 = quantile(count, 0.25, na.rm = TRUE),
                           Q3 = quantile(count, 0.75, na.rm = TRUE)
                         ),
                         by = list(gene, condition)
  ]
  
  return(dt_result)
}

#' @title Task 8.2: Filter Genes based on Mean Difference (Treated > Control)
#'
#' @description Takes the robust statistics table (from T8.1), pivots it to a wide 
#' format using `dcast`, and filters for genes where the mean count in the 'treated' 
#' condition is strictly greater than the mean count in the 'control' condition.
#'
#' @param dt_stats_long A data.table in long format (from T8.1) containing metrics for 
#' each (gene, condition) pair.
#'
#' @return A data.table in wide format containing only the filtered genes.
#' @import data.table
task8_filter_mean_difference <- function(dt_stats_long) {
  
  # 1. Pivot the table to wide format (one row per gene)
  dt_wide <- dcast(dt_stats_long, 
                   gene ~ condition, 
                   value.var = c("mean_count", "median_count", "Q1", "Q3"))
  
  # 2. Filter the result
  dt_result <- dt_wide[mean_count_treated > mean_count_control]
  
  return(dt_result)
}

#' @title Task 9: Wide-to-Long-to-Wide Transformation Pipeline
#'
#' @description Transforms a wide count matrix (gene x sample) into a long format, 
#' calculates per-sample totals, and then re-pivots the data to a wide summary 
#' format (gene x condition mean counts). This demonstrates the efficiency of 
#' `melt` and `dcast`.
#'
#' @param dt_wide_counts A data.table in wide format (gene x samples).
#' @param dt_metadata A data.table containing 'sample_id' and 'condition'.
#'
#' @return A data.table in wide format with mean counts grouped by gene and condition.
#' @import data.table
task9_wide_long_wide_transform <- function(dt_wide_counts, dt_metadata) {
  
  # Identify index of sample columns
  sample_col_indices <- 2:ncol(dt_wide_counts)
  
  # 1. Wide -> Long (melt)
  dt_long <- melt(dt_wide_counts,
                  id.vars = "gene",
                  measure.vars = sample_col_indices,
                  variable.name = "sample_id",
                  value.name = "count")
  
  # 2. Add per-sample totals 
  dt_sample_totals <- dt_long[, .(sample_total = sum(count, na.rm = TRUE)), by = sample_id]
  dt_long_with_totals <- dt_sample_totals[dt_long, on = "sample_id"]
  
  # 3. Aggregate (gene, condition) and Long -> Wide (dcast)
  
  # Join with metadata 
  dt_joined <- dt_metadata[dt_long_with_totals, on = "sample_id"]
  
  # Pivot and aggregation
  dt_result <- dcast(dt_joined, 
                     gene ~ condition, 
                     value.var = "count", 
                     fun.aggregate = mean)
  
  return(dt_result)
}

#' @title Helper Function: Calculate Interval Overlap Length
#'
#' @description Calculates the non-negative length of the overlap between two intervals.
#' This function is essential for Task 10 and 11 when converting genomic coordinates 
#' (start1, end1) and (start2, end2) into a physical length.
#'
#' @param start1 Start position of the first interval.
#' @param end1 End position of the first interval.
#' @param start2 Start position of the second interval.
#' @param end2 End position of the second interval.
#'
#' @return Numeric vector representing the overlap length (or 0 if no overlap).
calc_overlap_length <- function(start1, end1, start2, end2) {
  pmax(0, pmin(end1, end2) - pmax(start1, start2))
}

#' @title Task 10.1: Map ATAC Peaks to Genes via Spatial Join (foverlaps)
#'
#' @description Performs an optimized genomic interval join using `data.table::foverlaps` 
#' to map ATAC peaks to overlapping gene regions. It counts the number of overlapping 
#' peaks per gene and calculates the overlap length per pair.
#'
#' @param dt_peaks A data.table containing ATAC peaks with 'chr', 'start', 'end'.
#' @param dt_genes A data.table containing gene annotations with 'chr', 'start', 'end', and 'gene' ID.
#'
#' @return A data.table containing all overlapping pairs and the calculated overlap length ('overlap_len').
#' @import data.table
task10_calculate_overlap_stats <- function(dt_peaks, dt_genes) {
  
  # 1. Preparation and Keys (necessary for foverlaps)
  # Copies created here are local to the function call
  dt_peaks_local <- copy(dt_peaks)
  dt_genes_local <- copy(dt_genes)
  
  # Rename the gene file columns to avoid conflicts and prepare keys
  setnames(dt_genes_local, old = c("start", "end"), new = c("gene_start", "gene_end"))
  
  # Set keys for the spatial join (mandatory for foverlaps)
  setkey(dt_peaks_local, chr, start, end)
  setkey(dt_genes_local, chr, gene_start, gene_end)
  
  # 2. Spatial join (foverlaps)
  dt_overlapped_final <- foverlaps(
    dt_peaks_local,
    dt_genes_local,
    by.x = c("chr", "start", "end"),
    by.y = c("chr", "gene_start", "gene_end"),
    nomatch = 0,
    type = "any"
  )
  
  # 3. Calculate the overlap length (reusing the helper function)
  dt_overlapped_final[, overlap_len := calc_overlap_length(start, end, gene_start, gene_end)]
  
  # 4. Select essential columns for the output (prepares the data for aggregation)
  dt_result <- dt_overlapped_final[, .(
    gene, 
    chr = i.chr, 
    peak_start = start, 
    peak_end = end, 
    overlap_len
  )]
  
  return(dt_result)
}

#' @title Task 10.2: Summarize Overlap Length and Find Top 20 Genes
#'
#' @description Aggregates the results from the spatial join to calculate the total 
#' base-pair overlap length per gene and returns the top 20 genes with the longest 
#' total overlap.
#'
#' @param dt_overlapped_stats A data.table containing 'gene' and 'overlap_len' (from T10.1).
#'
#' @return A data.table with 'gene' and 'total_overlap_bp' for the top 20 genes.
#' @import data.table
task10_summarize_and_find_top_n <- function(dt_overlapped_stats) {
  
  # 1. Aggregate gene and sum overlap lenght
  dt_overlap_sum <- dt_overlapped_stats[
    , .(total_overlap_bp = sum(overlap_len)),
    by = gene
  ]
  
  # 2. Order and select Top 20 genes
  dt_result <- dt_overlap_sum[order(-total_overlap_bp)][1:20] 
  
  return(dt_result)
}

#' @title Task 11.1/11.2: Map SNPs to Genes and Summarize HIGH Impact Counts
#'
#' @description Converts single SNP positions to 1-bp intervals, performs an optimized 
#' Non-Equi Join to map HIGH impact variants to overlapping genes, and summarizes 
#' the count of HIGH impact variants per gene and per sample.
#'
#' @param dt_variants A data.table containing 'chr', 'pos', 'impact', and 'sample_id'.
#' @param dt_genes A data.table containing gene annotations ('chr', 'start', 'end', 'gene').
#'
#' @return A data.table summarized by 'gene' and 'sample_id' with 'high_impact_count'.
#' @import data.table
task11_map_and_summarize_high_impact <- function(dt_variants, dt_genes) {
  
  # 1. Preparation: Local copies and convert SNP position to 1-bp interval
  dt_variants_local <- copy(dt_variants)
  dt_genes_local <- copy(dt_genes)
  
  dt_variants_local[, variant_start := pos]
  dt_variants_local[, variant_end := pos]
  
  # 2. Filter HIGH impact variants BEFORE the join for efficiency
  dt_variants_high <- dt_variants_local[impact == 'HIGH']
  
  # 3. Set keys for spatial join (genes has standard start/end names)
  setkey(dt_variants_high, chr, variant_start, variant_end)
  setkey(dt_genes_local, chr, start, end)
  
  # 4. Execute the Non-Equi Join (Spatial): GENES[HIGH_VARIANTS, on=...]
  dt_overlapped_snps <- dt_genes_local[dt_variants_high, 
                                       on = .(chr, start <= variant_end, end >= variant_start),
                                       nomatch = 0,
                                       allow.cartesian = TRUE
  ]
  
  # 5. Summarize (Task 11.2): Count variants per gene and sample
  dt_result <- dt_overlapped_snps[, 
                                  .(high_impact_count = .N), 
                                  by = .(gene, sample_id)
  ]
  
  return(dt_result)
}

#' @title Task 11.3: Find Genes Hit by HIGH Variants Across All Samples
#'
#' @description Identifies which genes in the summary table (from T11.2) have at least 
#' one HIGH impact variant in every single sample present in the input.
#'
#' @param dt_high_impact_summary A data.table summarized by 'gene' and 'sample_id' 
#' (from T11.1/11.2).
#' @param dt_all_variants The original dt_variants table, used only to determine 
#' the total number of unique samples (total universe).
#'
#' @return A data.table containing only the 'gene' IDs that were hit in all samples.
#' @import data.table
task11_find_genes_in_all_samples <- function(dt_high_impact_summary, dt_all_variants) {
  
  # 1. Determine the total number of unique samples (the universe)
  total_unique_samples <- uniqueN(dt_all_variants$sample_id)
  
  # 2. Count how many unique samples hit each gene
  dt_hit_counts <- dt_high_impact_summary[, 
                                          .(n_samples_hit = uniqueN(sample_id)), 
                                          by = gene
  ]
  
  # 3. Filter only genes where the hit count equals the total unique sample count
  dt_result <- dt_hit_counts[n_samples_hit == total_unique_samples][, .(gene)]
  
  return(dt_result)
}

#' @title Task 12.1: Combine Cohorts Safely and Order
#'
#' @description Combines two potentially heterogeneous cohort data.tables using 
#' `rbindlist` with safe options (`use.names=TRUE, fill=TRUE`) and orders the result.
#'
#' @param dt_cohortA A data.table containing sample metadata for Cohort A.
#' @param dt_cohortB A data.table containing sample metadata for Cohort B.
#'
#' @return A combined data.table, ordered by 'cohort', 'condition', and 'sample_id'.
#' @import data.table
task12_combine_cohorts <- function(dt_cohortA, dt_cohortB) {
  
  # 1. Combine using rbindlist for safe column alignment and fill
  dt_combined <- rbindlist(list(dt_cohortA, dt_cohortB), 
                           use.names = TRUE, 
                           fill = TRUE)
  
  # 2. Order the result in-place
  setorder(dt_combined, cohort, condition, sample_id)
  
  return(dt_combined)
}

#' @title Task 12.2: Calculate Mean Counts for Top 100 Most Variable Genes
#'
#' @description Calculates the Top 100 most variable genes across the dataset, 
#' filters the bulk counts for these genes, joins with the combined cohort metadata,
#' and computes the mean count grouped by cohort, condition, and gene.
#'
#' @param dt_counts A data.table containing bulk counts ('gene', 'count', 'sample_id').
#' @param dt_combined_cohorts A data.table containing combined cohort metadata (from T12.1).
#'
#' @return A data.table with mean counts for the Top 100 genes, grouped by cohort, condition, and gene.
#' @import data.table
task12_calculate_top_n_means <- function(dt_counts, dt_combined_cohorts) {
  
  # 1. Find Top 100 genes most variable (using variance)
  top100_genes <- dt_counts[, .(variance = var(count, na.rm = TRUE)), by = gene][
    order(-variance)
  ][1:100, gene]
  
  # 2. Join e Calculate mean 
  dt_result <- dt_counts[gene %in% top100_genes][ 
    dt_combined_cohorts, on = "sample_id", nomatch = 0 
  ][,
    .(mean_count = mean(count, na.rm = TRUE)), 
    by = .(cohort, condition, gene) 
  ]
  
  return(dt_result)
}


#' @title Final Revision: Combine and Clean Single-Cell Metadata
#'
#' @description Performs necessary cleaning on the 'cell' IDs (trimming, removing suffixes like _N/_T)
#' and executes an inner join between the clusters data and the cell types data to create a master table.
#'
#' @param dt_clusters A data.table containing 'cell' ID and 'integration_cluster' (from SeuratIntegration file).
#' @param dt_celltypes A data.table containing 'cell' ID, 'cell_type', and 'sample_type' (N/T).
#'
#' @return A combined data.table ('master table') containing 'cell', 'integration_cluster', 'cell_type', and 'sample_type'.
#' @import data.table
finaltask_combine_and_clean <- function(dt_clusters, dt_celltypes) {
  
  # 1. Preparation: Remove whitespace and problematic suffix (_X_)
  dt_clusters_local <- copy(dt_clusters)
  dt_celltypes_local <- copy(dt_celltypes)
  
  # Remove _T_, _N_, etc., to ensure cell IDs match.
  dt_clusters_local[, cell := trimws(as.character(cell))]
  dt_celltypes_local[, cell := trimws(as.character(cell))]
  dt_clusters_local[, cell := gsub("_[A-Z]_", "", cell)] 
  
  # 2. Join (inner join: nomatch = 0)
  setkey(dt_clusters_local, cell)
  setkey(dt_celltypes_local, cell)
  
  dt_master_table <- dt_clusters_local[dt_celltypes_local, nomatch = 0]
  
  return(dt_master_table)
}

#' @title Final Revision: Count Cell Types per Cluster (FR-Task 2)
#'
#' @description Counts the total number of cells for each combination of 'integration_cluster' 
#' and 'cell_type' in the master table.
#'
#' @param dt_master_table The combined single-cell data.table (from FR-T1).
#'
#' @return A data.table with columns 'integration_cluster', 'cell_type', and 'count'.
#' @import data.table
finaltask_count_per_cluster <- function(dt_master_table) {
  
  # Counts the grouped rows by cluster and cell type
  dt_counts <- dt_master_table[, .N, by = .(integration_cluster, cell_type)]
  setnames(dt_counts, "N", "count")
  
  return(dt_counts)
}

#' @title Final Revision: Normalize Cell Counts by Cluster and Tissue Type (FR-Task 3 & 5)
#'
#' @description Calculates the cell count per (Cluster, Cell Type, Tissue Type) and 
#' computes the normalized percentage of each cell type within its respective 
#' (Cluster x Tissue) group.
#'
#' @param dt_master_table The combined single-cell data.table (from FR-T1).
#'
#' @return A data.table with counts, total_in_group, and 'percentage' (normalized).
#' @import data.table
finaltask_normalize_by_tissue <- function(dt_master_table) {
  
  # 1. Counting: Group by the three factors
  dt_counts_tissue <- dt_master_table[, .N, by = .(integration_cluster, cell_type, sample_type)]
  setnames(dt_counts_tissue, "N", "count")
  
  # 2. Normalization: Calculate total per group and percentage (in-place)
  dt_counts_tissue[, total_in_group := sum(count), by = .(integration_cluster, sample_type)]
  dt_counts_tissue[, percentage := (count / total_in_group) * 100]
  
  return(dt_counts_tissue)
}