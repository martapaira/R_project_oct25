# Project: Bioinformatics Performance Analysis: data.table vs. data.frame

# This script implements the core tasks, comparing the performance of three approaches:
# 1. data.frame (dplyr)
# 2. data.table
# 3. sqldf (SQL)

# The script is run within the Docker container to generate reproducible results.

# ----------------------------------------------------
# 1. INITIAL CONFIGURATION AND LIBRARIES
# ----------------------------------------------------
library(data.table)
library(dplyr)      # Used for the data.frame approach
library(sqldf)      # Third package: SQL syntax
library(tidyr)
library(ggplot2)
library(scales)

# --- LOAD DEVELOPED R PACKAGE ---
# Loads the functions created in R/task_functions.R
library(ProgettoBioinfo)

# --- OUTPUT CONFIGURATION ---
output_dir <- "results"
# Creates the 'results' folder if it doesn't exist. showWarnings = FALSE prevents needless messages.
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE)
}

# Initializes the table for performance comparison
performance_comparison <- data.table(Task = character(), 
                                     `Time_DataFrame (s)` = numeric(), 
                                     `Time_DataTable (s)` = numeric(),
                                     `Time_SQL (s)` = numeric())

# Robust function to extract time from system.time()
extract_time <- function(time_obj) {
  # Extracts the "elapsed" element (total time elapsed, rounded to six decimal places)
  return(round(time_obj["elapsed"], 6))
}

# Data Loading
# Files are located directly in /app
dt_counts <- fread("bulk_counts_long.csv")
dt_metadata <- fread("sample_metadata.csv")

# Conversion to data.frame for dplyr and sqldf benchmark
df_counts <- as.data.frame(dt_counts) 
df_metadata <- as.data.frame(dt_metadata)


# ----------------------------------------------------
# TASK 1: Filter and Summarize Bulk Counts
# Goal: Filter, summarize, and group bulk counts.
# --------------------------------------------------
cat("\n--- START TASK 1: PERFORMANCE COMPARISON ---\n")

# --- DATA PREPARATION FOR TASK 1.1 and 1.2 ---
# Preliminary join to get the 'condition' column necessary for filtering
df_unita <- df_counts %>% left_join(df_metadata, by = "sample_id")
dt_unita <- dt_counts[dt_metadata, on = "sample_id"] 

# ----------------------------------------------------
# TASK 1.1: Filter and Summarize (mean/median)
# ----------------------------------------------------
# Keep only treated samples and genes starting with GENE_00, groups by gene and calculate mean and median counts.

cat("\nRunning Task 1.1: Filter & Summarize...\n")

# 1.1.1: Data.frame / Dplyr 
time_df_1_1_result <- system.time({
  df_risultato_1_1 <- df_unita %>% 
    filter(grepl("^GENE_00", gene) & condition == "treated") %>% 
    group_by(gene) %>% 
    summarise(media_counts = mean(count), median_counts = median(count), .groups = 'drop')
})
time_df_1_1 <- extract_time(time_df_1_1_result)
cat("DataFrame T1.1 Time:", time_df_1_1, "seconds.\n")

# 1.1.2: Data.table
time_dt_1_1_result <- system.time({
  dt_risultato_1_1 <- dt_unita[ 
    grepl("^GENE_00", gene) & condition == "treated", 
    list(media_counts = mean(count), median_counts = median(count)),
    by = gene 
  ]
})
time_dt_1_1 <- extract_time(time_dt_1_1_result)
cat("data.table T1.1 Time:", time_dt_1_1, "seconds.\n")

# 1.1.3: SQL (sqldf)
time_sql_1_1_result <- system.time({
  sql_risultato_1_1 <- sqldf::sqldf("
    SELECT T1.gene, AVG(T1.count) AS media_counts, MEDIAN(T1.count) AS median_counts
    FROM df_counts AS T1  
    LEFT JOIN df_metadata AS T2 ON T1.sample_id = T2.sample_id  
    WHERE T2.condition = 'treated' AND T1.gene LIKE 'GENE_00%'
    GROUP BY T1.gene
  ")
})
time_sql_1_1 <- extract_time(time_sql_1_1_result)
cat("SQL T1.1 Time:", time_sql_1_1, "seconds.\n")

# Update comparative table after Task 1.1
performance_comparison <- rbindlist(list(performance_comparison, 
                                         list("T1.1 Filter & Summarize", time_df_1_1, time_dt_1_1, time_sql_1_1)),
                                    use.names = TRUE, fill = TRUE)
print(head(df_risultato_1_1, 5))


# ----------------------------------------------------
# TASK 1.2: Join metadata and summarize by condition
# ----------------------------------------------------
# Group by both gene and condition to calculate the mean count for every gene in both the 'treated' and 'control' groups.

cat("\nRunning Task 1.2: Join & Summarize...\n")

# 1.2.1: Data.frame / Dplyr
time_df_1_2_result <- system.time({
  df_risultato_1_2 <- df_counts %>%
    left_join(df_metadata, by = "sample_id") %>%
    group_by(gene, condition) %>%
    summarise(media_counts = mean(count), .groups = 'drop') #create a nex column "media_counts" with all the calculated means 
})
time_df_1_2 <- extract_time(time_df_1_2_result)
cat("DataFrame T1.2 Time:", time_df_1_2, "seconds.\n")

# 1.2.2: Data.table
time_dt_1_2_result <- system.time({
  dt_risultato_1_2 <- dt_counts[dt_metadata, on = "sample_id"][,
                                                               list(media_counts = mean(count)),
                                                               by = list(gene, condition)
  ]
})
time_dt_1_2 <- extract_time(time_dt_1_2_result)
cat("data.table T1.2 Time:", time_dt_1_2, "seconds.\n")

# 1.2.3: SQL (sqldf)
time_sql_1_2_result <- system.time({
  sql_risultato_1_2 <- sqldf::sqldf("
    SELECT T1.gene, T2.condition, AVG(T1.count) AS media_counts
    FROM df_counts AS T1
    LEFT JOIN df_metadata AS T2 ON T1.sample_id = T2.sample_id
    GROUP BY T1.gene, T2.condition
  ")
})
time_sql_1_2 <- extract_time(time_sql_1_2_result)
cat("SQL T1.2 Time:", time_sql_1_2, "seconds.\n")

# Update comparative table after Task 1.2
performance_comparison <- rbindlist(list(performance_comparison, 
                                         list("T1.2 Join & Summarize", time_df_1_2, time_dt_1_2, time_sql_1_2)),
                                    use.names = TRUE, fill = TRUE)
print(head(df_risultato_1_2, 5))

# ----------------------------------------------------
# FINAL RESULTS SAVING
# ----------------------------------------------------
# Output files are saved into the 'results' subfolder
fwrite(dt_risultato_1_1, file.path(output_dir, "Task1_FilterSummarize_Results.csv"))
fwrite(dt_risultato_1_2, file.path(output_dir,"Task1_JoinSummarize_Results.csv"))

cat("\nTask 1 analysis completed. Results and performance saved.\n")



# ----------------------------------------------------
# TASK 2: Add QC-Style Derived Columns (No Copying)
# Goal: Add log2_counts and overwrite the 'high' flag based on a gene-wise threshold.
# ----------------------------------------------------
cat("\n--- START TASK 2: QC DERIVED COLUMNS ---\n")

# Initialize a local performance table for Task 2
performance_comparison_task2 <- data.table(Task = character(), 
                                           `Time_DataFrame (s)` = numeric(), 
                                           `Time_DataTable (s)` = numeric())
                                           
# Create clean copies of the ORIGINAL data for in-place operations.
dt_counts_t2 <- fread("bulk_counts_long.csv")
df_counts_t2 <- as.data.frame(dt_counts_t2) 

# ----------------------------------------------------
# TASK 2.1: Add log2 and fixed flag (count > 100)
# ----------------------------------------------------
# Add two new columns: a standard log2(count + 1) transformation and a simple binary flag for any count greater than 100.
cat("\nRunning Task 2.1: Log2 & Fixed Flag (3-way comparison)...\n")

# --- 2.1.1: Data.frame / Dplyr ---
time_df_2_1_result <- system.time({
  df_counts_t2_result <- df_counts_t2 %>%
    mutate(log2_counts = log2(count + 1), # Mutate adds two new columns (log2_counts and high). (count+1) to avoid problems if any count=0
           high = ifelse(count > 100, 1, 0))  #If count > 0 is TRUE write 1, otherwise 0
})
time_df_2_1 <- extract_time(time_df_2_1_result)
cat("DataFrame T2.1 Time:", time_df_2_1, "seconds.\n")


# --- 2.1.2: Data.table (Operation In-Place :=) ---
time_dt_2_1_result <- system.time({
  dt_counts_t2[, log2_counts := log2(count + 1)]
  dt_counts_t2[, high := ifelse(count > 100, 1, 0)]
})
time_dt_2_1 <- extract_time(time_dt_2_1_result)
cat("data.table T2.1 Time:", time_dt_2_1, "seconds.\n")


# Update performance table Task 2.1
performance_comparison_task2 <- rbindlist(list(performance_comparison_task2, 
                                               list("T2.1 Log2 & Fixed Flag", time_df_2_1, time_dt_2_1)), use.names = TRUE, fill =TRUE)

print(head(dt_counts_t2, 5))

# ----------------------------------------------------
# TASK 2.2: Overwrite Flag by Gene Threshold (count > median(count) by=gene)
# ----------------------------------------------------
# Overwrites the high column created in T2.1, flagging counts only if they are greater than the median count for that specific gene.

cat("\nRunning Task 2.2: Overwrite Flag by Gene (3-way comparison)...\n")

# --- 2.2.1: Data.frame / Dplyr ---
time_df_2_2_result <- system.time({
  df_risultato_2_2 <- df_counts_t2 %>%
    group_by(gene) %>%
    mutate(high = ifelse(count > median(count), 1, 0)) %>%
    ungroup()
})
time_df_2_2 <- extract_time(time_df_2_2_result)
cat("DataFrame T2.2 Time:", time_df_2_2, "seconds.\n")


# --- 2.2.2: Data.table ---
time_dt_2_2_result <- system.time({
  dt_counts_t2[, high := ifelse(count > median(count), 1, 0), by = gene]
})
time_dt_2_2 <- extract_time(time_dt_2_2_result)
cat("data.table T2.2 Time:", time_dt_2_2, "seconds.\n")

# Update performance table Task 2.2
performance_comparison_task2 <- rbindlist(list(performance_comparison_task2, 
                                               list("T2.2 Overwrite Flag by Gene", time_df_2_2, time_dt_2_2)), 
                                          use.names = TRUE, fill = TRUE)

# Save the final result modified in dt_counts_t2
# dt_risultato_2 now contains log2_counts AND the overwritten 'high' column
dt_risultato_2 <- dt_counts_t2

print(head(dt_risultato_2, 5))

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Update the main performance table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task2), use.names = TRUE, fill = TRUE)

# Complete Result Task 2
fwrite(dt_risultato_2, file.path(output_dir,"Task2_DerivedColumns_Results.csv"))
cat("\nTask 2 analysis completed. CSV file now contains log2_counts and the overwritten 'high' column.\n")



# ----------------------------------------------------
# TASK 3: Speed Up Joins/Lookups (setkey and Indexes)
# Goal: Benchmark join times on tables with and without key/index.
# ----------------------------------------------------
cat("\n--- START TASK 3: SPEED UP JOINS/LOOKUPS ---\n")

# Initialize a local performance table for Task 3
performance_comparison_task3 <- data.table(Task = character(), 
                                           `Time_DataFrame (s)` = numeric(), 
                                           `Time_DT_NoIndex (s)` = numeric(), 
                                           `Time_DT_Keyed/Index (s)` = numeric()) 

# Re-load clean data (start from original data to measure only the join)
dt_counts_t3 <- fread("bulk_counts_long.csv")
dt_metadata_t3 <- fread("sample_metadata.csv")
df_counts_t3 <- as.data.frame(dt_counts_t3)
df_metadata_t3 <- as.data.frame(dt_metadata_t3)

# ----------------------------------------------------
# TASK 3.1: Equi-Join with setkey()
# Objective: Join metadata into counts using sample_id as the key.
# ----------------------------------------------------
# Compare standard dplyr::left_join against data.table join that uses setkey()
# setkey() pre-sorts the smaller metadata table, allowing data.table to use a much faster "binary search" algorithm for the join

cat("\nRunning Task 3.1: Equi-Join Comparison...\n")

# --- 3.1.1: Data.frame / Dplyr (Join Standard) ---
time_df_3_1_result <- system.time({
  df_risultato_3_1 <- df_counts_t3 %>%
    left_join(df_metadata_t3, by = "sample_id")
})
time_df_3_1 <- extract_time(time_df_3_1_result)
cat("DataFrame T3.1 Time:", time_df_3_1, "seconds.\n")


# --- 3.1.2: Data.table (Join con setkey()) ---
# 1. Define the key on the smaller table (metadata) for optimization.
setkey(dt_metadata_t3, sample_id)

time_dt_key_3_1_result <- system.time({
  # Execute the join with the optimized DT[i] syntax
  dt_risultato_3_1_key <- dt_metadata_t3[dt_counts_t3, on = "sample_id"]
})
time_dt_key_3_1 <- extract_time(time_dt_key_3_1_result)
cat("data.table T3.1 (Keyed) Time:", time_dt_key_3_1, "seconds.\n")


# Update performance table Task 3.1
performance_comparison_task3 <- rbindlist(list(performance_comparison_task3, 
                                               list("T3.1 Equi-Join (Key vs No-Key)", time_df_3_1, NA_real_ , time_dt_key_3_1)), 
                                          use.names = TRUE, fill= TRUE)
print(head(dt_risultato_3_1_key, 5))


# ----------------------------------------------------
# TASK 3.2: Secondary Index and Subset (Lookups)
# Objective: Add a secondary index and measure the time of a query.
# ----------------------------------------------------
# Find a specific set of records before an index is created (forcing a "full table scan") and after.

cat("\nRunning Task 3.2: Secondary Index and Lookup...\n")

# Measure lookup WITHOUT INDEX (uses normal DT subsetting)
gene_subset <- dt_counts_t3[1:10, gene] # Take 10 random genes for the lookup
sample_target <- "S01"
time_no_index_3_2_result <- system.time({
  dt_lookup_no_index <- dt_counts_t3[gene %in% gene_subset & sample_id == sample_target]
})
time_dt_no_index <- extract_time(time_no_index_3_2_result)
cat("data.table T3.2 (No Index) Time:", time_dt_no_index, "seconds.\n")

# Add the secondary index on (gene, sample_id)
setindex(dt_counts_t3, gene, sample_id)

# Measure lookup WITH INDEX
time_with_index_3_2_result <- system.time({
  dt_lookup_with_index <- dt_counts_t3[gene %in% gene_subset & sample_id == sample_target]
})
time_dt_with_index <- extract_time(time_with_index_3_2_result)
cat("data.table T3.2 (With Index) Time:", time_dt_with_index, "seconds.\n")

# Update performance table Task 3.2
performance_comparison_task3 <- rbindlist(list(performance_comparison_task3, 
                                               list("T3.2 Lookup (Index vs No-Index)", NA, time_dt_no_index, time_dt_with_index)), 
                                          fill = TRUE, use.names = TRUE)

print(head(dt_lookup_with_index, 5))

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Task 3 results into the main table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task3), use.names = TRUE, fill = TRUE)

# Saving the complete join result
fwrite(dt_risultato_3_1_key, file.path(output_dir,"Task3_JoinResults_Complete.csv"))
cat("\nTask 3 analysis completed. Cumulative performance files saved.\n")



# ----------------------------------------------------
# TASK 4: Annotate Counts with Patient Info
# Goal: Join Metadata and Find Top 10 Genes per Condition
# ----------------------------------------------------

cat("\n--- START TASK 4: ANNOTATION AND TOP GENES ---\n")

# Initialize a local performance table for Task 4
performance_comparison_task4 <- data.table(Task = character(), 
                                           `Time_DataFrame (s)` = numeric(), 
                                           `Time_DataTable (s)` = numeric(),
                                           `Time_SQL (s)` = numeric()) 

# Re-load data for Task 4 for a clean benchmark
dt_counts_t4 <- fread("bulk_counts_long.csv")
dt_metadata_t4 <- fread("sample_metadata.csv")
df_counts_t4 <- as.data.frame(dt_counts_t4) 
df_metadata_t4 <- as.data.frame(dt_metadata_t4)

# ----------------------------------------------------
# TASK 4.1: Join and Calculate Total Counts Per Patient
# ----------------------------------------------------
# Link samples to patients, then calculates the sum of all counts for each patient.

cat("\nRunning Task 4.1: Total Counts per Patient (3-way comparison)...\n")

# --- 4.1.1: Data.frame / Dplyr ---
time_df_4_1_result <- system.time({
  df_risultato_4_1 <- df_counts_t4 %>%
    # Join to get the patient_id
    left_join(df_metadata_t4, by = "sample_id") %>%
    # Calculate the total sum of counts for each patient
    group_by(patient_id) %>%
    summarise(total_counts = sum(count), 
              .groups = 'drop')
})
time_df_4_1 <- extract_time(time_df_4_1_result)
cat("DataFrame T4.1 Time:", time_df_4_1, "seconds.\n")


# --- 4.1.2: Data.table ---
time_dt_4_1_result <- system.time({
  # Join and sum in a single data.table chain
  dt_risultato_4_1 <- dt_counts_t4[dt_metadata_t4, on = "sample_id"][,
                                                                     list(total_counts = sum(count)),
                                                                     by = patient_id 
  ]
})
time_dt_4_1 <- extract_time(time_dt_4_1_result)
cat("data.table T4.1 Time:", time_dt_4_1, "seconds.\n")


# --- 4.1.3: SQL (sqldf) ---
time_sql_4_1_result <- system.time({
  sql_risultato_4_1 <- sqldf::sqldf("
    SELECT T2.patient_id, SUM(T1.count) AS total_counts
    FROM df_counts_t4 AS T1
    LEFT JOIN df_metadata_t4 AS T2 ON T1.sample_id = T2.sample_id
    GROUP BY T2.patient_id
  ")
})
time_sql_4_1 <- extract_time(time_sql_4_1_result)
cat("SQL T4.1 Time:", time_sql_4_1, "seconds.\n")

# Update performance table Task 4.1
performance_comparison_task4 <- rbindlist(list(performance_comparison_task4, 
                                               list("T4.1 Total Counts per Patient", time_df_4_1, time_dt_4_1, time_sql_4_1)), 
                                          use.names = TRUE, fill= TRUE)

print(head(df_risultato_4_1, 5))


# ----------------------------------------------------
# TASK 4.2: Find Top 10 Genes by Mean Count per Condition
# ----------------------------------------------------
# Find the mean count for each gene within each condition, then rank these means to find the top 10 highest-expressing genes.

cat("\nRunning Task 4.2: Top 10 Genes by Condition...\n")

# --- 4.2.1: Data.frame / Dplyr ---
time_df_4_2_result <- system.time({
  df_risultato_4_2 <- df_counts_t4 %>%
    left_join(df_metadata_t4, by = "sample_id") %>%
    # 1. Calculate mean per (gene, condition)
    group_by(gene, condition) %>%
    summarise(mean_count = mean(count), .groups = 'drop_last') %>%
    # 2. Order by condition and then by count in descending order
    arrange(condition, desc(mean_count)) %>%
    # 3. Take only the first 10 for each condition group
    slice_head(n = 10) %>%
    ungroup()
})
time_df_4_2 <- extract_time(time_df_4_2_result)
cat("DataFrame T4.2 Time:", time_df_4_2, "seconds.\n")


# --- 4.2.2: Data.table ---
time_dt_4_2_result <- system.time({
  # 1. Calculate mean per (gene, condition) and save it in a temporary variable
  dt_ranked <- dt_counts_t4[dt_metadata_t4, on = "sample_id"][,
                                                              list(mean_count = mean(count)),
                                                              by = list(gene, condition)
  ]
  # 2. Assign rank and filter
  dt_risultato_4_2 <- dt_ranked[
    # Assign a rank for each condition group, rank is a function of data.table
    , rank := rank(-mean_count, ties.method = "min"), 
    by = condition
  ][
    # Filter only those with rank <= 10
    rank <= 10
  ]
})
time_dt_4_2 <- extract_time(time_dt_4_2_result)
cat("data.table T4.2 Time:", time_dt_4_2, "seconds.\n")

# --- 4.2.3: SQL (sqldf) ---
# SQL implementation uses Window Functions for ranking
time_sql_4_2_result <- system.time({
  sql_risultato_4_2 <- sqldf::sqldf("
    WITH RankedGenes AS (
        SELECT 
            T1.gene, 
            T2.condition, 
            AVG(T1.count) AS mean_count,
            -- Assegna un rank per ogni condizione, ordinando per mean_count decrescente
            RANK() OVER (PARTITION BY T2.condition ORDER BY AVG(T1.count) DESC) as rn
        FROM df_counts_t4 AS T1
        LEFT JOIN df_metadata_t4 AS T2 ON T1.sample_id = T2.sample_id
        GROUP BY T1.gene, T2.condition
    )
    SELECT gene, condition, mean_count 
    FROM RankedGenes
    WHERE rn <= 10;
  ")
})
time_sql_4_2 <- extract_time(time_sql_4_2_result)
cat("SQL T4.2 Time:", time_sql_4_2, "seconds.\n")


# Update performance table Task 4.2
performance_comparison_task4 <- rbindlist(list(performance_comparison_task4, 
                                               list("T4.2 Top 10 Genes by Condition", time_df_4_2, time_dt_4_2, time_sql_4_2)), 
                                          use.names = TRUE, fill= TRUE)

print(dt_risultato_4_2, nrows = 20)

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Task 4 performance table into the main cumulative table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task4), use.names = TRUE, fill = TRUE)

# Saving the complete result of Task 4 (Top 10)
fwrite(dt_risultato_4_2, file.path(output_dir, "Task4_Top10GenesByCondition_Results.csv"))

cat("\nTask 4 analysis completed. Cumulative performance files saved.\n")



# ----------------------------------------------------
# TASK 5: Classify Values Against Reference Ranges
# Goal: Perform Non-Equi Join and Count Abnormal Rates.
# ----------------------------------------------------

cat("\n--- START TASK 5: CLINICAL NON-EQUI JOIN ---\n")

# Initialize a local performance table for Task 5
performance_comparison_task5 <- data.table(Task = character(), 
                                           `Time_DataFrame (s)` = numeric(), 
                                           `Time_DataTable (s)` = numeric(),
                                           `Time_SQL (s)` = numeric()) 
# Load data for Task 5
dt_labs_t5 <- fread("clinical_labs.csv")
dt_ranges_t5 <- fread("lab_reference_ranges.csv")

# Rename the 'lab' column to 'lab_name' in both tables for a successful join
setnames(dt_labs_t5, old = "lab", new = "lab_name", skip_absent = TRUE)
setnames(dt_ranges_t5, old = "lab", new = "lab_name", skip_absent = TRUE)

# Conversion to data.frame for dplyr and sqldf benchmark
df_labs_t5 <- as.data.frame(dt_labs_t5)
df_ranges_t5 <- as.data.frame(dt_ranges_t5)

# ----------------------------------------------------
# TASK 5.1: Non-Equi Join and Labeling (Base of Analysis)
# ----------------------------------------------------
# Label each lab as "normal" vs "out_of_range" using a non-equi join

cat("\nRunning Task 5.1: Non-Equi Join and Labeling (3-way comparison)...\n")

# --- 5.1.1: Data.frame / Dplyr ---
time_df_5_1_result <- system.time({
  # Step 1: Join on lab_name (We add relationship = "many-to-many" to silence the warning)
  df_joined_t5 <- df_labs_t5 %>%
    left_join(df_ranges_t5, by = "lab_name", relationship = "many-to-many") %>%
    # Step 2: Classification (Non-Equi Logic)
    mutate(is_normal = (value >= lower & value <= upper),
           status = ifelse(is_normal, "Normal", "Out_of_Range")) %>%
    select(patient_id, lab_name, value, status) 
})
time_df_5_1 <- extract_time(time_df_5_1_result)
cat("DataFrame T5.1 Time:", time_df_5_1, "seconds.\n")

# --- 5.1.2: Data.table (Non-Equi Join) ---
time_dt_5_1_result <- system.time({
  # 1. Define the key on the range limits for efficient alignment
  setkey(dt_labs_t5, lab_name)
  setkey(dt_ranges_t5, lab_name)
  
  # 2. Execute the Non-Equi Join in a single line: 
  # The 'on' conditions define the intervals [lower, upper]
  dt_risultato_5_1 <- dt_labs_t5[dt_ranges_t5, 
                                 on = .(lab_name, value >= lower, value <= upper),
                                 # nomatch = 0 takes only results that fall within the interval, if no match it is discarded
                                 nomatch = 0
  ][,
    .(patient_id, lab_name, value, status = "Normal")
  ]
})
time_dt_5_1 <- extract_time(time_dt_5_1_result)
cat("data.table T5.1 Time:", time_dt_5_1, "seconds.\n")


# --- 5.1.3: SQL (sqldf) ---
# SQL supports the non-equi join using the WHERE clause
# CASE is similar to ifelse
time_sql_5_1_result <- system.time({
  sql_risultato_5_1 <- sqldf::sqldf("
    SELECT T1.patient_id, T1.lab_name, T1.value, 
           CASE 
               WHEN T1.value >= T2.lower AND T1.value <= T2.upper THEN 'Normal'
               ELSE 'Out_of_Range'
           END AS status
    FROM df_labs_t5 AS T1
    LEFT JOIN df_ranges_t5 AS T2 ON T1.lab_name = T2.lab_name
  ")
})
time_sql_5_1 <- extract_time(time_sql_5_1_result)
cat("SQL T5.1 Time:", time_sql_5_1, "seconds.\n")

# Update performance table Task 5.1
performance_comparison_task5 <- rbindlist(list(performance_comparison_task5, 
                                               list("T5.1 Non-Equi Join Label", time_df_5_1, time_dt_5_1, time_sql_5_1)), 
                                          use.names = TRUE, fill= TRUE)
print(head(df_joined_t5, 5))


# ----------------------------------------------------
# TASK 5.2: Count Abnormal Rates Per Patient and Per Lab
# ----------------------------------------------------
# Calculate the total and abnormal test rates for each patient/lab combination.

cat("\nRunning Task 5.2: Summarize Abnormal Rates...\n")

# We create a Base table (dt_final_t5_base) for Task 5.2 from the SQL result (which contains the status), ensuring the same starting point for all three comparisons.
dt_final_t5_base <- as.data.table(sql_risultato_5_1) 
# Copy of the Base result for the data.frame (dplyr) approach
df_final_t5_base_df <- as.data.frame(dt_final_t5_base)

# --- 5.2.1: Data.frame / Dplyr ---
time_df_5_2_result <- system.time({
  df_risultato_5_2 <- df_final_t5_base_df %>% 
    # Calculate the total tests (TotalTests) and count of Out_of_Range tests (AbnormalTests)
    group_by(patient_id, lab_name) %>%
    summarise(
      TotalTests = n(),
      AbnormalTests = sum(status == "Out_of_Range"),
      AbnormalRate = AbnormalTests / TotalTests,
      .groups = 'drop'
    )
})
time_df_5_2 <- extract_time(time_df_5_2_result)
cat("DataFrame T5.2 Time:", time_df_5_2, "seconds.\n")


# --- 5.2.2: Data.table ---
time_dt_5_2_result <- system.time({
  dt_risultato_5_2 <- dt_final_t5_base[, 
                                       list(
                                         TotalTests = .N, # .N is the equivalent of n()
                                         AbnormalTests = sum(status == "Out_of_Range"),
                                         AbnormalRate = sum(status == "Out_of_Range") / .N
                                       ),
                                       by = list(patient_id, lab_name)
  ]
})
time_dt_5_2 <- extract_time(time_dt_5_2_result)
cat("data.table T5.2 Time:", time_dt_5_2, "seconds.\n")


# --- 5.2.3: SQL (sqldf) ---
time_sql_5_2_result <- system.time({
  sql_risultato_5_2 <- sqldf::sqldf("
    SELECT patient_id, lab_name, 
           COUNT(*) AS TotalTests,
           SUM(CASE WHEN status = 'Out_of_Range' THEN 1 ELSE 0 END) AS AbnormalTests,
           CAST(SUM(CASE WHEN status = 'Out_of_Range' THEN 1 ELSE 0 END) AS REAL) / COUNT(*) AS AbnormalRate
    FROM sql_risultato_5_1
    GROUP BY patient_id, lab_name
  ")
})
time_sql_5_2 <- extract_time(time_sql_5_2_result)
cat("SQL T5.2 Time:", time_sql_5_2, "seconds.\n")

# Update performance table Task 5.2
performance_comparison_task5 <- rbindlist(list(performance_comparison_task5, 
                                               list("T5.2 Summarize Abnormal Rates", time_df_5_2, time_dt_5_2, time_sql_5_2)), 
                                          use.names = TRUE, fill= TRUE)

print(head(df_risultato_5_2, 5))

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Task 5 results into the main cumulative table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task5), use.names = TRUE, fill = TRUE)

# Saving the complete result (Task 5.2)
fwrite(dt_risultato_5_2, file.path(output_dir, "Task5_AbnormalRates_Results.csv"))

cat("\nTask 5 analysis completed. Cumulative performance files saved.\n")



# ----------------------------------------------------
# TASK 6: Nearest-Time Matching (Rolling Join)
# Goal: Match Lab vs Vitals by nearest time, calculate lag and correlation.
# ----------------------------------------------------

cat("\n--- START TASK 6: NEAREST-TIME ROLLING JOIN & CORRELATION ---\n")

# Load necessary library for Task 6
library(lubridate)  # Contains the ymd_hms() function used for date/time conversion.

# Initialize a local performance table for Task 6
performance_comparison_task6 <- data.table(Task = character(),
                                           `Time_DT_Rolling (s)` = numeric(),
                                           `Time_DataFrame (s)` = numeric(), 
                                           `Time_DataTable (s)` = numeric()
                                           ) 
# Data Loading for Task 6 
# Note: dt_labs_t5 is used as the base and copied, it was loaded in Task 5 prep.
dt_labs_t6 <- copy(dt_labs_t5) 
dt_vitals_t6 <- fread("vitals_time_series.csv")

# --- TEMPORAL DATA PREPARATION ---
cat("Preparing temporal data and Pivot from Long to Wide...\n")
  # 1. Convert text to datetime object and sets the timezone as UTC (Universal Coordinated Time)
  if ("time_iso" %in% names(dt_labs_t6)) { dt_labs_t6[, time_iso := ymd_hms(as.character(time_iso), tz = "UTC")] }
  if ("time_iso" %in% names(dt_vitals_t6)) { dt_vitals_t6[, time_iso := ymd_hms(as.character(time_iso), tz = "UTC")] }
  
  # 2. PIVOT: Transform vitals from 'long' to 'wide' format
  # fun.aggregate = mean added to handle duplicates
  # dcast is data.table's "pivot wide" command.
  dt_vitals_t6_wide <- dcast(dt_vitals_t6, patient_id + time_iso ~ vital, value.var = "value", fun.aggregate = mean, na.rm = TRUE) 
  # Order data, which is mandatory for rolling join to work
  setorder(dt_vitals_t6_wide, patient_id, time_iso) 

# Expected vital columns after Pivot
expected_vitals_cols <- c("HR", "SBP") 

# ----------------------------------------------------
# TASK 6.1: Rolling Join (Nearest Match and Lag Calculation)
# ----------------------------------------------------
# Take the lab draw times and 'rolls' them back to find the last available vital sign time, then calculate lag_minutes.

cat("\nRunning Task 6.1: Rolling Join\n")

# Rename the time column in 'x' (vitals) BEFORE the rolling join to preserve both times (from "time_iso" to "vitals_time") ---
setnames(dt_vitals_t6_wide, "time_iso", "vitals_time", skip_absent = TRUE)

# --- 6.1.2: Data.table (Optimized Rolling Join) ---
time_dt_rolling_6_1_result <- system.time({
  
  # Preserve 'vitals_time' in a copy, because the join will overwrite the key.
  dt_vitals_t6_wide[, vitals_time_matched := vitals_time]
  
  # Execute Join (dt_vitals_t6_wide[dt_labs_t6])
  # x (vitals) has: [patient_id, vitals_time (key), vitals_time_matched (copy), HR, SBP]
  # i (labs)   has: [patient_id, time_iso (key), lab_name, value]
  # on = .(col_x = col_i)
  # Join patient_id (from x) with patient_id (from i), and join vitals_time (from x) with time_iso (from i).
  
  dt_join_result <- dt_vitals_t6_wide[dt_labs_t6, 
                                      on = .(patient_id, vitals_time = time_iso), 
                                      roll = TRUE, mult = "last"] # mult = "last" takes only the last match (nearest in time)
  # Rename TIME columns after join
  setnames(dt_join_result, 
           old = c("vitals_time", "vitals_time_matched"), 
           new = c("lab_time", "vitals_time"), 
           skip_absent = TRUE)
  
  # Calculate the lag
  if ("lab_time" %in% names(dt_join_result) && "vitals_time" %in% names(dt_join_result)) {
    # Ensure both are POSIXct for safe subtraction
    dt_join_result[, lab_time := ymd_hms(lab_time, tz = "UTC")]
    dt_join_result[, vitals_time := ymd_hms(vitals_time, tz = "UTC")]
    
    dt_join_result[, lag_minutes := as.numeric(difftime(lab_time, vitals_time, units = "mins"))]
  } else { 
    dt_join_result[, lag_minutes := NA_real_] 
    warning("Lag calculation failed: 'lab_time' or 'vitals_time' missing after join.")
  }
  
  # Rename final columns (HR, SBP -> nearest_hr/sbp, value -> lab_value)
  value_col_name <- if ("value" %in% names(dt_join_result)) "value" else (if ("i.value" %in% names(dt_join_result)) "i.value" else NULL)
  cols_to_rename_old <- c("HR", "SBP")
  cols_to_rename_new <- c("nearest_hr", "nearest_sbp")
  if (!is.null(value_col_name) && value_col_name %in% names(dt_join_result) ) {
    cols_to_rename_old <- c(cols_to_rename_old, value_col_name)
    cols_to_rename_new <- c(cols_to_rename_new, "lab_value")
  }
  setnames(dt_join_result, old = cols_to_rename_old, new = cols_to_rename_new, skip_absent=TRUE)
  dt_risultato_6_1 <- dt_join_result
})

time_dt_rolling_6_1 <- extract_time(time_dt_rolling_6_1_result)
cat("data.table T6.1 Time (Rolling Join):", time_dt_rolling_6_1, "seconds.\n")

print(head(dt_risultato_6_1, 5))

# Update performance table Task 6.1 
performance_comparison_task6 <- rbindlist(list(performance_comparison_task6,
                                               list(Task= "T6.1 Rolling Join", `Time_DT_Rolling (s)` = time_dt_rolling_6_1)),
                                          fill = TRUE, use.names = TRUE)


# ----------------------------------------------------
# TASK 6.2: Summarize CRP vs HR/SBP Correlation Per Patient (DT/Dplyr)
# ----------------------------------------------------
# Filter for 'CRP' labs and correlate the lab values against the nearest Heart Rate (HR) and Systolic Blood Pressure (SBP) found in T6.1, grouped by each patient.

cat("\nRunning Task 6.2: Summarize Correlation (data.table focus)...\n")

# Data base for Task 6.2
required_cols_cor <- c("nearest_hr", "nearest_sbp", "lab_name", "lab_value")
if (all(required_cols_cor %in% names(dt_risultato_6_1))) {
  # Filter for rows with non-missing vitals and lab_name == "CRP"
  dt_cor_base <- dt_risultato_6_1[!is.na(nearest_hr) & !is.na(nearest_sbp) & lab_name == "CRP" & !is.na(lab_value)] # !is.na means that is not NA
  df_cor_base <- as.data.frame(dt_cor_base)
} else {
  warning("Required columns for Task 6.2 (nearest_hr/sbp) missing from join result. Skipping correlation.")
  dt_cor_base <- data.table() #create a table with filtered dates 
  df_cor_base <- data.frame() #create a table in data.frame format
}

# --- 6.2.1: Data.frame / Dplyr ---
time_df_6_2_result <- system.time({
  if (nrow(df_cor_base) > 0 && "lab_value" %in% names(df_cor_base) && "nearest_hr" %in% names(df_cor_base)) { 
    # Counts how many rows (CRP exams) each patient has and keeps only patients with more than two exams
    valid_patients_df <- dt_cor_base %>% count(patient_id) %>% filter(n >= 2) %>% pull(patient_id) 
    if (length(valid_patients_df) > 0) {
      df_risultato_6_2 <- df_cor_base %>%
        filter(patient_id %in% valid_patients_df) %>%
        group_by(patient_id) %>%
        summarise(
          # Calculate correlation between lab_value and nearest HR
          # pairwise.complete.obs ignores NA pairs
          cor_CRP_HR = cor(lab_value, nearest_hr, use = "pairwise.complete.obs"), 
          # Calculate correlation between lab_value and nearest SBP
          cor_CRP_SBP = cor(lab_value, nearest_sbp, use = "pairwise.complete.obs"), 
          .groups = 'drop'
        )
    } else { df_risultato_6_2 <- data.frame(patient_id=character(), cor_CRP_HR=numeric(), cor_CRP_SBP=numeric()); cat("Nessun paziente con dati sufficienti (Dplyr).\n") }
  } else { df_risultato_6_2 <- data.frame(patient_id=character(), cor_CRP_HR=numeric(), cor_CRP_SBP=numeric()); cat("Dati base per correlazione mancanti (Dplyr).\n") }
})
time_df_6_2 <- extract_time(time_df_6_2_result)
cat("DataFrame T6.2 Time:", time_df_6_2, "seconds.\n")


# --- 6.2.2: Data.table ---
time_dt_6_2_result <- system.time({
  if (nrow(dt_cor_base) > 0 && "lab_value" %in% names(dt_cor_base) && "nearest_hr" %in% names(dt_cor_base)) {
    dt_risultato_6_2 <- dt_cor_base[, .N, by = patient_id][N >= 2][dt_cor_base, on = "patient_id"][, #.N to count rows and >=2 to keep only rows > 2
                                                                                                   .(cor_CRP_HR = cor(lab_value, nearest_hr, use = "pairwise.complete.obs"),
                                                                                                     cor_CRP_SBP = cor(lab_value, nearest_sbp, use = "pairwise.complete.obs")),
                                                                                                   by = patient_id
    ]
    if (nrow(dt_risultato_6_2) == 0) { cat("No patients with sufficient data (data.table).\n") }
  } else { dt_risultato_6_2 <- data.table(); cat("Base data for correlation missing (data.table).\n") }
})
time_dt_6_2 <- extract_time(time_dt_6_2_result)
cat("data.table T6.2 Time:", time_dt_6_2, "seconds.\n")

# Update performance table Task 6.2
performance_comparison_task6 <- rbindlist(list(performance_comparison_task6,
                                               list(Task= "T6.2 Summarize Correlation", 
                                                    `Time_DataFrame (s)` = time_df_6_2, 
                                                    `Time_DataTable (s)` = time_dt_6_2)),
                                          fill = TRUE, use.names = TRUE) 
print(head(df_risultato_6_2, 5))

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Task 6 results into the main cumulative table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task6), use.names = TRUE, fill = TRUE)

fwrite(dt_risultato_6_1, file.path(output_dir,"Task6_RollingJoinResults_Complete.csv"))
fwrite(dt_risultato_6_2, file.path(output_dir,"Task6_CorrelationResults_Complete.csv"))

cat("\nTask 6 analysis completed. Cumulative performance files saved.\n")



# ----------------------------------------------------
# TASK 7: Efficient Genomic Window Filtering
# Goal: Extract peaks on chr2 (2-4Mb) and find the Top 50 by score.
# ----------------------------------------------------
cat("\n--- START TASK 7: ATAC-SEQ GENOMIC FILTERING ---\n")

# Initialize a local performance table for Task 7
performance_comparison_task7 <- data.table(Task = character(),
                                           `Time_DataFrame (s)` = numeric(),
                                           `Time_DataTable (s)` = numeric()
                                           ) 
# Data Loading for Task 7
dt_peaks_t7 <- fread("atac_peaks.bed.csv")
df_peaks_t7 <- as.data.frame(dt_peaks_t7)
  
# Define the limits for the filter as variables
chr_target <- "chr2"
start_min <- 2000000
start_max <- 4000000

# ----------------------------------------------------
# TASK 7.1: Filter by Genomic Region
# ----------------------------------------------------
# Find all peaks that fall on chr2 and have a start position between 2,000,000 and 4,000,000.

cat("\nRunning Task 7.1: Filter by Genomic Region (3-way comparison)...\n")

# --- 7.1.1: Data.frame / Dplyr ---
time_df_7_1_result <- system.time({
  df_filtered_t7 <- df_peaks_t7 %>%
    filter(chr == chr_target & start >= start_min & start <= start_max)
})
time_df_7_1 <- extract_time(time_df_7_1_result)
cat("DataFrame T7.1 Time:", time_df_7_1, "seconds.\n")


# --- 7.1.2: Data.table ---
time_dt_7_1_result <- system.time({
  # Use the df_peaks_t7 copy
  dt_filtered_t7 <- dt_peaks_t7[chr == chr_target & start >= start_min & start <= start_max]
})
time_dt_7_1 <- extract_time(time_dt_7_1_result)
cat("data.table T7.1 Time:", time_dt_7_1, "seconds.\n")

# Update performance table Task 7.1
performance_comparison_task7 <- rbindlist(list(performance_comparison_task7,
                                               list("T7.1 Filter Genomic Region", time_df_7_1, time_dt_7_1)),
                                          use.names = TRUE, fill = TRUE) 

print(head(dt_filtered_t7, 5))

# ----------------------------------------------------
# TASK 7.2: Find Top 50 Peaks by Score
# ----------------------------------------------------

# Take results from T7.1, order them by the score and select the top 50.
cat("\nRunning Task 7.2: Find Top 50 Peaks by Score (3-way comparison)...\n")

# --- 7.2.1: Data.frame / Dplyr ---
# Use the df_filtered_t7 result
time_df_7_2_result <- system.time({
  df_top50_t7 <- df_filtered_t7 %>%
    arrange(desc(score)) %>%
    slice_head(n = 50)
})
time_df_7_2 <- extract_time(time_df_7_2_result)
cat("DataFrame T7.2 Time:", time_df_7_2, "seconds.\n")


# --- 7.2.2: Data.table ---
# Use the dt_filtered_t7 result
dt_filtered_t7_copy <- copy(dt_filtered_t7) # Copy is needed because setorder modifies in-place
time_dt_7_2_result <- system.time({
  setorder(dt_filtered_t7_copy, -score) # Order in-place; -score means descending
  dt_top50_t7 <- head(dt_filtered_t7_copy, 50) # Take the first 50
})
time_dt_7_2 <- extract_time(time_dt_7_2_result)
cat("data.table T7.2 Time (setorder):", time_dt_7_2, "seconds.\n")

# Update performance table Task 7.2
performance_comparison_task7 <- rbindlist(list(performance_comparison_task7,
                                               list("T7.2 Top 50 by Score", time_df_7_2, time_dt_7_2)),
                                          use.names = TRUE, fill = TRUE)

print(head(dt_top50_t7, 5))

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Task 7 results into the main cumulative table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task7), use.names = TRUE, fill = TRUE)

# Saving the complete resulto of Task 7 
fwrite(dt_top50_t7, file.path(output_dir, "Task7_Top50Peaks_Results.csv"))

cat("\nTask 7 analysis completed. Cumulative performance files saved.\n")



# ----------------------------------------------------
# TASK 8: Multi-Column Operations per Group (Robust Statistics)
# Goal: Calculate mean, median, Q1/Q3 per (gene, condition) and filter.
# ----------------------------------------------------
cat("\n--- START TASK 8: ROBUST STATISTICS AND FILTER ---\n")

# Initialize a local performance table for Task 8
performance_comparison_task8 <- data.table(Task = character(),
                                           `Time_DataFrame (s)` = numeric(),
                                           `Time_DataTable (s)` = numeric()) 

# Use data already loaded at the start of the script (dt_counts, dt_metadata)
dt_counts_t8 <- copy(dt_counts)
dt_metadata_t8 <- copy(dt_metadata)
df_counts_t8 <- as.data.frame(dt_counts_t8)
df_metadata_t8 <- as.data.frame(dt_metadata_t8)

# Preparation: Join data to group by 'condition'
df_joined_t8 <- df_counts_t8 %>% left_join(df_metadata_t8, by = "sample_id")
dt_joined_t8 <- dt_counts_t8[dt_metadata_t8, on = "sample_id"]

# ----------------------------------------------------
# TASK 8.1: Calculate Robust Statistics (Mean, Median, Q1/Q3)
# ----------------------------------------------------
# Group by both gene and condition and calculate mean, median, Q1, and Q3.

cat("\nRunning Task 8.1: Calculate Robust Stats (3-way comparison)...\n")

# --- 8.1.1: Data.frame / Dplyr ---
time_df_8_1_result <- system.time({
  df_stats_t8 <- df_joined_t8 %>%
    group_by(gene, condition) %>%
    summarise(
      mean_count = mean(count),
      median_count = median(count),
      Q1 = quantile(count, 0.25),
      Q3 = quantile(count, 0.75),
      .groups = 'drop'
    )
})
time_df_8_1 <- extract_time(time_df_8_1_result)
cat("DataFrame T8.1 Time:", time_df_8_1, "seconds.\n")


# --- 8.1.2: Data.table ---
time_dt_8_1_result <- system.time({
  dt_stats_t8 <- dt_joined_t8[,
                              list(
                                mean_count = mean(count),
                                median_count = median(count),
                                Q1 = quantile(count, 0.25),
                                Q3 = quantile(count, 0.75)
                              ),
                              by = list(gene, condition)
  ]
})
time_dt_8_1 <- extract_time(time_dt_8_1_result)
cat("data.table T8.1 Time:", time_dt_8_1, "seconds.\n")

# Update performance table Task 8.1
performance_comparison_task8 <- rbindlist(list(performance_comparison_task8,
                                               list("T8.1 Calculate Robust Stats", time_df_8_1, time_dt_8_1)),
                                          use.names = TRUE, fill = TRUE)

print(head(dt_stats_t8, 5))


# ----------------------------------------------------
# TASK 8.2: Filter by Mean Difference (treated > control)
# ----------------------------------------------------
# Filter for genes whith treated mean higher than control mean.

cat("\nRunning Task 8.2: Filter by Treated > Control Mean...\n")

# --- 8.2.1: Data.frame / Dplyr ---
time_df_8_2_result <- system.time({
  df_risultato_8_2 <- df_stats_t8 %>%
    pivot_wider(names_from = condition,  # Names of condition ("treated" and "control") are now names of new columns
                values_from = c(mean_count, median_count, Q1, Q3)) %>%  # Take the values of these four columns and place them under Treated and Control (new columns)
    filter(mean_count_treated > mean_count_control)
})
time_df_8_2 <- extract_time(time_df_8_2_result)
cat("DataFrame T8.2 Time:", time_df_8_2, "seconds.\n")


# --- 8.2.2: Data.table ---
time_dt_8_2_result <- system.time({
  # 1. Transform the table into wide format (pivot)
  dt_wide <- dcast(dt_stats_t8, 
                   gene ~ condition, # one row per gene and columns created by condition values
                   value.var = c("mean_count", "median_count", "Q1", "Q3")) # Tell 'value.var' to rotate ALL statistic columns.
  # 2. Filter for treated > control
  dt_risultato_8_2 <- dt_wide[mean_count_treated > mean_count_control]
})
time_dt_8_2 <- extract_time(time_dt_8_2_result)
cat("data.table T8.2 Time:", time_dt_8_2, "seconds.\n")


# Update performance table Task 8.2
performance_comparison_task8 <- rbindlist(list(performance_comparison_task8,
                                               list("T8.2 Filter (Treated > Control)", time_df_8_2, time_dt_8_2)),
                                          use.names = TRUE, fill = TRUE)

print(head(dt_risultato_8_2, 5))

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Task 8 results into the main cumulative table 
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task8), fill = TRUE)

fwrite(dt_risultato_8_2, file.path(output_dir, "Task8_TreatedVsControl_Results.csv"))
cat("\nTask 8 analysis completed. Cumulative performance files saved.\n")



# ----------------------------------------------------
# TASK 9: Data Transformation (Wide -> Long -> Totals -> Wide)
# Goal: Pivot the matrix, calculate per-sample totals, and aggregate by condition.
# ----------------------------------------------------
cat("\n--- START TASK 9: WIDE -> LONG -> WIDE TRANSFORMATION ---\n")

# Initialize a local performance table for Task 9
performance_comparison_task9 <- data.table(Task = character(),
                                           `Time_DataFrame (s)` = numeric(), 
                                           `Time_DataTable (s)` = numeric()) 

# --- DATA LOADING AND LOCAL PREPARATION FOR TASK 9 ---
dt_counts_wide_t9_local <- fread("bulk_counts_wide.csv")
dt_metadata_t9_local <- fread("sample_metadata.csv")

df_counts_wide_t9_local <- as.data.frame(dt_counts_wide_t9_local)
df_metadata_t9_local <- as.data.frame(dt_metadata_t9_local) 

setkey(dt_metadata_t9_local, sample_id) # Prepare the key for the data.table join

# Indices of sample columns (all columns except the first 'gene')
SAMPLE_COLS_INDICES <- 2:ncol(df_counts_wide_t9_local)  # Creates a vector from 2 to the number of columns (i.e., samples)

# --- TASK 9.1: Wide -> Long -> Totals -> Wide ---
# 3-step process (Wide-to-Long, Add Totals, Aggregate-to-Wide)

cat("\n--- Running Task 9.1: Transformation Pipeline...\n")

# --- 9.1.1: Data.frame / Tidyverse (pivot_longer + group_by + pivot_wider) ---
time_df_9_1_result <- system.time({
  # 1. Wide -> Long
  df_long <- df_counts_wide_t9_local %>% 
    pivot_longer(cols = all_of(SAMPLE_COLS_INDICES), # Select sample columns (S01, S02...)
                 names_to = "sample_id", 
                 values_to = "count")
  
    # 2. Add per-sample totals
  df_with_totals <- df_long %>%
    group_by(sample_id) %>%
    # Sum all counts for each sample
    mutate(sample_total = sum(count, na.rm = TRUE)) %>%
    ungroup()
  
    # 3. Aggregate by (gene, condition) and Long -> Wide
  df_final_wide <- df_with_totals %>%
    left_join(df_metadata_t9_local, by = "sample_id") %>%   #Join on sample_id 
    group_by(gene, condition) %>%                           # Aggregate by (gene, condition)
    summarise(mean_count = mean(count, na.rm = TRUE), .groups = 'drop') %>%  
    # Long -> Wide (gene x condition)
    pivot_wider(names_from = condition, values_from = mean_count)
  
  df_risultato_9_1 <- df_final_wide
})

time_df_9_1 <- extract_time(time_df_9_1_result)
cat("DataFrame T9.1 Time (Tidyverse):", time_df_9_1, "seconds.\n")


# --- 9.1.2: Data.table (melt + dcast) ---
time_dt_9_1_result <- system.time({
  # 1. Wide -> Long (melt)
  dt_long <- melt(dt_counts_wide_t9_local, 
                  id.vars = "gene",        # Keep "gene" as fixed column
                  measure.vars = SAMPLE_COLS_INDICES, 
                  variable.name = "sample_id", 
                  value.name = "count")
  
  # 2. Add per-sample totals
  dt_sample_totals <- dt_long[, .(sample_total = sum(count, na.rm = TRUE)), by = sample_id]
  dt_long_with_totals <- dt_sample_totals[dt_long, on = "sample_id"]
  
  # 3. Aggregate by (gene, condition) and Long -> Wide
  dt_joined <- dt_metadata_t9_local[dt_long_with_totals, on = "sample_id"] #Join 
  #fun.aggregate = mean to calculate the mean and dcast for the pivot
  dt_risultato_9_1 <- dcast(dt_joined, gene ~ condition, value.var = "count", fun.aggregate = mean)  #fun.aggregate = mean to calculate the mean and dcast for the pivot
})
time_dt_9_1 <- extract_time(time_dt_9_1_result)
cat("data.table T9.1 Time (melt/dcast):", time_dt_9_1, "seconds.\n")

print(head(dt_long, 5))
print(head(dt_long_with_totals, 5))
print(head(dt_risultato_9_1, 5))

# Update performance table Task 9.1
performance_comparison_task9 <- rbindlist(list(performance_comparison_task9,
                                         list("T9.1 Wide->Long->Wide", time_df_9_1, time_dt_9_1)),
                                    use.names = TRUE, fill = TRUE)

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Task 9 results into the main cumulative table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task9), fill = TRUE)

# Saving complete result of Task 9
fwrite(dt_risultato_9_1, file.path(output_dir, "Task9_GeneXCondition_Results.csv"))

cat("\nTask 9 analysis completed. Cumulative performance files saved.\n")




# ----------------------------------------------------
# TASK 10: ATAC-to-Gene Mapping (Spatial Join)
# Goal: Map ATAC peaks to gene regions and calculate overlap length.
# ----------------------------------------------------

cat("\n--- START TASK 10: ATAC-TO-GENE SPATIAL JOIN ---\n")

# Initialize a local performance table for Task 10
performance_comparison_task10 <- data.table(Task = character(),
                                            `Time_DataFrame (s)` = numeric(),
                                            `Time_DataTable (s)` = numeric(),
                                            `Time_SQL (s)` = numeric()) 

# --- DATA LOADING AND LOCAL PREPARATION FOR TASK 10 ---
dt_peaks_t10 <- fread("atac_peaks.bed.csv")
dt_genes_t10 <- fread("gene_annotation.bed.csv") 

# Assume the correct order for peaks: chr, start, end, peaks_id, score
setnames(dt_peaks_t10, c("V1", "V2", "V3", "V4", "V5"), c("chr", "start", "end", "peaks_id", "score"), skip_absent = TRUE)

# Assume dt_genes_t10 has columns: chr, start, end, gene
setnames(dt_genes_t10, c("V1", "V2", "V3", "V4"), c("chr", "start", "end", "gene"), skip_absent = TRUE)

df_genes_t10 <- as.data.frame(dt_genes_t10)
df_peaks_t10 <- as.data.frame(dt_peaks_t10)


# Function to calculate the overlap length
calc_overlap_length <- function(start1, end1, start2, end2) {
  pmax(0, pmin(end1, end2) - pmax(start1, start2))
}


# ----------------------------------------------------
# TASK 10.1: Spatial Join (Intersection) and Peak Count per Gene
# ----------------------------------------------------
cat("\nRunning Task 10.1: Spatial Join & Peak Count (DT vs Dplyr)...\n")

# --- 10.1.1: Data.frame / Dplyr  ---
time_df_10_1_result <- system.time({
  # Execute join on chr AND the overlap condition (non-equi join)
  df_overlapped <- df_peaks_t10 %>%
    left_join(df_genes_t10, by = "chr", relationship = "many-to-many") %>% # Join on chromosome
    filter(start.x <= end.y & end.x >= start.y) # Overlap condition
  
  # Count peaks per gene
  df_count <- df_overlapped %>%
    group_by(gene) %>% 
    summarise(peak_count = n(), .groups = 'drop') # Group results by gene and count overlaps found
  df_risultato_10_1 <- df_count
})
time_df_10_1 <- extract_time(time_df_10_1_result)
cat("DataFrame T10.1 Time:", time_df_10_1, "seconds.\n")


# --- 10.1.2: Data.table ---
time_dt_10_1_result <- system.time({
  # 1. Prepare keys/indices
  setkey(dt_peaks_t10, chr, start, end)
  setkey(dt_genes_t10, chr, start, end)
  # 2. Execute the Non-Equi (Spatial) Join with foverlaps
  dt_overlapped <- foverlaps(dt_peaks_t10, dt_genes_t10, nomatch = 0)
  # 3. Count peaks per gene
  dt_count <- dt_overlapped[, .(peak_count = .N), by = .(gene)]
  
  dt_risultato_10_1 <- dt_count
})
time_dt_10_1 <- extract_time(time_dt_10_1_result)
cat("data.table T10.1 Time:", time_dt_10_1, "seconds.\n")

# --- 10.1.3: SQL (sqldf) ---
df_peaks_t10_sql <- df_peaks_t10
df_genes_t10_sql <- df_genes_t10
time_sql_10_1_result <- system.time({
  sql_risultato_10_1 <- sqldf::sqldf("
    SELECT T2.gene, COUNT(T1.chr) AS peak_count
    FROM df_peaks_t10_sql AS T1
    INNER JOIN df_genes_t10_sql AS T2
    ON T1.chr = T2.chr
    WHERE T1.start <= T2.end AND T1.end >= T2.start
    GROUP BY T2.gene
  ")
})
time_sql_10_1 <- extract_time(time_sql_10_1_result)
cat("SQL T10.1 Time:", time_sql_10_1, "seconds.\n")


# Update performance table Task 10.1
performance_comparison_task10 <- rbindlist(list(performance_comparison_task10,
                                                list("T10.1 Spatial Join & Peak Count", time_df_10_1, time_dt_10_1, time_sql_10_1)),
                                           use.names = TRUE, fill = TRUE)

print(head(dt_risultato_10_1, 5))

# ----------------------------------------------------
# TASK 10.2: Calculate Overlap Length and Top 20 Genes
# ----------------------------------------------------
cat("\nRunning Task 10.2: Overlap Length and Top 20 Geni...\n")

# --- 10.2.1: Data.frame / Dplyr ---
time_df_10_2_result <- system.time({
  # Re-execute spatial join
  df_overlapped_base <- df_peaks_t10 %>%
    left_join(df_genes_t10, by = "chr", relationship = "many-to-many") %>%
    filter(start.x <= end.y & end.x >= start.y)
  
  df_overlap_sum <- df_overlapped_base %>%
    mutate(
      # Calculate overlap
      overlap_len = calc_overlap_length(start.x, end.x, start.y, end.y)
    ) %>%
    group_by(gene) %>% # Group by 'gene'
    summarise(total_overlap_bp = sum(overlap_len), .groups = 'drop') %>%
    # Top 20
    arrange(desc(total_overlap_bp)) %>% # Order from largest to smallest
    slice_head(n = 20)
  
  df_risultato_10_2 <- df_overlap_sum
})
time_df_10_2 <- extract_time(time_df_10_2_result)
cat("DataFrame T10.2 Time:", time_df_10_2, "seconds.\n")


# --- 10.2.2: Data.table ---
time_dt_10_2_result <- system.time({
  # Prepare keys (needed for foverlaps)
  setkey(dt_peaks_t10, chr, start, end)
  # Rename columns of the gene file to avoid conflicts
  setnames(dt_genes_t10, old = c("start", "end"), new = c("gene_start", "gene_end"))
  setkey(dt_genes_t10, chr, gene_start, gene_end)

  # Execute spatial join (foverlaps automatically handles intervals)
  dt_overlapped_final <- foverlaps(
    dt_peaks_t10,
    dt_genes_t10,
    by.x = c("chr", "start", "end"),
    by.y = c("chr", "gene_start", "gene_end"),
    nomatch = 0
  )

  # Calculate the overlap length between each peak and the corresponding gene
  dt_overlapped_final[, overlap_len := calc_overlap_length(start, end, gene_start, gene_end)]

  # Aggregate by gene and select the 20 genes with the highest overall overlap
  dt_overlap_sum <- dt_overlapped_final[
    , .(total_overlap_bp = sum(overlap_len)),  # Group by gene and sum the lengths
    by = gene
    ][order(-total_overlap_bp)][1:20] # Order descending and take the first 20 rows
  # Save the final result
  dt_risultato_10_2 <- dt_overlap_sum
})

time_dt_10_2 <- extract_time(time_dt_10_2_result)

cat("data.table T10.2 Time:", time_dt_10_2, "seconds.\n")
print(head(dt_risultato_10_2, 5))

# Update performance table
performance_comparison_task10 <- rbindlist(list(performance_comparison_task10,
                                                list("T10.2 Overlap Length & Top 20", time_df_10_2, time_dt_10_2, NA)), 
                                                use.names = TRUE, fill = TRUE)

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Task 10 results into the main cumulative table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task10), use.names = TRUE, fill = TRUE)

# Saving complete result of Task 10
fwrite(dt_risultato_10_2, file.path(output_dir, "Task10_Top20Overlap_Results.csv"))

cat("\nTask 10 analysis completed. Cumulative performance files saved.\n")




# ----------------------------------------------------
# TASK 11: Map SNPs to Genes
# Goal: Convert variant positions to 1-bp intervals and find overlaps with gene intervals.
# Summarize HIGH impact variants by gene and by sample.
# ----------------------------------------------------

cat("\n--- START TASK 11: SNP-TO-GENE SPATIAL JOIN ---\n")

# Initialize a local performance table for Task 11
performance_comparison_task11 <- data.table(Task = character(),
                                            `Time_DataFrame (s)` = numeric(),
                                            `Time_DataTable (s)` = numeric(),
                                            `Time_SQL (s)` = numeric())

# --- DATA LOADING AND PREPARATION FOR TASK 11 ---
# Load files
dt_variants_t11 <- fread("variants.csv") 
dt_genes_t10 <- fread("gene_annotation.bed.csv") 

# Local copies for Task 11
dt_variants_t11 <- copy(dt_variants_t11)
dt_genes_t11 <- copy(dt_genes_t10)

# Convert variant positions to 1-bp intervals (pos, pos)
# Done in-place on the local copy
dt_variants_t11[, variant_start := pos] # Use 'pos' column
dt_variants_t11[, variant_end := pos]   # Use 'pos' column

# Create data.frame copies for comparisons
df_variants_t11 <- as.data.frame(dt_variants_t11)
df_genes_t11 <- as.data.frame(dt_genes_t10)


# ----------------------------------------------------
# TASK 11.1/11.2: Spatial Join, Filter 'HIGH', and Summarize by Gene/Sample
# ----------------------------------------------------
cat("\nRunning Task 11.1/11.2: Spatial Join, Filter HIGH, and Summarize (3-way comparison)...\n")

# --- 11.1.1: Data.frame / Dplyr ---
time_df_11_1_result <- system.time({
  df_overlapped_snps <- df_variants_t11 %>%
    filter(impact == 'HIGH') %>% 
    left_join(df_genes_t11, by = "chr", relationship = "many-to-many") %>%
    # Spatial join (non-equi) condition: SNP interval must overlap gene interval
    filter(variant_start <= end & variant_end >= start) # Usa 'start' e 'end' dal file geni
  
  # Summary (Task 11.2)
  df_summary_11_2 <- df_overlapped_snps %>%
    group_by(gene, sample_id) %>%
    summarise(high_impact_count = n(), .groups = 'drop')
  
  df_risultato_11_2 <- df_summary_11_2
})
time_df_11_1 <- extract_time(time_df_11_1_result)
cat("DataFrame T11.1/2 Time:", time_df_11_1, "seconds.\n")


# --- 11.1.2: Data.table ---
time_dt_11_1_result <- system.time({
  # 1. Filter HIGH impact first
  dt_variants_high <- dt_variants_t11[impact == 'HIGH']
  # 2. Set keys for spatial join
  setkey(dt_variants_high, chr, variant_start, variant_end)
  setkey(dt_genes_t11, chr, start, end) # Use 'start' and 'end'
  
  ## 3. Execute the Non-Equi (Spatial) Join (GENES[HIGH_VARIANTS, on=...])
  dt_overlapped_snps <- dt_genes_t11[dt_variants_high, 
                                     on = .(chr, start <= variant_end, end >= variant_start),
                                     nomatch = 0]
  # 4. Summary (Task 11.2)
  dt_risultato_11_2 <- dt_overlapped_snps[, .(high_impact_count = .N), by = .(gene, sample_id)]
})
time_dt_11_1 <- extract_time(time_dt_11_1_result)
cat("data.table T11.1/2 Time:", time_dt_11_1, "seconds.\n")

# --- 11.1.3: SQL (sqldf) ---
df_variants_t11_sql <- df_variants_t11
df_genes_t11_sql <- df_genes_t11
time_sql_11_1_result <- system.time({
  sql_risultato_11_2 <- sqldf::sqldf("
    SELECT T2.gene, T1.sample_id, COUNT(T1.chr) AS high_impact_count
    FROM df_variants_t11_sql AS T1
    INNER JOIN df_genes_t11_sql AS T2
    ON T1.chr = T2.chr
    WHERE 
      T1.impact = 'HIGH' AND
      T1.variant_start <= T2.end AND T1.variant_end >= T2.start -- Usa 'start' e 'end'
    GROUP BY T2.gene, T1.sample_id
  ")
})
time_sql_11_1 <- extract_time(time_sql_11_1_result)
cat("SQL T11.1/2 Time:", time_sql_11_1, "seconds.\n")

# Update performance table Task 11.1/2
performance_comparison_task11 <- rbindlist(list(performance_comparison_task11,
                                                list("T11.1/2 Join, Filter & Summary", time_df_11_1, time_dt_11_1, time_sql_11_1)),
                                           use.names = TRUE, fill = TRUE)

print(head(dt_risultato_11_2, 5))

# ----------------------------------------------------
# TASK 11.3: List Genes with HIGH variants in ALL samples
# ----------------------------------------------------
cat("\nRunning Task 11.3: Find Genes with HIGH variants in ALL samples...\n")

# Calculate the total number of unique samples present in the summary data
total_samples_with_high_variants <- dt_risultato_11_2[, uniqueN(sample_id)]
cat("Numero totale di campioni:", total_samples_with_high_variants, "\n")

# Base data: Use the result from the T11.2 summary 

# --- 11.3.1: Data.frame / Dplyr ---
time_df_11_3_result <- system.time({
  df_genes_all_samples <- df_risultato_11_2 %>%
    group_by(gene) %>%
    summarise(n_samples_hit = n_distinct(sample_id), .groups = 'drop') %>%
    filter(n_samples_hit == total_samples_with_high_variants)
})
time_df_11_3 <- extract_time(time_df_11_3_result)
cat("DataFrame T11.3 Time:", time_df_11_3, "seconds.\n")


# --- 11.3.2: Data.table ---
time_dt_11_3_result <- system.time({
  dt_genes_all_samples <- dt_risultato_11_2[,
                                            .(n_samples_hit = uniqueN(sample_id)), 
                                            by = gene
  ][
    n_samples_hit == total_samples_with_high_variants
  ]
})
time_dt_11_3 <- extract_time(time_dt_11_3_result)
cat("data.table T11.3 Time:", time_dt_11_3, "seconds.\n")


# --- 11.3.3: SQL (sqldf) ---
sql_risultato_11_2_sql <- dt_risultato_11_2 # Usa il risultato DT come input
time_sql_11_3_result <- system.time({
  sql_genes_all_samples <- sqldf::sqldf(
    paste("SELECT gene, COUNT(DISTINCT sample_id) AS n_samples_hit
           FROM sql_risultato_11_2_sql
           GROUP BY gene
           HAVING n_samples_hit =", total_samples_with_high_variants)
  )
})
time_sql_11_3 <- extract_time(time_sql_11_3_result)
cat("SQL T11.3 Time:", time_sql_11_3, "seconds.\n")

# Update performance table Task 11.3
performance_comparison_task11 <- rbindlist(list(performance_comparison_task11,
                                                list("T11.3 Genes across all samples", time_df_11_3, time_dt_11_3, time_sql_11_3)),
                                           use.names = TRUE, fill = TRUE)

print(head(dt_genes_all_samples, 5))


# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Task 11 results into the main cumulative table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task11), use.names = TRUE, fill = TRUE)

# Saving complete result Task 11.2 (Gene/Sample Summary)
fwrite(dt_risultato_11_2, file.path(output_dir, "Task11_HighImpactSummary_Results.csv"))

# Saving complete result Task 11.3 (Genes in all samples)
fwrite(dt_genes_all_samples, file.path(output_dir, "Task11_GenesAcrossAllSamples_Results.csv"))

cat("\nTask 11 analysis completed. Cumulative performance files saved.\n")



# ----------------------------------------------------
# TASK 12: Combine Cohorts Safely
# Goal: rbindlist() and calculate means per cohort.
# ----------------------------------------------------
cat("\n--- ISTART TASK 12: COMBINE COHORTS ---\n")

# Initialize a local performance table for Task 12
performance_comparison_task12 <- data.table(Task = character(),
                                            `Time_DataFrame (s)` = numeric(), # Dplyr
                                            `Time_DataTable (s)` = numeric(), # data.table
                                            `Time_SQL (s)` = numeric())       # sqldf

# --- DATA LOADING FOR TASK 12 ---
cat("Loading Task 12 data (Cohort A and B)...\n")
# Ensure files are available (loaded globally or here)
  dt_cohortA <- fread("cohortA_samples.csv")
  df_cohortA <- as.data.frame(dt_cohortA)

  dt_cohortB <- fread("cohortB_samples.csv")
  df_cohortB <- as.data.frame(dt_cohortB)

  # Note: Global dt_counts is used later for variance calculation

# ----------------------------------------------------
# TASK 12.1: Combine Cohorts (rbindlist) and Order
# ----------------------------------------------------
cat("\nRunning Task 12.1: Combine Cohorts (3-way comparison)...\n")

# --- 12.1.1: Data.frame / Dplyr (bind_rows) ---
time_df_12_1_result <- system.time({
  # Tidyverse function to stack tables
  df_combined <- bind_rows(df_cohortA, df_cohortB) %>% 
    arrange(cohort, condition, sample_id) # Order the table
  
  df_risultato_12_1_df <- df_combined # Save Dplyr result
})
time_df_12_1 <- extract_time(time_df_12_1_result)
cat("DataFrame T12.1 Time (bind_rows):", time_df_12_1, "seconds.\n")


# --- 12.1.2: Data.table (rbindlist) ---
time_dt_12_1_result <- system.time({
  dt_combined <- rbindlist(list(dt_cohortA, dt_cohortB), use.names = TRUE, fill = TRUE)
  setorder(dt_combined, cohort, condition, sample_id) # Order using data.table's setorder
  
  dt_risultato_12_1_dt <- dt_combined # Save data.table result
})
time_dt_12_1 <- extract_time(time_dt_12_1_result)
cat("data.table T12.1 Time (rbindlist):", time_dt_12_1, "seconds.\n")


# --- 12.1.3: SQL (sqldf) ---
time_sql_12_1_result <- system.time({
  # SQL uses UNION ALL (equivalent of rbind)
  # fill=TRUE used by UNION ALL 
  sql_risultato_12_1_sql <- sqldf::sqldf("
        SELECT * FROM df_cohortA
        UNION ALL
        SELECT * FROM df_cohortB
        ORDER BY cohort, condition, sample_id
    ")
})
time_sql_12_1 <- extract_time(time_sql_12_1_result)
cat("SQL T12.1 Time (UNION ALL):", time_sql_12_1, "seconds.\n")


# Update performance table Task 12.1
performance_comparison_task12 <- rbindlist(list(performance_comparison_task12,
                                                list("T12.1 Combine & Sort Cohorts", time_df_12_1, time_dt_12_1, time_sql_12_1)),
                                           use.names = TRUE, fill = TRUE)

print(head(dt_risultato_12_1_dt, 5)) 


# ----------------------------------------------------
# TASK 12.2: Join and Calculate Means Per Cohort/Condition
# (Top 100 Most Variable Genes)
# ----------------------------------------------------
cat("\nRunning Task 12.2: Join & Calculate Cohort Means (3-way comparison)...\n")

# --- 12.2.1: Data.frame / Dplyr ---
time_df_12_2_result <- system.time({
  # 1. Find Top 100 most variable genes
  top100_genes_df <- df_counts %>%
    group_by(gene) %>%
    summarise(variance = var(count, na.rm = TRUE)) %>%
    arrange(desc(variance)) %>%
    slice_head(n = 100) %>%
    pull(gene) 
  
  # 2. Join and Calculate Means
  df_risultato_12_2 <- df_counts %>%
    filter(gene %in% top100_genes_df) %>% 
    left_join(df_risultato_12_1_df, by = "sample_id") %>% # Use combined data from 12.1.1
    filter(!is.na(cohort)) %>% 
    group_by(cohort, condition, gene) %>%
    summarise(mean_count = mean(count, na.rm = TRUE), .groups = 'drop') # na.rm = Removes any NA values
})
time_df_12_2 <- extract_time(time_df_12_2_result)
cat("DataFrame T12.2 Time:", time_df_12_2, "seconds.\n")


# --- 12.2.2: Data.table ---
time_dt_12_2_result <- system.time({
  # 1. Find Top 100 most variable genes
  top100_genes_dt <- dt_counts[, .(variance = var(count, na.rm = TRUE)), by = gene][order(-variance)][1:100, gene]
  
  # 2. Join and Calculate Means (in a chain)
  dt_risultato_12_2 <- dt_counts[gene %in% top100_genes_dt][ # Filter for Top 100
    dt_risultato_12_1_dt, on = "sample_id", nomatch = 0      # Join with combined data from 12.1.2
  ][,
    .(mean_count = mean(count, na.rm = TRUE)),# Calculate mean
    by = .(cohort, condition, gene)           # Group
  ]
})
time_dt_12_2 <- extract_time(time_dt_12_2_result)
cat("data.table T12.2 Time:", time_dt_12_2, "seconds.\n")


# Update performance table Task 12.2
performance_comparison_task12 <- rbindlist(list(performance_comparison_task12,
                                                list("T12.2 Join & Cohort Means (Top 100)", time_df_12_2, time_dt_12_2)),
                                           use.names = TRUE, fill = TRUE)

print(head(dt_risultato_12_2, 5))

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Task 12 results into the main cumulative table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_task12), use.names = TRUE, fill = TRUE)

# Saving complete result of Task 12
fwrite(dt_risultato_12_2, file.path(output_dir, "Task12_CohortMeans_Results.csv"))

cat("\nTask 12 analysis completed. Cumulative performance files saved.\n")
cat("\n*** TOTAL ANALYSIS COMPLETED UP TO TASK 12 ***\n")



# ----------------------------------------------------
# TASK FINAL REVISION: Cell Type Association to Clusters (N vs T)
# Goal: Associate cell type to integration clusters, calculating normalized percentages.
# ----------------------------------------------------

cat("\n--- START FINAL REVISION: SINGLE-CELL ANALYSIS ---\n")

# Initialize a local performance table for the Final Revision
performance_comparison_final <- data.table(Task = character(),
                                           `Time_DataFrame (s)` = numeric(),
                                           `Time_DataTable (s)` = numeric(),
                                           `Time_SQL (s)` = numeric()) 

# --- CARICAMENTO DATI LOCALI PER FINAL REVISION ---
cat("Loading Final Revision data...\n")
dt_clusters <- fread("annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv")
dt_celltypes <- fread("nt_combined_clustering.output.csv")

# --- FIX: CLEANUP AND STANDARDIZATION OF THE JOIN KEY 'cell' ---
# Remove leading/trailing spaces
dt_clusters[, cell := trimws(as.character(cell))]
dt_celltypes[, cell := trimws(as.character(cell))]

# CRUCIAL FIX (Regex): Remove _X_, _Y_, _T_, etc. to match formats
# Searches for an underscore, followed by ONE uppercase letter, followed by an underscore,
# and replaces it with "" (nothing). This is done because the two files present the data differently.
dt_clusters[, cell := gsub("_[A-Z]_", "", cell)]

# Create data.frame copies (already cleaned)
df_clusters <- as.data.frame(dt_clusters)
df_celltypes <- as.data.frame(dt_celltypes)

# --- CHECK, ensuring our "cleaning" was successful ---
common_cells <- length(intersect(dt_clusters$cell, dt_celltypes$cell))
cat("DEBUG: Found", common_cells, "common cells between the two files.\n")
if (common_cells == 0) {
  cat("WARNING: No common cells found. Joins will fail!\n")
}

# ----------------------------------------------------
# FR-TASK 1: Combine Cell Types and Clusters
# ----------------------------------------------------
cat("\nRunning FR-Task 1: Combine Cell Types and Clusters...\n")

# --- 1.1: Data.frame / Dplyr (inner_join) ---
time_df_fr_1_result <- system.time({
  # Use inner_join to keep only cells present in both sets
  df_master_table <- df_clusters %>%
    inner_join(df_celltypes, by = "cell")
})
time_df_fr_1 <- extract_time(time_df_fr_1_result)
cat("DataFrame FR-T1 Time (inner_join):", time_df_fr_1, "seconds.\n")


# --- 1.2: Data.table (Join) ---
time_dt_fr_1_result <- system.time({
  setkey(dt_clusters, cell) # Order both tables by cell
  setkey(dt_celltypes, cell)
  # Use an inner join (nomatch = 0 means ignore cells that do not match)
  dt_master_table <- dt_clusters[dt_celltypes, nomatch = 0]
})
time_dt_fr_1 <- extract_time(time_dt_fr_1_result)
cat("data.table FR-T1 Time (Keyed Join):", time_dt_fr_1, "seconds.\n")


# --- 1.3: SQL (sqldf) ---
time_sql_fr_1_result <- system.time({
  # Use the cleaned dataframes
  sql_master_table <- sqldf::sqldf("
    SELECT T1.*, T2.cell_type, T2.sample_type
    FROM df_clusters AS T1
    INNER JOIN df_celltypes AS T2 ON T1.cell = T2.cell
  ")
})
time_sql_fr_1 <- extract_time(time_sql_fr_1_result)
cat("SQL FR-T1 Time (INNER JOIN):", time_sql_fr_1, "seconds.\n")

# Update performance table FR-Task 1
performance_comparison_final <- rbindlist(list(performance_comparison_final,
                                               list("FR-T1 Combine Tables", time_df_fr_1, time_dt_fr_1, time_sql_fr_1)),
                                          use.names = TRUE, fill = TRUE)

print(head(dt_master_table, 5))


# ----------------------------------------------------
# FR-TASK 2: Count Cell Types per Cluster (Total)
# ----------------------------------------------------
cat("\nRunning FR-Task 2: Count Cell Types per Cluster (Total)...\n")

# --- 2.1: Data.frame / Dplyr (count) ---
time_df_fr_2_result <- system.time({
  df_counts_per_cluster <- df_master_table %>%
    count(integration_cluster, cell_type, name = "count")
})
time_df_fr_2 <- extract_time(time_df_fr_2_result)
cat("DataFrame FR-T2 Time (count):", time_df_fr_2, "seconds.\n")


# --- 2.2: Data.table (.N by) ---
time_dt_fr_2_result <- system.time({
  # Count rows within integration cluster and cell type
  dt_counts_per_cluster <- dt_master_table[, .N, by = .(integration_cluster, cell_type)]
  # Rename column N to count
  setnames(dt_counts_per_cluster, "N", "count") # Rename for consistency
})
time_dt_fr_2 <- extract_time(time_dt_fr_2_result)
cat("data.table FR-T2 Time (.N by):", time_dt_fr_2, "seconds.\n")


# --- 2.3: SQL (sqldf) ---
time_sql_fr_2_result <- system.time({
  sql_counts_per_cluster <- sqldf::sqldf("
    SELECT integration_cluster, cell_type, COUNT(*) AS count
    FROM sql_master_table
    GROUP BY integration_cluster, cell_type
  ")
})
time_sql_fr_2 <- extract_time(time_sql_fr_2_result)
cat("SQL FR-T2 Time (GROUP BY):", time_sql_fr_2, "seconds.\n")

# Update performance table FR-Task 2
performance_comparison_final <- rbindlist(list(performance_comparison_final,
                                               list("FR-T2 Count per Cluster", time_df_fr_2, time_dt_fr_2, time_sql_fr_2)),
                                          use.names = TRUE, fill = TRUE)

print(head(dt_counts_per_cluster, 5))


# ----------------------------------------------------
# FR-TASK 3 & 5: Count and Normalize % by Cluster and Tissue
# ----------------------------------------------------
cat("\nRunning FR-Task 3 & 5: Count & Normalize % by Cluster and Tissue...\n")

# --- 3&5.1: Data.frame / Dplyr ---
time_df_fr_3_5_result <- system.time({
  # Task 3: Counting
  df_counts_tissue <- df_master_table %>%
    count(integration_cluster, cell_type, sample_type, name = "count")
  
  # Task 5: % Normalization
  df_normalized_tissue <- df_counts_tissue %>%
    group_by(integration_cluster, sample_type) %>%
    mutate(
      # Calculate total per group (total number of cells for each cluster and condition)
      total_in_group = sum(count),
      percentage = (count / total_in_group) * 100
    ) %>%
    ungroup()
})
time_df_fr_3_5 <- extract_time(time_df_fr_3_5_result)
cat("DataFrame FR-T3/5 Time (count + mutate):", time_df_fr_3_5, "seconds.\n")


# --- 3&5.2: Data.table ---
time_dt_fr_3_5_result <- system.time({
  # Task 3: Counting
  dt_counts_tissue <- dt_master_table[, .N, by = .(integration_cluster, cell_type, sample_type)]
  setnames(dt_counts_tissue, "N", "count") # Rename
  
  # Task 5: % Normalization (in-place)
  dt_counts_tissue[, total_in_group := sum(count), by = .(integration_cluster, sample_type)]
  dt_counts_tissue[, percentage := (count / total_in_group) * 100]
  
  dt_normalized_tissue <- dt_counts_tissue
})
time_dt_fr_3_5 <- extract_time(time_dt_fr_3_5_result)
cat("data.table FR-T3/5 Time (.N by + :=):", time_dt_fr_3_5, "seconds.\n")


# --- 3&5.3: SQL (sqldf) ---
time_sql_fr_3_5_result <- system.time({
  sql_normalized_tissue <- sqldf::sqldf("
    WITH Counts AS (
      SELECT 
        integration_cluster, 
        cell_type, 
        sample_type, 
        COUNT(*) AS count
      FROM sql_master_table
      GROUP BY integration_cluster, cell_type, sample_type
    )
    SELECT 
      *,
      (CAST(count AS REAL) / SUM(count) OVER (PARTITION BY integration_cluster, sample_type)) * 100 AS percentage
    FROM Counts
  ")
})
time_sql_fr_3_5 <- extract_time(time_sql_fr_3_5_result)
cat("SQL FR-T3/5 Time (Window Function):", time_sql_fr_3_5, "seconds.\n")

# Update performance table FR-Task 3/5
performance_comparison_final <- rbindlist(list(performance_comparison_final,
                                               list("FR-T3/5 Count & Normalize", time_df_fr_3_5, time_dt_fr_3_5, time_sql_fr_3_5)),
                                          use.names = TRUE, fill = TRUE)

print(head(dt_normalized_tissue, 5))


# ----------------------------------------------------
# FR-TASK 6: Plot Generation (Cell Type Distribution)
# ----------------------------------------------------

# --- 1. Plot Generation (Assignment to variable) ---
dt_plot_data <- dt_normalized_tissue 

plot_final_revision <- ggplot(dt_plot_data, 
                              aes(x = integration_cluster, 
                                  y = percentage, 
                                  fill = cell_type)) +
  
  geom_bar(stat = "identity", position = "stack") +
  
  facet_wrap(~ sample_type) +
  
  labs(title = "Cell Type Distribution within Integration Clusters, by Tissue Type",
       y = "Normalized Percentage (%)",
       x = "Integration Cluster",
       fill = "Cell Type") +
  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14)) +
  
  scale_y_continuous(labels = scales::percent_format(scale = 1))

# --- 2. Saving the Plot (PNG file) ---
# This will save the file to the Docker working directory output folder
ggsave(filename = file.path(output_dir, "FinalTask_CellTypeDistribution.png"), 
       plot = plot_final_revision, 
       width = 10, 
       height = 6, 
       units = "in", 
       dpi = 300)

# Print the plot (useful for immediate visualization in R/RStudio)
print(plot_final_revision)

# ----------------------------------------------------
# FINAL RESULTS SAVING AND CUMULATIVE UPDATE
# ----------------------------------------------------
# Merge Final Revision results into the main cumulative table
performance_comparison <- rbindlist(list(performance_comparison, performance_comparison_final), use.names = TRUE,fill = TRUE)

print(performance_comparison)

# Saving complete Final Revision results
if (exists("dt_master_table") && nrow(dt_master_table) > 0) {
  fwrite(dt_master_table, file.path(output_dir, "FinalTask_MasterTable.csv"))
  fwrite(dt_counts_per_cluster, file.path(output_dir, "FinalTask_CountsPerCluster.csv"))
  fwrite(dt_normalized_tissue, file.path(output_dir, "FinalTask_NormalizedCounts_Tissue.csv"))
} else {
  cat("FINAL REVISION RESULTS EMPTY - NO FILES SAVED (JOIN FAILED).\n")
}

cat("\nFinal Revision analysis completed. Cumulative performance files saved.\n")