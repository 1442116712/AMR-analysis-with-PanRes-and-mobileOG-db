library(readxl)
library(openxlsx)
library(dplyr)

# Overlap MGE de-duplication ####
df <- read_excel('C:/Users/CFL/Desktop/UK/PhD/MGE.xlsx', sheet = "MGE")

# Calculate coverage and filter (Threshold: 90%)
# Updated variable names to lowercase: query_sequence_length, subject_sequence_length
df$coverage <- (df$end_of_alignment_in_query - df$start_of_alignment_in_query + 1) / df$subject_sequence_length 
df <- df[df$coverage >= 0.9, ]
df$coverage <- ifelse(df$coverage > 1, 1, df$coverage)

# Calculate the midpoint location for MGEs
# Logic: (start + stop) / 2
df$midpoint_MGE <- df$start + (df$stop - df$start) / 2

# De-duplication Function
remove_overlapping_MGEs <- function(df) {
  # Pre-calculate original indices to map back from sub-groups to the main data frame
  df$orig_idx <- seq_len(nrow(df))
  # Initialize a logical vector: TRUE means keep, FALSE means remove
  keep <- rep(TRUE, nrow(df))
  
  # Optimization: Identify groups (specific_contig + strand) that have at least 2 entries
  groups <- df %>% 
    group_by(specific_contig, strand) %>% 
    summarise(n = n(), .groups = 'drop') %>% 
    filter(n >= 2)
  
  # Iterate through groups with potential overlaps
  for (g in seq_len(nrow(groups))) {
    # Extract global indices belonging to the current group
    idx <- which(df$specific_contig == groups$specific_contig[g] & 
                   df$strand == groups$strand[g])
    
    # Create a temporary subset for the current group
    df_sub <- df[idx, ]
    n_sub <- nrow(df_sub)
    
    # Pairwise comparison within the group
    for (i in 1:(n_sub - 1)) {
      for (j in (i + 1):n_sub) {
        
        # Get the global row indices for the two MGEs being compared
        idx_i <- idx[i]
        idx_j <- idx[j]
        
        # CRITICAL: If either MGE has already been marked for removal, skip this pair
        if (!keep[idx_i] | !keep[idx_j]) next
        
        # Calculate the physical distance between midpoints
        dist <- abs(df_sub$midpoint_MGE[i] - df_sub$midpoint_MGE[j])
        
        # Calculate the allowed distance (Threshold for overlap: 2 bp buffer)
        allowed <- (df_sub$query_sequence_length[i] + df_sub$query_sequence_length[j]) / 2 - 2
        
        # If an overlap is detected (distance < allowed)
        if (!is.na(dist) && dist < allowed) {
          # Competition: Retain the hit with the higher bitscore (lowercase)
          if (df_sub$bitscore[i] < df_sub$bitscore[j]) {
            keep[idx_i] <- FALSE # Mark i for removal
          } else {
            keep[idx_j] <- FALSE # Mark j for removal
          }
        }
      }
    }
  }
  
  # Return the filtered data frame, removing the helper index column
  return(df[keep, -which(names(df) == "orig_idx")])
}

# Execute Filtering
df_final <- remove_overlapping_MGEs(df)

# Final cleanup: Remove any entirely empty or NA rows
df_final <- df_final %>% filter_all(any_vars(. !="" & !is.na(.)))

write.xlsx(df_final, 'C:/Users/CFL/Desktop/UK/PhD/output_file_MGE.xlsx')





# ARG - MGE Co-localization Analysis ####
library(readxl)
library(openxlsx)
library(dplyr)

# 1. Load data
df_arg <- read_excel('C:/Users/CFL/Desktop/UK/PhD/MGE.xlsx', sheet = "ARG")
df_mge <- read_excel('C:/Users/CFL/Desktop/UK/PhD/MGE.xlsx', sheet = "MGE_unique")

# 2. Define the linkage function with specific MGE column selection
link_MGE_ARG <- function(df_mge, df_arg) {
  
  # Step 1: Pre-select specific columns from the MGE table
  # This ensures we only carry relevant MGE metadata
  mge_subset <- df_mge %>%
    select(
      specific_contig, pident, bitscore, subject_sequence_length, 
      e_value, query_sequence_length, mobileOG_ID, gene_name, 
      best_hit_accession_ID, major_mobileOG_category, minor_mobileOG_category, 
      source_database, start, stop, strand, gc_content, coverage, midpoint_MGE
    )
  
  # Step 2: Prepare the ARG data for linkage (make sure not replicate column name)
  # We extract identity info and taxonomy to attach to the MGEs
  arg_subset <- df_arg %>%
    select(
      Contig, class, fa_name, midpoint_ARG, 
      dataset, Phylum, Class, Order, Family, Genera, Species, Genome
    )
  
  # Step 3: Join tables and calculate distance
  # Using many-to-many relationship to capture all ARG neighbors for each MGE
  linked_df <- mge_subset %>%
    inner_join(arg_subset, by = c("specific_contig" = "Contig"), relationship = "many-to-many") %>%
    # Calculate physical gap between midpoints
    mutate(gap_distance = abs(midpoint_MGE - midpoint_ARG)) %>%
    # Filter for proximity (within 10kb)
    filter(gap_distance <= 10000)
  
  return(linked_df)
}

# 3. Execute the linkage
df_final_linked <- link_MGE_ARG(df_mge, df_arg)

# 4. Export the results
write.xlsx(df_final_linked, 'C:/Users/CFL/Desktop/UK/PhD/output_file_MGE_ARG_links.xlsx')
