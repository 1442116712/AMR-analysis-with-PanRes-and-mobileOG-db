library(readxl)
library(openxlsx)
library(dplyr)

#Necessary information in excel
#MGE - Specific_Contig,Gene_Name,ORF_End,ORF_Start,End/Start_of_Alignment_in_Subject/Query,Query_Sequence_Length,Subject_Sequence_Length,Bitscore
#ARG - Contig,class,fa_name,qend,qstart,send,start

# Overlap gene de-duplication ####
df <- read_excel('C:/Users/KEN/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "MGE")

# filter ARG coverage < 90%
df$coverage <- df$Query_Sequence_Length / df$Subject_Sequence_Length 
df <- df[df$coverage >= 0.9, ]
df$coverage <- ifelse(df$coverage > 1, 1, df$coverage)

# determine strand and midpoint
df$strand <- ifelse(df$End_of_Alignment_in_Query - df$Start_of_Alignment_in_Query > 0, "F", "R")
df$q_direction <- NULL
df$s_direction <- NULL
df$center_MGE <- df$ORF_Start + (df$ORF_End - df$ORF_Start) / 2

# De-duplication
remove_overlapping_MGEs <- function(df) {
  
  # Initialise the logic vector to indicate retention of all genes
  keep <- rep(TRUE, nrow(df))
  
  # Group comparable vectors based on the same contig + class + strand
  for (contig in unique(df$Specific_Contig)) {
    for (gene in unique(df$Gene_Name)) {
      for (strand in unique(df$strand)) {  
        df_contig_gene_strand <- df[df$Specific_Contig == contig & df$Gene_Name == gene & df$strand == strand, ]
        
        # Skip groups with only one vector
        if (nrow(df_contig_gene_strand) < 2) {
          next
        }
        
        # The vectors of the group compare each other in turn 
        for (i in 1:(nrow(df_contig_gene_strand) - 1)) {
          for (j in (i + 1):nrow(df_contig_gene_strand)) {
            # Get the center and length of all genes
            center_i <- df_contig_gene_strand$center_MGE[i]
            center_j <- df_contig_gene_strand$center_MGE[j]
            length_i <- df_contig_gene_strand$Query_Sequence_Length[i]
            length_j <- df_contig_gene_strand$Query_Sequence_Length[j]
            
            # Check for NA values and skip 
            if (is.na(center_i) || is.na(center_j) || is.na(length_i) || is.na(length_j) || 
                is.na(df_contig_gene_strand$Bitscore[i]) || is.na(df_contig_gene_strand$Bitscore[j])) {
              next
            }
            
            # Calculate the distance from the centers
            actual_distance <- abs(center_i - center_j)
            # Calculate the allowed distance from length (less than 4 units are allowed overlapping, can change to 12 if using nucleotide data)
            allowed_distance <- (length_i + length_j) / 2 - 12
            
            # If the overlap area is greater than or equal to 12 units, it is considered to be an overlap gene
            if (!is.na(actual_distance) && actual_distance < allowed_distance) {
              # label and delete overlap genes with smaller bitscore, the labelled objects will not be compared next turn
              if (df_contig_gene_strand$Bitscore[i] < df_contig_gene_strand$Bitscore[j]) {
                keep[which(df$Specific_Contig == contig & df$Gene_Name == gene & df$strand == strand)[i]] <- FALSE
              } else {
                keep[which(df$Specific_Contig == contig & df$Gene_Name == gene & df$strand == strand)[j]] <- FALSE
              }
            }
          }
        }
      }
    }
  }
  
  # Return of retained genes
  return(df[keep, ])
}

# Output the filter result
filtered_df <- remove_overlapping_MGEs(df)
df_final <- filtered_df %>% filter_all(any_vars(. !="" & !is.na(.)))
# Backup the data
df_backup <- df_final

# ARG Grouping ####
df1 <- read_excel('C:/Users/KEN/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "ARG")

# Define the grouping function
group_ARGs <- function(df1) {
  df1$Group <- NA  # 创建 Group 列
  
  # Grouped by Contig
  for (contig in unique(df1$Contig)) {
    subset_df <- df1[df1$Contig == contig, ]  # 选出当前 Contig 的数据
    
    # Grouped by class
    for (class in unique(subset_df$class)) {
      class_rows <- subset_df[subset_df$class == class, ]  # 选出当前 class 的数据
      
      # Assign group numbers by row
      if (nrow(class_rows) > 0) {
        df1$Group[df1$Contig == contig & df1$class == class] <- 1:nrow(class_rows)
      }
    }
  }
  
  return(df1)
}

df1_grouped <- group_ARGs(df1)

# ARG and MGE were matched, and results with a midpoint distance of less than 10 kb were filtered and populated with the corresponding ARG information ####
MGE_compilation <- function(df_final, df1_grouped) {
  # Retrieve all unique class values
  unique_classes <- unique(df1_grouped$class)
  
  # In df_final, create a column for each class, appending the column name with `_gap`
  df_final[, paste0(unique_classes, "_gap")] <- NA
  
  # Create a fa_name column for each class
  for (class_name in unique_classes) {
    fa_name_col <- paste0(class_name, "_fa_name")
    df_final[[fa_name_col]] <- NA
  }
  
  for (group in unique(df1_grouped$Group)) {
    df_grouped_1 <- df1_grouped %>% filter(Group == group)
    
    for (class_name in unique(df_grouped_1$class)) {
      df_class <- df_grouped_1 %>% filter(class == class_name)
      
      # Match using the centre_ARG column, and incidentally match fa_name
      matched_data <- df_class %>% select(Contig, center_ARG, fa_name)
      
      # Perform a left join, populating `df_final` with `center_ARG` and `fa_name`
      df_final <- df_final %>%
        left_join(matched_data, by = c("Specific_Contig" = "Contig"))
      
      if ("center_ARG" %in% colnames(df_final)) {
        # Fill in centre_ARG and fa_name
        df_final <- df_final %>%
          mutate(
            !!paste0(class_name, "_gap") := ifelse(!is.na(center_ARG),
                                                   ifelse(is.na(!!sym(paste0(class_name, "_gap"))), center_ARG, paste(!!sym(paste0(class_name, "_gap")), center_ARG, sep = ";")),
                                                   !!sym(paste0(class_name, "_gap"))),
            !!paste0(class_name, "_fa_name") := ifelse(!is.na(fa_name),
                                                       ifelse(is.na(!!sym(paste0(class_name, "_fa_name"))), fa_name, paste(!!sym(paste0(class_name, "_fa_name")), fa_name, sep = ";")),
                                                       !!sym(paste0(class_name, "_fa_name")))
          ) %>%
          select(-center_ARG, -fa_name) 
      }
    }
  }
  
  # Calculate the absolute value of the difference between centre_ARG and centre_MGE, and update the centre_ARG column
  df_final <- df_final %>%
    mutate(across(all_of(paste0(unique_classes, "_gap")), 
                  ~ sapply(seq_along(.), function(i) {
                    val <- .[i]
                    if (!is.na(val)) {  
                      # Ensure that val is of character type
                      val <- as.character(val)
                      
                      # Check for multiple values separated by semicolons
                      values <- unlist(strsplit(val, ";"))
                      
                      # Retrieve the corresponding centre_MGE value for this row
                      center_MGE_value <- center_MGE[i]
                      
                      # Calculate the absolute value of the difference between each value and the corresponding centre_MGE for that row
                      abs_diffs <- sapply(values, function(v) abs(as.numeric(v) - center_MGE_value))
                      
                      # Only retain values where the difference is less than or equal to 10,000
                      filtered_values <- values[abs_diffs <= 10000]
                      filtered_diffs <- abs_diffs[abs_diffs <= 10000]
                      
                      # If there are multiple differences, merge them into a single string
                      if (length(filtered_diffs) == 1) {
                        return(filtered_diffs[1])  # If there is only one difference, return directly
                      } else if (length(filtered_diffs) > 1) {
                        return(paste(filtered_diffs, collapse = ";"))  # If there are multiple differences, combine them into a single string
                      } else {
                        return(NA)  # If no qualifying value is found, return NA
                      }
                    } else {
                      return(NA) 
                    }
                  })))
  
  # Filter out all rows where the difference exceeds 10,000 corresponding to fa_name
  df_final <- df_final %>%
    filter(if_any(all_of(paste0(unique_classes, "_gap")), ~ !is.na(.))) %>%
    mutate(across(ends_with("_fa_name"), ~ ifelse(is.na(.), "NA", .))) %>%
    mutate(across(ends_with("_gap"), ~ ifelse(is.na(.), "NA", .))) %>% 
    select_if(~ !all(. == "NA")) 
  
  df_final <- df_final %>%
    left_join(
      df1 %>%
        select(Contig, Kingdom, Phylum, Class, Order, Family, Genera, Species, Genome) %>%
        distinct(Contig, .keep_all = TRUE), 
      by = c("Specific_Contig" = "Contig")
    )
  return(df_final)
}

df_final <- MGE_compilation(df_final, df1_grouped)

# Print results
write.xlsx(df_final, 'C:/Users/KEN/Desktop/UK/PhD/output_file_MGE.xlsx')

