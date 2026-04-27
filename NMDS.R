# ============================================
# NMDS Analysis Script
# ============================================

# Load required libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================
# FILE READING FUNCTIONS
# ============================================
read_uploaded_file <- function(filepath, filename) {
  ext <- tolower(tools::file_ext(filename))
  
  result <- tryCatch({
    if (ext %in% c("csv")) {
      df <- tryCatch(read.csv(filepath, check.names = FALSE, stringsAsFactors = FALSE), 
                     error = function(e) NULL)
      if (is.null(df)) {
        df <- read.csv2(filepath, check.names = FALSE, stringsAsFactors = FALSE)
      }
    } else if (ext %in% c("txt", "tsv")) {
      df <- tryCatch(read.delim(filepath, check.names = FALSE, stringsAsFactors = FALSE), 
                     error = function(e) NULL)
      if (is.null(df)) {
        df <- tryCatch(read.csv(filepath, check.names = FALSE, stringsAsFactors = FALSE), 
                       error = function(e) NULL)
      }
      if (is.null(df)) {
        df <- read.csv2(filepath, check.names = FALSE, stringsAsFactors = FALSE)
      }
    } else {
      stop(paste("Unsupported file type:", ext))
    }
    return(df)
  }, error = function(e) {
    stop(paste("Error reading file:", e$message))
  })
  return(result)
}

# ============================================
# DATA TRANSFORMATION FUNCTION
# ============================================
apply_transformation <- function(otu_mat, method) {
  otu_mat[is.na(otu_mat)] <- 0
  otu_mat[is.infinite(otu_mat)] <- 0
  
  result <- switch(method,
                   "none" = otu_mat,
                   "wisconsin" = wisconsin(otu_mat),
                   "hellinger" = decostand(otu_mat, "hellinger"),
                   "pa" = decostand(otu_mat, "pa"),
                   "chi.square" = decostand(otu_mat, "chi.square"),
                   "log" = log1p(otu_mat),
                   "relabund" = decostand(otu_mat, "total"),
                   "sqrt" = sqrt(otu_mat),
                   "standardize" = decostand(otu_mat, "standardize"),
                   otu_mat)
  
  result[is.na(result)] <- 0
  result[is.infinite(result)] <- 0
  return(result)
}

# ============================================
# DATA FILTERING FUNCTION
# ============================================
filter_otu_data <- function(otu_mat, 
                            filter_samples = FALSE, min_reads = 100,
                            filter_otus = TRUE, min_occurrence = 2,
                            filter_low_abundance = FALSE, min_abundance = 10) {
  
  filter_messages <- c()
  
  # Filter samples by read depth
  if(filter_samples && min_reads > 0) {
    sample_reads <- rowSums(otu_mat)
    before <- nrow(otu_mat)
    otu_mat <- otu_mat[sample_reads >= min_reads, , drop = FALSE]
    after <- nrow(otu_mat)
    if(after < before) {
      filter_messages <- c(filter_messages, paste("Removed", before - after, "samples with <", min_reads, "reads"))
    }
  }
  
  if(nrow(otu_mat) < 3) {
    stop("Too few samples after filtering - need at least 3")
  }
  
  # Filter rare OTUs (by occurrence)
  if(filter_otus && min_occurrence > 0) {
    otu_occurrence <- colSums(otu_mat > 0)
    before <- ncol(otu_mat)
    otu_mat <- otu_mat[, otu_occurrence >= min_occurrence, drop = FALSE]
    after <- ncol(otu_mat)
    if(after < before) {
      filter_messages <- c(filter_messages, paste("Removed", before - after, "rare OTUs (present in <", min_occurrence, "samples)"))
    }
  }
  
  # Filter low abundance OTUs
  if(filter_low_abundance && min_abundance > 0) {
    otu_abundance <- colSums(otu_mat)
    before <- ncol(otu_mat)
    otu_mat <- otu_mat[, otu_abundance >= min_abundance, drop = FALSE]
    after <- ncol(otu_mat)
    if(after < before) {
      filter_messages <- c(filter_messages, paste("Removed", before - after, "low-abundance OTUs (<", min_abundance, "total reads)"))
    }
  }
  
  if(ncol(otu_mat) < 2) {
    stop("Too few OTUs after filtering - need at least 2")
  }
  
  # Remove zero-sum samples
  sample_sums <- rowSums(otu_mat)
  otu_mat <- otu_mat[sample_sums > 0, , drop = FALSE]
  
  if(nrow(otu_mat) < 3) {
    stop("Too few samples with non-zero reads after filtering")
  }
  
  # Print filtering summary
  if(length(filter_messages) > 0) {
    cat("\nFiltering Summary:\n")
    cat(paste(filter_messages, collapse = "\n"))
    cat(sprintf("\nFinal: %d samples, %d OTUs\n", nrow(otu_mat), ncol(otu_mat)))
  }
  
  return(otu_mat)
}

# ============================================
# NMDS ANALYSIS FUNCTION
# ============================================
run_nmds_analysis <- function(otu_file, metadata_file, 
                              group_var = NULL,
                              transformation = "wisconsin",
                              distance_method = "bray",
                              k = 2, try = 40, trymax = 200,
                              filter_samples = FALSE, min_reads = 100,
                              filter_otus = TRUE, min_occurrence = 2,
                              filter_low_abundance = FALSE, min_abundance = 10) {
  
  # Read files
  cat("Reading OTU table...\n")
  otu_data <- read_uploaded_file(otu_file, basename(otu_file))
  
  cat("Reading metadata...\n")
  metadata <- read_uploaded_file(metadata_file, basename(metadata_file))
  
  # Check sample matching
  metadata_samples <- as.character(metadata[[1]])
  metadata_samples <- trimws(metadata_samples)
  otu_colnames <- names(otu_data)[-1]
  sample_cols <- intersect(otu_colnames, metadata_samples)
  
  if(length(sample_cols) == 0) {
    stop("Sample IDs don't match between OTU table and metadata!")
  }
  
  cat(sprintf("Found %d matching samples\n", length(sample_cols)))
  
  # Extract and process OTU data
  otu_samples <- otu_data[, sample_cols, drop = FALSE]
  otu_samples <- as.data.frame(lapply(otu_samples, function(x) {
    x <- suppressWarnings(as.numeric(as.character(x)))
    x[is.na(x)] <- 0
    return(x)
  }))
  
  # Transpose: rows = samples, columns = OTUs
  otu_mat <- as.matrix(t(otu_samples))
  rownames(otu_mat) <- sample_cols
  colnames(otu_mat) <- make.names(otu_data[[1]], unique = TRUE)
  
  # Apply filtering
  cat("\nApplying filters...\n")
  otu_mat <- filter_otu_data(otu_mat, 
                             filter_samples, min_reads,
                             filter_otus, min_occurrence,
                             filter_low_abundance, min_abundance)
  
  # Apply transformation
  cat(sprintf("\nApplying transformation: %s\n", transformation))
  otu_transformed <- apply_transformation(otu_mat, transformation)
  
  # Determine distance method
  if(transformation == "pa") {
    distance_method <- "jaccard"
    cat("Note: Using Jaccard distance for presence-absence data\n")
  }
  
  # Calculate distance matrix
  cat(sprintf("\nCalculating %s distance matrix...\n", distance_method))
  dist_matrix <- vegdist(otu_transformed, method = distance_method)
  
  if(length(dist_matrix) == 0 || any(is.na(dist_matrix))) {
    stop("Distance matrix invalid - try different settings")
  }
  
  # Run NMDS
  cat("\nRunning NMDS...\n")
  set.seed(123)
  nmds_result <- metaMDS(otu_transformed, 
                         distance = distance_method, 
                         k = k,
                         try = try, 
                         trymax = trymax,
                         autotransform = FALSE, 
                         trace = 1)
  
  # Print NMDS summary
  cat("\n========== NMDS Results ==========\n")
  cat(sprintf("Stress: %.4f\n", nmds_result$stress))
  cat(sprintf("Converged: %s\n", nmds_result$converged))
  cat(sprintf("Iterations: %d\n", nmds_result$iterations))
  
  # Statistical tests if grouping variable provided
  if(!is.null(group_var) && group_var %in% names(metadata)) {
    cat("\n========== Statistical Tests ==========\n")
    
    # Prepare metadata
    otu_samples <- rownames(otu_transformed)
    meta_df <- metadata
    names(meta_df)[1] <- "SampleID"
    meta_df$SampleID <- as.character(meta_df$SampleID)
    
    meta_filtered <- meta_df[meta_df$SampleID %in% otu_samples, ]
    meta_filtered <- meta_filtered[match(otu_samples, meta_filtered$SampleID), ]
    
    if(nrow(meta_filtered) > 0) {
      grouping_var <- meta_filtered[[group_var]]
      
      # Remove NA values
      na_idx <- is.na(grouping_var)
      if(any(na_idx)) {
        grouping_var <- grouping_var[!na_idx]
        keep_idx <- which(!na_idx)
        if(length(keep_idx) >= 3) {
          dist_subset <- as.dist(as.matrix(dist_matrix)[keep_idx, keep_idx])
        } else {
          warning("Not enough samples after removing NAs for statistical tests")
        }
      } else {
        dist_subset <- dist_matrix
      }
      
      if(exists("dist_subset") && length(unique(grouping_var)) >= 2) {
        grouping_factor <- as.factor(grouping_var)
        
        # PERMANOVA
        set.seed(123)
        permanova_result <- adonis2(dist_subset ~ grouping_factor, permutations = 999)
        cat("\nPERMANOVA Results:\n")
        print(permanova_result)
        
        # ANOSIM
        anosim_result <- anosim(dist_subset, grouping_factor, permutations = 999)
        cat("\nANOSIM Results:\n")
        print(anosim_result)
        
        # PERMDISP2
        betadisper_result <- betadisper(dist_subset, grouping_factor)
        permdisp_result <- permutest(betadisper_result, permutations = 999)
        cat("\nPERMDISP2 (Homogeneity of Dispersion) Results:\n")
        print(permdisp_result)
      }
    }
  }
  
  return(list(nmds = nmds_result, 
              distance_matrix = dist_matrix,
              transformed_data = otu_transformed,
              metadata = metadata,
              group_var = group_var,
              transformation = transformation,
              distance_method = distance_method))
}

# ============================================
# PLOTTING FUNCTION
# ============================================
plot_nmds <- function(nmds_results, 
                      point_size = 3,
                      show_ellipses = TRUE,
                      ellipse_alpha = 0.1,
                      color_palette = "Default") {
  
  # Extract NMDS scores
  scores_df <- as.data.frame(scores(nmds_results$nmds, display = "sites"))
  scores_df$SampleID <- rownames(scores_df)
  
  # Merge with metadata
  meta_df <- nmds_results$metadata
  names(meta_df)[1] <- "SampleID"
  meta_df$SampleID <- as.character(meta_df$SampleID)
  
  merged <- merge(scores_df, meta_df, by = "SampleID", all.x = TRUE)
  
  if(!is.null(nmds_results$group_var) && nmds_results$group_var %in% names(merged)) {
    merged$GroupingFactor <- as.factor(merged[[nmds_results$group_var]])
    
    # Remove NA groups
    merged <- merged[!is.na(merged$GroupingFactor), ]
    
    if(nrow(merged) == 0) {
      stop("No valid data to plot after removing NA groups")
    }
    
    # Create plot
    p <- ggplot(merged, aes(x = NMDS1, y = NMDS2, color = GroupingFactor)) +
      geom_point(size = point_size, alpha = 0.7) +
      theme_minimal() +
      labs(title = paste("NMDS Ordination - Grouped by", nmds_results$group_var),
           subtitle = paste("Transformation:", nmds_results$transformation, 
                            "| Distance:", nmds_results$distance_method),
           x = "NMDS1", y = "NMDS2", 
           color = nmds_results$group_var, 
           fill = nmds_results$group_var) +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
            legend.title = element_text(face = "bold"),
            panel.border = element_rect(fill = NA, color = "grey70"),
            panel.grid = element_line(color = "grey90"))
    
    # Add ellipses
    if(show_ellipses) {
      group_counts <- table(merged$GroupingFactor)
      valid_groups <- names(group_counts[group_counts >= 3])
      if(length(valid_groups) > 0) {
        merged_filtered <- merged[merged$GroupingFactor %in% valid_groups, ]
        p <- p + stat_ellipse(data = merged_filtered,
                              aes(x = NMDS1, y = NMDS2, fill = GroupingFactor),
                              level = 0.95, alpha = ellipse_alpha, 
                              linewidth = 0, geom = "polygon",
                              show.legend = FALSE)
      }
    }
    
    # Add color palette
    if(color_palette != "Default") {
      p <- p + scale_color_brewer(palette = color_palette) +
        scale_fill_brewer(palette = color_palette)
    }
    
    # Add stress value
    stress_val <- round(nmds_results$nmds$stress, 3)
    p <- p + annotate("text", x = Inf, y = Inf, 
                      label = paste("Stress:", stress_val),
                      hjust = 1.1, vjust = 1.1, 
                      size = 5, fontface = "bold")
    
    return(p)
  } else {
    stop("Grouping variable not found in metadata")
  }
}

# ============================================
# SHEPARD PLOT FUNCTION
# ============================================
plot_shepard <- function(nmds_results) {
  stressplot(nmds_results$nmds, 
             main = "Shepard Plot - Original vs. NMDS Distances")
}

# ============================================
# EXAMPLE USAGE
# ============================================

# Replace these file paths with your actual data files
otu_file <- "PATH/TO/OTUtable.csv"
metadata_file <- "PATH/TO/Metadata.csv"

# Run NMDS analysis
results <- run_nmds_analysis(
  otu_file = otu_file,
  metadata_file = metadata_file,
  group_var = "VARIABLE",  # Replace with your grouping column name
  transformation = "wisconsin", #metaMDS transformation method
  distance_method = "bray", # metaMDS distance calculation method
  k = 2, # metaMDS dimensions
  try = 40, 
  trymax = 200,
  filter_samples = FALSE, # filter samples by total read count
  min_reads = 100, # if filter = TRUE, minumum total reads per sample
  filter_otus = TRUE, # filter rare OTUs
  min_occurrence = 2, # minimum number of samples an OTU appears in
  filter_low_abundance = FALSE, # filter OTUs by abundance
  min_abundance = 10 # minimum total reads per OTU
)

# Create NMDS plot
nmds_plot <- plot_nmds(results, 
                       point_size = 3,
                       show_ellipses = TRUE,
                       ellipse_alpha = 0.1,
                       color_palette = "Default")

# Display the plot
print(nmds_plot)

# Create Shepard plot
shepard_plot=plot_shepard(results)

# Save plot to file
ggsave("nmds_plot.pdf", nmds_plot, width = 10, height = 8, dpi = 300)
ggsave("shepard_plot.pdf", shepard_plot, width = 10, height = 8, dpi = 300)