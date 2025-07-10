# #################################################################
# 
#        Presence or absence of differential expressed 
#            genes in different comparisons
#              
# #################################################################

# Avoid crashing Galaxy with a UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Debug: Print all received arguments
write(paste("Arguments reçus:", paste(args, collapse=", ")), file="debug_args.log")

# Function to get the value of an option from command line arguments
getOptionValue <- function(option) {
  idx <- match(option, args)
  if (!is.na(idx) && (idx + 1) <= length(args)) {
    return(args[idx + 1])
  }
  return(NULL)
}

# Input files
input_files_up_all <- getOptionValue("--galaxy_input_files_1")
input_files_down_all <- getOptionValue("--galaxy_input_files_2")
input_files_stat_all <- getOptionValue("--galaxy_input_files_3")
name_file_all_genes_stats <- getOptionValue("--galaxy_param1")
name_file_down <- getOptionValue("--galaxy_param2")
name_file_up <- getOptionValue("--galaxy_param3")
coldata <- getOptionValue("--galaxy_param4")

# Parameters from sections (modified names)
comp1_ref <- getOptionValue("--galaxy_input1")         # Ancien: input_file_1_ref
comp1_mutant <- getOptionValue("--galaxy_input2")      # Ancien: input_file_1_mutant 
comp1_regulation <- getOptionValue("--galaxy_input3")  # Ancien: input_up_or_down_1
comp2_ref <- getOptionValue("--galaxy_input4")         # Ancien: input_file_2_ref
comp2_mutant <- getOptionValue("--galaxy_input5")      # Ancien: input_file_2_mutant
comp2_regulation <- getOptionValue("--galaxy_input6")  # Ancien: input_up_or_down_2

output_path1 <- getOptionValue("--output_path1")
output_path2 <- getOptionValue("--output_path2")

# Debug: Print parameter values
write(paste(
  "Paramètres reçus:",
  "\nComp1 - Ref:", comp1_ref,
  "\nComp1 - Mutant:", comp1_mutant,
  "\nComp1 - Regulation:", comp1_regulation,
  "\nComp2 - Ref:", comp2_ref,
  "\nComp2 - Mutant:", comp2_mutant,
  "\nComp2 - Regulation:", comp2_regulation
), file="debug_params.log")

# Create output directories if they don't exist
if (!dir.exists(output_path1)) {
  dir.create(output_path1, recursive = TRUE)
}

if (!dir.exists(output_path2)) {
  dir.create(output_path2, recursive = TRUE)
}

# Process input file lists
input_files_up <- strsplit(input_files_up_all, ",")[[1]]
input_files_down <- strsplit(input_files_down_all, ",")[[1]]
input_files_stat <- strsplit(input_files_stat_all, ",")[[1]]

names_files_stats <- readLines(name_file_all_genes_stats)
names_files_down <- readLines(name_file_down)
names_files_up <- readLines(name_file_up)

create_verified_matrix <- function(file_names, file_list) {
  if (length(file_names) != length(file_list)) {
    stop("Error: The list of file names and the provided list have different lengths.")
  }
  matrix <- cbind(file_names, file_list)
  return(matrix)
}

# Create matrices for each list of files with length verification
tryCatch({
  up_matrix <- create_verified_matrix(names_files_up, input_files_up)
}, error = function(e) {
  print(paste("Error for UP matrix:", e$message))
})

tryCatch({
  down_matrix <- create_verified_matrix(names_files_down, input_files_down)
}, error = function(e) {
  print(paste("Error for DOWN matrix:", e$message))
})

tryCatch({
  stat_matrix <- create_verified_matrix(names_files_stats, input_files_stat)
}, error = function(e) {
  print(paste("Error for STAT matrix:", e$message))
})

# First comparison
name_input_file_1 <- paste0(comp1_ref, "_vs_", comp1_mutant, "_signif_", comp1_regulation, "_regulated.txt")

if (name_input_file_1 %in% up_matrix[, 1]) {
  index_file_1 <- which(up_matrix[, 1] == name_input_file_1)
  input_file_1 <- normalizePath(up_matrix[index_file_1, 2])
} else if (name_input_file_1 %in% down_matrix[, 1]) {
  index_file_1 <- which(down_matrix[, 1] == name_input_file_1)
  input_file_1 <- normalizePath(down_matrix[index_file_1, 2])
} else {
  stop(paste("Fichier", name_input_file_1, "non trouvé dans les matrices up/down"))
}

name_input_stat_file_1 <- paste0(comp1_ref, "_vs_", comp1_mutant, "_all_gene_stats.tsv")
if (name_input_stat_file_1 %in% stat_matrix[, 1]) {
  index_stat_file_1 <- which(stat_matrix[, 1] == name_input_stat_file_1)
  stat_input_file_1 <- normalizePath(stat_matrix[index_stat_file_1, 2])
} else {
  stop(paste("Fichier stats", name_input_stat_file_1, "non trouvé"))
}

# Second comparison  
name_input_file_2 <- paste0(comp2_ref, "_vs_", comp2_mutant, "_signif_", comp2_regulation, "_regulated.txt")

if (name_input_file_2 %in% up_matrix[, 1]) {
  index_file_2 <- which(up_matrix[, 1] == name_input_file_2)
  input_file_2 <- normalizePath(up_matrix[index_file_2, 2])
} else if (name_input_file_2 %in% down_matrix[, 1]) {
  index_file_2 <- which(down_matrix[, 1] == name_input_file_2)
  input_file_2 <- normalizePath(down_matrix[index_file_2, 2])
} else {
  stop(paste("Fichier", name_input_file_2, "non trouvé dans les matrices up/down"))
}

name_input_stat_file_2 <- paste0(comp2_ref, "_vs_", comp2_mutant, "_all_gene_stats.tsv")
if (name_input_stat_file_2 %in% stat_matrix[, 1]) {
  index_stat_file_2 <- which(stat_matrix[, 1] == name_input_stat_file_2)
  stat_input_file_2 <- normalizePath(stat_matrix[index_stat_file_2, 2])
} else {
  stop(paste("Fichier stats", name_input_stat_file_2, "non trouvé"))
}

# Process gene lists
file_1 = basename(name_input_file_1)
file_2 = basename(name_input_file_2)

gene_list_1 = read.delim(file = input_file_1,
                         header = FALSE,
                         stringsAsFactors = FALSE)

gene_list_2 = read.delim(file = input_file_2,
                         header = FALSE,
                         stringsAsFactors = FALSE)

## Data frame with common genes and corresponding fold change
gene_list_1 = gene_list_1[,-2]  # Remove the second column (basemean)
gene_list_2 = gene_list_2[,-2]

gene_list_1$foldchange_2 = gene_list_2[,2][match(gene_list_1[,1], gene_list_2[,1])]
gene_list_1$pvalue_2 = gene_list_2[,5][match(gene_list_1[,1], gene_list_2[,1])]
gene_list_1[,3] = gene_list_1[,5]
gene_list_1[,4] = gene_list_1$foldchange_2
gene_list_1[,5] = gene_list_1$pvalue_2

gene_list = gene_list_1[!is.na(gene_list_1[,4]),]
gene_list <- gene_list[-1,]
colnames(gene_list)[1] = "gene_id"
colnames(gene_list)[2] = paste("log2foldchange", file_1, sep = " ")
colnames(gene_list)[3] = paste("padj", file_1, sep = " ")
colnames(gene_list)[4] = paste("log2foldchange", file_2, sep = " ")
colnames(gene_list)[5] = paste("padj", file_2, sep = " ")

# Write common genes output
name1 <- paste0("common_genes_", comp1_ref, "_vs_", comp1_mutant, "_", comp1_regulation, 
                "_vs_", comp2_ref, "_vs_", comp2_mutant, "_", comp2_regulation)
output_file1 <- file.path(output_path1, paste0(name1, ".txt"))
write.table(gene_list[,1:5], file = output_file1, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Process statistics files
stat_file_1 = basename(name_input_stat_file_1)
stat_file_2 = basename(name_input_stat_file_2)

stat_1 = read.delim(file = stat_input_file_1,
                    header = FALSE,
                    stringsAsFactors = FALSE)

stat_2 = read.delim(file = stat_input_file_2,
                    header = FALSE,
                    stringsAsFactors = FALSE)

stat_1 = stat_1[,-2]
stat_2 = stat_2[,-2]

gene_list$foldchange_stat_1 = stat_1[,2][match(gene_list$gene_id, stat_1[,1])]
gene_list$pvalue_1 = stat_1[,5][match(gene_list$gene_id, stat_1[,1])]
gene_list$foldchange_stat_2 = stat_2[,2][match(gene_list$gene_id, stat_2[,1])]
gene_list$pvalue_2 = stat_2[,5][match(gene_list$gene_id, stat_2[,1])]

gene_list[,2] = gene_list$foldchange_stat_1
gene_list[,3] = gene_list$pvalue_1
gene_list[,4] = gene_list$foldchange_stat_2
gene_list[,5] = gene_list$pvalue_2

colnames(gene_list)[2] = paste("log2foldchange", stat_file_1, sep = " ")
colnames(gene_list)[3] = paste("padj", stat_file_1, sep = " ")
colnames(gene_list)[4] = paste("log2foldchange", stat_file_2, sep = " ")
colnames(gene_list)[5] = paste("padj", stat_file_2, sep = " ")

# Write statistics output
name2 <- paste0("stat_common_genes_", comp1_ref, "_vs_", comp1_mutant, "_", comp1_regulation,
                "_vs_", comp2_ref, "_vs_", comp2_mutant, "_", comp2_regulation)
output_file2 <- file.path(output_path2, paste0(name2, ".txt"))
write.table(gene_list[,1:5], file = output_file2, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Debug: Final confirmation
write("Traitement terminé avec succès", file="debug_completion.log")
