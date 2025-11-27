library(readr)

temp_matrix <- read_tsv("/output/hla_ms_seqs.tsv")
temp_matrix_other <- read_tsv("/user_matrix/output_matrix/mbp/background.tsv") # , show_col_types = FALSE

temp_matrix <- as.data.frame(temp_matrix)
temp_matrix_other <- as.data.frame(temp_matrix_other)

rownames(temp_matrix) <- temp_matrix[[1]]
temp_matrix <- temp_matrix[,-1]

rownames(temp_matrix_other) <- temp_matrix_other[[1]]
temp_matrix_other <- temp_matrix_other[,-1]

temp_matrix[] <- lapply(temp_matrix, as.numeric)
temp_matrix_other[] <- lapply(temp_matrix_other, as.numeric)

temp_matrix_other <- temp_matrix_other[, colnames(temp_matrix)]

for (i in 1:ncol(temp_matrix)) {
  bg_mean <- mean(temp_matrix_other[,i], na.rm = TRUE)
  temp_matrix[,i] <- temp_matrix[,i] - bg_mean
}

for (i in 1:nrow(temp_matrix)) {
  row_vals <- as.numeric(temp_matrix[i, ])    # Force to numeric vector
  row_mean <- mean(row_vals, na.rm = TRUE)
  temp_matrix[i,] <- temp_matrix[i,] - row_mean
}

#print(temp_matrix)
write.csv(temp_matrix, '/twins_binding_score.csv', row.names=TRUE, col.names = TRUE) #, quote=FALSE
colMeans(temp_matrix)
