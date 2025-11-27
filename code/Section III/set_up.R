library( plyr )
library("scales")   
library(dplyr)
library(lsr)
library(stringr)

hla_df <- read.csv('E:/Document/meta_data/twins_hla_typing.csv', sep='\t')
subject_list <- hla_df$subject_id

res <- list()

for (subject in subject_list) {
  file_path <- paste0("D:/large_files/PBMC_TWINS/TWIN_", subject, ".csv")
  df <- read.table(file_path, sep = "\t", header = TRUE)
  df1 <- df %>%
    filter(!grepl("OR", v_family)) %>%
    filter(!grepl("N", v_family)) %>%
    filter(!grepl("Out", frame_type)) %>%
    filter(!grepl("Stop", frame_type)) %>%
    filter(!grepl("\\*", amino_acid))
  df1 <- df1[!(is.na(df1$v_resolved) | df1$v_resolved==""), ]
  
  column_names <- c('seq', 'v_resolved', 'amino_acid', 'j_resolved', 
                    'productive_frequency', 'templates')
  df1 <- df1 %>%
    select(all_of(column_names))
  
  sampled_data <- head(df1, 50000)
  sampled_data$subject_id <- subject
  res <- rbind(res, sampled_data)
}
hla_df <- read.csv('E:/Document/meta_data/twins_hla_typing.csv', sep='\t')
pos_df <- hla_df[hla_df[['HLA_DRB1.15']] == 1, ]
neg_df <- hla_df[hla_df[['HLA_DRB1.15']] == 0, ]

pos_list <- pos_df$subject_id
neg_list <- neg_df$subject_id

pos_res <- res[res$subject_id %in% pos_list, ]
neg_res <- res[res$subject_id %in% neg_list, ]

pos_res1 <- pos_res[pos_res$v_resolved == 'TRBV05-01*01',]
# pos_res1 <- pos_res1[pos_res1$j_resolved == 'TRBJ02-07*01',]

neg_res1 <- neg_res[neg_res$v_resolved == 'TRBV05-01*01',]
# neg_res1 <- neg_res1[neg_res1$j_resolved == 'TRBJ02-07*01',]

#initialize matrix
min_clength <- max(min(nchar(pos_res1$amino_acid)), min(nchar(neg_res1$amino_acid)))
max_clength <- min(max(nchar(pos_res1$amino_acid)), max(nchar(neg_res1$amino_acid)))

cmat_chi <- matrix(, nrow = max_clength, ncol = max_clength)
cmat_cramer <- matrix(, nrow = max_clength, ncol = max_clength)

residuals <- data.frame(matrix(ncol = 6, nrow = 0))

n_pos <- NULL
n_neg <- NULL

for (n in min_clength:max_clength){
  print(n)
  n_length_aa_pos <- pos_res1$amino_acid[nchar(pos_res1$amino_acid) == n]
  n_length_aa_neg <- neg_res1$amino_acid[nchar(neg_res1$amino_acid) == n]
  n_neg <- append(n_neg, length(n_length_aa_neg))
  n_pos <- append(n_pos, length(n_length_aa_pos))
  for(pos in 1:n){
    # print(pos)
    aa_count_pos_n <- str_sub(n_length_aa_pos,pos,pos) %>% 
      paste(collapse="") %>% 
      strsplit(split="") %>% unlist %>%
      table %>% data.frame
    aa_count_pos_n2 <- data.frame(aa_count_pos_n[,-1])
    rownames(aa_count_pos_n2) <- aa_count_pos_n[,1]
    
    aa_count_neg_n <- str_sub(n_length_aa_neg,pos,pos) %>% 
      paste(collapse="") %>% 
      strsplit(split="") %>% unlist %>%
      table %>% data.frame
    aa_count_neg_n2 <- data.frame(aa_count_neg_n[,-1])
    rownames(aa_count_neg_n2) <- aa_count_neg_n[,1]
    if (nrow(aa_count_neg_n2) == 0 || nrow(aa_count_pos_n2) == 0){
      cmat_chi[pos, n] <- 1
      cmat_cramer[pos, n] <- 0
      print("empty df")
    }
    else{
      df3<- merge(aa_count_pos_n2, aa_count_neg_n2, by = "row.names", all=T)
      df <- data.frame(df3[,-1])
      rownames(df) <- df3[,1]
      df[is.na(df)] = 0
      colnames(df) <- c("pos", "neg")
      
      chi <- chisq.test(as.table(t(df)))
      cramers_v <- cramersV(as.table(t(df)))
      if(pos<=2 | pos >= (n-2)){
        cmat_chi[pos, n] <- 1
        cmat_cramer[pos, n] <- 0
      }
      else{    
        cmat_chi[pos, n] <- chi$p.value
        cmat_cramer[pos, n] <- cramers_v
      }
      temp_df <- data.frame("length" =n, "position"=pos)
      if(length(chi$residuals) <= 2){
        #print(c(pos,n))
        r <- cbind(names = rownames(t(df)), t(df))
        r[,2] <- colnames(r)[2]
        n_residuals <- cbind(r, chi$stdres)
        colnames(n_residuals) <- c("celltype", "aa", "residual")
      }else{
        n_residuals <- data.frame(chi$stdres)
        colnames(n_residuals) <- c("celltype", "aa", "residual")
      }
      temp_df <- cbind(cbind(cbind(temp_df, n_residuals), chi$p.value),cramers_v)
      residuals <- rbind(residuals, temp_df)
    }}
}

write.table(residuals, "/twins_residuals.txt")
write.table(cmat_chi, "/mat_chi_all.txt")
write.table(cmat_cramer, "/mat_cramer_all.txt")
write.table(n_neg, "/neg_n_all.txt")
write.table(n_pos, "/pos_n_all.txt")

