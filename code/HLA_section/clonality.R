subject_df <- read.csv('/subject_list/ms03.csv', sep='\t')

subject_list <- subject_df$subject_id 


rslt_df <- data.frame(stringsAsFactors = FALSE)

for (subject in subject_list) {
  temp_data <- ms03[[paste0(subject)]]
  
  clonality <- sqrt(sum((table(temp_data$CDR3.nt)/sum(table(temp_data$CDR3.nt)))^2))

  model <- data.frame(
    subject_id = subject,
    clonality = clonality,
    MS_status = 'MS',
    stringsAsFactors = FALSE
  )
  
  rslt_df <- rbind(rslt_df, model)
}

write.csv(rslt_df, '/ms03.csv', row.names=FALSE,  quote=FALSE)


## twins samples
subject_df <- read.csv('/meta_data/twins_therapy.csv', sep='\t')

subject_list <- subject_df$subject_id 

rslt_df <- data.frame(stringsAsFactors = FALSE)

for (subject in subject_list) {
  temp_data <- twins$data[[paste0(subject)]]
  
  clonality <- sqrt(sum((table(temp_data$CDR3.nt)/sum(table(temp_data$CDR3.nt)))^2))
  
  model <- data.frame(
    subject_id = subject,
    clonality = clonality,
    MS_status = 'MS',
    stringsAsFactors = FALSE
  )
  
  rslt_df <- rbind(rslt_df, model)
}

write.csv(rslt_df, '/twins.csv', row.names=FALSE,  quote=FALSE)


