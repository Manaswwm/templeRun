#this function takes in the crm_coordinates of a single CRM and constructs a region file that will be used by TEMPLE
make_regionfile_for_temple = function(crm_coord){
  
  #extracting the reference sequence again - I need to be able to take this from other function - too lazy
  crm_seq = system(paste("samtools faidx ", ingroup_ref, " ",crm_coord$seqnames,":",crm_coord$start,"-",crm_coord$end, sep = ""), intern = TRUE)
  
  #have to adjust the output from samtools - in case of longer sequences, breaks it down into 60 sequence list
  crm_seq = paste(crm_seq[-1], collapse = "")
  
  #making a dataframe of the coordinate information - formatting and adding other information that is compatible with TEMPLE
  region_file_df = data.frame(set = "pop_IB&pop_NS", chr = crm_coord$seqnames, start = 1, stop = nchar(crm_seq), id = "crm", pwm = "all")
  
  #writing to a file now
  write.table(x = region_file_df, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", file = paste("backup_files/region_files/crm_",crm_coord$seqnames,"_",crm_coord$start,"_",crm_coord$end,"_regionfile.txt", sep = ""))
  
}