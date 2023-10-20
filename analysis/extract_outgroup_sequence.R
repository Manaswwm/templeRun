#short script that takes in CRM coordinates and gives out already constructed CRM sequence for the outgroup species
extract_outgroup_sequence = function(crm_filename){
  
  #this function takes in an already constructed CRM file name that is based on the chr and start and stop positions of the CRM
  #applying a dirty trick here of going into the sequence file directories of every backup file and asking if the file exists within one of them
  #in the batch where the condition is true - the sequence is extracted and returned
  
  ### messy if statements to extract the sequence file directly - checking all the batch directories
  
  #batch1
  if(file.exists(paste("/netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project_new/aratha/temple_batches_nopval/batch_1/sequence_files/", crm_filename, sep = ""))){
    
    #if you are inside it means the sequence file is found - now reading it
    seq_file = readLines(con = paste("/netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project_new/aratha/temple_batches_nopval/batch_1/sequence_files/", crm_filename, sep = ""))
    
    #the outgroup sequence is always the 184th entry
    out_seq = seq_file = seq_file[184]
    
  }
  #batch1
  if(file.exists(paste("/netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project_new/aratha/temple_batches_nopval/batch_2/sequence_files/", crm_filename, sep = ""))){
    
    #if you are inside it means the sequence file is found - now reading it
    seq_file = readLines(con = paste("/netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project_new/aratha/temple_batches_nopval/batch_2/sequence_files/", crm_filename, sep = ""))
    
    #the outgroup sequence is always the 184th entry
    out_seq = seq_file = seq_file[184]
    
  }
  #batch1
  if(file.exists(paste("/netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project_new/aratha/temple_batches_nopval/batch_3/sequence_files/", crm_filename, sep = ""))){
    
    #if you are inside it means the sequence file is found - now reading it
    seq_file = readLines(con = paste("/netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project_new/aratha/temple_batches_nopval/batch_3/sequence_files/", crm_filename, sep = ""))
    
    #the outgroup sequence is always the 184th entry
    out_seq = seq_file = seq_file[184]
    
  }
  
  #returning the outgroup sequence
  return(out_seq)
  
}