#function that will take in filtered vcf file and construct sequences per individual and sequence for the outgroup species
make_sequences_for_temple = function(path_to_filtered_vcf_ib, path_to_filtered_vcf_ns, crm_coord){

  #using vk phylo to construct the fasta sequence for all accessions for both pops at the variant sites only
  #Note - if there are no variants detected in a population then the vk output will be null
  #I will not store the vk output per population in a file, instead will use intern
  
  #issuing vk phylo command for the two pops
  #CEU
  vk_phylo_ib = system(paste("vk phylo fasta ", path_to_filtered_vcf_ib, sep = ""), intern = TRUE)
  
  #CHS
  vk_phylo_ns = system(paste("vk phylo fasta ", path_to_filtered_vcf_ns, sep = ""), intern = TRUE)
  
  ### !!!!!!!! - check if the vcf file (both) are empty before proceeding - both need to have information and not just one and the other not !!!!!!!! ###
  if(length(vk_phylo_ib) > 0 & length(vk_phylo_ns) > 0){
    
    ## from this point on - only the CRMs for which there exists atleast one variant in both will be processed ##
    
    #extracting the reference sequence for the crm region now - not storing in the local storage
    crm_seq = system(paste("samtools faidx ",ingroup_ref," ",crm_coord$seqnames,":",crm_coord$start,"-",crm_coord$end, sep = ""), intern = TRUE)
    
    #have to adjust the output from samtools - in case of longer sequences, breaks it down into 60 sequence list
    crm_seq = paste(crm_seq[-1], collapse = "")
    
    ## now getting the outgroup species sequence for the alignment with the CRM ##
    alyra_seq = get_outgroup_sequence(crm_coord = crm_coord)
      
    #going ahead only if for the CRM I have atleast 20% information available in the alignment
    if(lengths(regmatches(alyra_seq, gregexpr("-", alyra_seq)))/nchar(alyra_seq) < 0.8){

      ## declaring an empty list which will fill up all the contents that are to be put in the fasta file ##
      temple_fasta_write = c()
      
      ##### sub-section for CEU pop #####
      
      ## pre-processing ##
      #extracting the individual ids
      pop_ids_ib = vk_phylo_ib[seq_len(length(vk_phylo_ib)) %%2 == 1]
      
      #processing to only keep the ids and remove the ">" sign
      pop_ids_ib = gsub(x = pop_ids_ib, pattern = ">", replacement = "")
      
      #extracting the varying site sequences
      var_seq_ib = vk_phylo_ib[seq_len(length(vk_phylo_ib)) %%2 == 0]
      
      #extracting the positions of these variants from the filtered vcf files (pop-specific)
      var_pos_ib = read.table(file = path_to_filtered_vcf_ib, skip = 16, header = FALSE, sep = "\t")
      
      #processing to keep only the second column containing the positions
      var_pos_ib = var_pos_ib$V2
      
      #converting these from genomic coordinates to crm-specific coordinates
      var_pos_ib = var_pos_ib - crm_coord$start + 1
      
      #constructing a matrix from the variant position and bases
      #BUG - I have to mention the matrix contruction has to be by row else it messes up a lot
      #when it goes by column by default the variant assignment per individual is lost and gives sign of variable sites
      var_info_ib = matrix(unlist(str_split(var_seq_ib, pattern = "")), nrow = length(pop_ids_ib), ncol = length(var_pos_ib), byrow = TRUE)
      
      #converting this into a dataframe - easier for me to access
      var_info_ib = as.data.frame(var_info_ib)
      
      #giving accession names as rownames and crm-specific positions as column names
      rownames(var_info_ib) = pop_ids_ib
      colnames(var_info_ib) = var_pos_ib
      
      #declaring the population
      temple_fasta_write = append(temple_fasta_write, "//pop_IB")
      
      #going over all accessions individually in a for loop and changing the variant position per individual
      #note - in case of missing data - vk phylo inserts an N - the missing data will be penalised by TEMPLE and all the inferences will be made by considering the missing data too
      for(ind in 1:nrow(var_info_ib)){
        
        #accessing variant sequence information of one individual from the population here
        ind = var_info_ib[ind, , drop = FALSE]
        
        #accessing the accession of this individual
        acc_ind = rownames(ind)
        
        #giving this accessing the wild type sequence, bases that are variants will be replaced from this sequence
        acc_seq = crm_seq
        
        #going over all the variant positions - going in a column-wise manner
        for(var_pos in 1:ncol(ind)){
          
          #creating the mutant sequence here - replacing the variant positions
          substr(acc_seq, start = as.numeric(colnames(ind)[var_pos]), stop = as.numeric(colnames(ind)[var_pos])) = ind[,var_pos]
          
        }
        
        #saving individual-specific id in the fasta file list
        temple_fasta_write = append(temple_fasta_write, paste(">",crm_coord$seqnames,"_",crm_coord$start, "_", crm_coord$end, "_IB",acc_ind, sep = ""))
        
        #saving individual-specific sequence in the fasta file list
        temple_fasta_write = append(temple_fasta_write, acc_seq)
        
      }
      
      #### sub-section for NS pop ####
      
      #first declaring that the file writing has entered the next pop
      temple_fasta_write = append(temple_fasta_write, "//pop_NS")
      
      ## pre-processing ##
      #extracting the individual ids
      pop_ids_ns = vk_phylo_ns[seq_len(length(vk_phylo_ns)) %%2 == 1]
      
      #processing to only keep the ids and remove the ">" sign
      pop_ids_ns = gsub(x = pop_ids_ns, pattern = ">", replacement = "")
      
      #extracting the varying site sequences
      var_seq_ns = vk_phylo_ns[seq_len(length(vk_phylo_ns)) %%2 == 0]
      
      #extracting the positions of these variants from the filtered vcf files (pop-specific)
      var_pos_ns = read.table(file = path_to_filtered_vcf_ns, skip = 16, header = FALSE, sep = "\t")
      
      #processing to keep only the second column containing the positions
      var_pos_ns = var_pos_ns$V2
      
      #converting these from genomic coordinates to crm-specific coordinates
      var_pos_ns = var_pos_ns - crm_coord$start + 1
      
      #constructing a matrix from the variant position and bases
      #BUG - I have to mention the matrix contruction has to be by row else it messes up a lot
      #when it goes by column by default the variant assignment per individual is lost and gives sign of variable sites
      var_info_ns = matrix(unlist(str_split(var_seq_ns, pattern = "")), nrow = length(pop_ids_ns), ncol = length(var_pos_ns), byrow = TRUE)
      
      #converting this into a dataframe - easier for me to access
      var_info_ns = as.data.frame(var_info_ns)
      
      #giving accession names as rownames and crm-specific positions as column names
      rownames(var_info_ns) = pop_ids_ns
      colnames(var_info_ns) = var_pos_ns
      
      #going over all accessions individually in a for loop and changing the variant position per individual
      #note - in case of missing data - vk phylo inserts an N - the missing data will be penalised by TEMPLE and all the inferences will be made by considering the missing data too
      for(ind in 1:nrow(var_info_ns)){
        
        #accessing variant sequence information of one individual from the population here
        ind = var_info_ns[ind, , drop = FALSE]
        
        #accessing the accession of this individual
        acc_ind = rownames(ind)
        
        #giving this accessing the wild type sequence, bases that are variants will be replaced from this sequence
        acc_seq = crm_seq
        
        #going over all the variant positions - going in a column-wise manner
        for(var_pos in 1:ncol(ind)){
          
          #creating the mutant sequence here - replacing the variant positions
          substr(acc_seq, start = as.numeric(colnames(ind)[var_pos]), stop = as.numeric(colnames(ind)[var_pos])) = ind[,var_pos]
          
        }
        
        #saving individual-specific id in the fasta file list
        temple_fasta_write = append(temple_fasta_write, paste(">",crm_coord$seqnames,"_",crm_coord$start, "_", crm_coord$end, "_NS",acc_ind, sep = ""))
        
        #saving individual-specific sequence in the fasta file list
        temple_fasta_write = append(temple_fasta_write, acc_seq)
        
      }
      
      #writing the outgroup sequence now - as per TEMPLE manual the outgroup sequence should always be inserted in the second pop
      #inserting it at the end of second pop - this is how it is shown in the example within the TEMPLE manual page12
      temple_fasta_write = append(temple_fasta_write, paste(">OUTGR",crm_coord$seqnames,"_",crm_coord$start,"_",crm_coord$end,"_WT", sep = ""))
      temple_fasta_write = append(temple_fasta_write, alyra_seq)
      
      #closing the fasta file writing now
      temple_fasta_write = append(temple_fasta_write, "//")
      
      #finally, writing the fasta file
      writeLines(text = temple_fasta_write, paste("backup_files/sequence_files/crm_", crm_coord$seqnames,"_",crm_coord$start, "_", crm_coord$end, ".fa", sep = ""))
      
      #return statement
      return(df = data.frame())
      
    }else{
      
      #if you are inside it means that you have stumbled on a CRM that has a gap ratio of 0.8 or more - these are relevant so listing them
      df = crm_coord
      
      #return statement
      return(df)
      
    }
  }else{
    
    #if you are inside it means than there were no variants in one or both populations
    #here too I return an empty dataframe - have to do it for purpose of consistency
    return(df = data.frame())
    
  }
}