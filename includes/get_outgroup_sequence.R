#function that takes in the CRM coordinates and get the corresponding outgroup sequences
#in case the alignment is absent for the part of the CRM, gaps are introduced artifically
get_outgroup_sequence = function(crm_coord){
  
  ### subsection 1 - creating the outgroup sequence first ###
  ## REST-API section - giving the coordinates of the CRM and obtaining the alignment to lyrata directly
  #building query
  api_query_pt1 = paste("https://rest.ensembl.org/alignment/region/arabidopsis_thaliana/",crm_coord$seqnames,":",crm_coord$start,"-",crm_coord$end,"?",sep = "")
  api_query_pt2 = paste("method=LASTZ_NET&compact=1;species_set=arabidopsis_thaliana;species_set=arabidopsis_lyrata;compara=plants")
  
  #firing query
  api_response = GET(paste(api_query_pt1, api_query_pt2, sep = ""), content_type("application/json"))
  
  #status check
  if(api_response$status_code == 200){
    
    #filtering
    alignments = lapply(content(api_response), list.remove, "tree")
    
    #### Problem - in an ideal world - for the CRM region in A. thaliana - there will be an end to end alignment in A. lyarata, but since we dont have unicorns that is not the case ####
    ## Presumption - for a majority of CRMs there will be an imperfect alignment(s), potentially multiple alignments of different lengths ##
    ## I will have to choose the alignment that is the longest ##
    
    ## resolving the choice of alignment here - its a bit tedious as the output from ensembl has a complex data structure
    
    #creating a function that could go over all the possible alignments and subset them and convert to a dataframe
    clean_alignments = function(alignment){
      
      #accessing only specific information per-alignment - lot going on in one line
      #in a nutshell - per alignment I am keeping track of the start pos, end pos and sequence for both species and binding two as rows in a dataframe
      df = as.data.frame(do.call(rbind, lapply(alignment$alignments, function(x) {x[c("species", "start", "end", "seq")]})))
      
      #print check
      return(df)
      
    }
    
    #getting a list of all the alignments as dataframe
    alignments_df = lapply(alignments, function(x) {clean_alignments(x)})
    
    #getting the alignment with the longest length - lot going on in one line
    longest_alignment = alignments_df[which.max(lapply(alignments_df, function(x){abs(unlist(x$end[x$species == "arabidopsis_thaliana"]) - unlist(x$start[x$species == "arabidopsis_thaliana"]))}))][[1]]
    
    #here columns and elements are list - converting them to elements
    longest_alignment = as.data.frame(t(apply(longest_alignment, 1, unlist)))
    
    #converting columns to respective data-type
    longest_alignment$start = as.numeric(longest_alignment$start)
    longest_alignment$end = as.numeric(longest_alignment$end)
    
    #processing sequences, chopping them into individual characters so that I can make matrix from them
    aratha_seq = unlist(strsplit(longest_alignment$seq[longest_alignment$species == "arabidopsis_thaliana"], split = ""))
    alyra_seq = unlist(strsplit(longest_alignment$seq[longest_alignment$species == "arabidopsis_lyrata"], split = ""))
    
    #making alignment matrix
    m = matrix(c(aratha_seq, alyra_seq), nrow = 2, byrow = TRUE)
    
    #checking if there are gaps in thaliana
    if(str_count(paste(m[1,], collapse = ""), pattern = "-") > 0){
      
      #deleting gaps in thaliana
      m_no_gaps_aratha = m[,-which(m[1,] == "-")]
      
    }else{m_no_gaps_aratha = m} #if no gaps
    
    #adding rownames and column names
    rownames(m_no_gaps_aratha) = c("arabidopsis_thaliana", "arabidopsis_lyrata")
    colnames(m_no_gaps_aratha) = seq(longest_alignment$start[longest_alignment$species == "arabidopsis_thaliana"], longest_alignment$end[longest_alignment$species == "arabidopsis_thaliana"], by = 1)
    
    #making an empty matrix whose dimensions are same as that of the CRM
    #I add plus 1 as this accounts for R not counting the first position as the actual first position
    sequence_matrix = matrix(nrow = 2, ncol = crm_coord$end - crm_coord$start + 1)
    
    #filling in gaps as all positions first
    sequence_matrix[,] = "-"
    
    #declaring row names and column names - have to be consistent with the other matrix
    rownames(sequence_matrix) = c("arabidopsis_thaliana", "arabidopsis_lyrata")
    colnames(sequence_matrix) = seq(crm_coord$start, crm_coord$end, by = 1)
    
    #since m_thaliana_no_gaps is the matrix with information, filling in the information into sequence matrix using a for loop
    #going over all column names - positions where there is alignment - the remaining shall remain gaps
    for(pos in colnames(m_no_gaps_aratha)){
      
      #replacing values in both species by specifically looking at individual positions
      sequence_matrix[1,pos] = m_no_gaps_aratha[1,pos]
      sequence_matrix[2,pos] = m_no_gaps_aratha[2,pos]
      
    }
    
    #re-accessing lyrata sequence - the final stage is accessing the outgroup sequence
    alyra_seq = paste(sequence_matrix[2,], collapse = "")
    
    #returning the outgroup sequence here
    return(alyra_seq)
    
  }else{
    
    #if you are inside then it means the output from ensembl was not as expected (status code was not 200)
    print(paste("For CRM coordinates : ",crm_coord$seqnames,":",crm_coord$start,"-",crm_coord$end," the return from Ensembl did not have a normal status code", sep = ""))
    
  }
  
}