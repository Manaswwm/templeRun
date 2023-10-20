### script for constructing lengths by employing the ancestral alleles and ratio_metric score per variant
## this script first calculate ratio-metric score per variant in both populations and then per position in the identified PWM per populations

#importing relevant libraries
library(ggplot2)
library(Biostrings)
library(TFBSTools)

#sourcing functions
source("extract_outgroup_sequence.R")

#importing the mutation output file from TEMPLE that is modified slightly in the output_track.R script
mutation_output_ib = read.delim(file = "pop_specific_variants/ib_pop_mutation_subset.txt", header = TRUE, sep = "\t")
mutation_output_ns = read.delim(file = "pop_specific_variants/ns_pop_mutation_subset.txt", header = TRUE, sep = "\t")

#importing the tfbs counts per pwm
tfbs_counts_info = read.delim(file = "pwm_info/pwm_counts.txt", header = TRUE, sep = "\t")

##### extracting the allele counts of all variants #####

### first pop is YRI

#creating an empty dataframe to fill in the annotations
mutation_annotate_ib = data.frame()

#going over every variant sequentially
for(row in 1:nrow(mutation_output_ib)){
  
  #@print check
  print("Going over all variants in the IB population ... ")
  
  #selecting a single variant
  var = mutation_output_ib[row,]
  
  #### from TEMPLE's manual - the PosInTFBS column specifically states the position of the variant with reference of the positive strand
  #### post confirming if the reported positions of the variants on the reverse strand actually contain the polarized alleles I found that they indeed were
  #### to summarise - TEMPLE's reporting of variants is in reference to the forward strand even if the variants occur on the reverse strand
  
  #extracting the ref and alt allele
  ref_allele = as.character(var$Al1)
  alt_allele = as.character(var$Al2)
  
  #### reverse-stranded alleles vs positive stranded alleles - in the output all the alleles that are listed under reverse = true are actually occuring on the REVERSE strand
  #### alleles that are listed under reverse = FALSE are actually occuring on the POSITIVE/SENSE strand
  #### for alleles that are occurring on the positive strand the logic of obtaining the ref and alt allele counts is straightforward - for reverse stranded alleles it is a bit COMPLICATED
  #### for reverse stranded variants - I will have to obtain the reverse complemented position and the allele counts will be the allele counts of the reverse complemented alleles on the forward strand
  #### so for variant on position 9 on the reverse strand with AL1 -> A; and AL2 -> G - I obtain the reverse complemented position on the forward strand (in a PWM of length 12 this will be 9 --> 4)
  #### then on the reverse complemented position the allele counts of reverse complemented alleles will be given; so count of allele A will be the count of allele T and count of allele G will be count of allele C
  if(var$Reverse == "true"){ #if you are inside it means you are on the reverse strand
    
    #extracting the coutns info the specific pwm of interest and changing the row names
    counts_pos = tfbs_counts_info[tfbs_counts_info$name == var$PWM,]
    rownames(counts_pos) = NULL
    
    #since the allele is on the reverse strand - extracting the reverse complemented position
    counts_pos_reverse = (dim(counts_pos)[1] - var$PosInTFBS) + 1
    counts_pos = counts_pos[counts_pos_reverse,]
    
    ## extracting the count of ref and alt allele
    ref_allele_count = as.numeric(counts_pos[colnames(counts_pos) == as.character(reverseComplement(DNAString(ref_allele)))]) + 1
    alt_allele_count = as.numeric(counts_pos[colnames(counts_pos) == as.character(reverseComplement(DNAString(alt_allele)))]) + 1
    
  }else{ #if you inside it means you are on the positive strand
    
    #extracting the coutns info the specific pwm of interest and changing the row names
    counts_pos = tfbs_counts_info[tfbs_counts_info$name == var$PWM,]
    rownames(counts_pos) = NULL
    
    #extracting the position within the TFBS
    counts_pos = counts_pos[var$PosInTFBS,]
    
    #extracting the count of ref and alt alleles
    ref_allele_count = as.numeric(counts_pos[colnames(counts_pos) == ref_allele]) + 1
    alt_allele_count = as.numeric(counts_pos[colnames(counts_pos) == alt_allele]) + 1
    
  }
  
  #writing the information into a dataframe and merging
  df = var
  
  #adding the counts of ref and alt alleles
  df$ref_allele_count = ref_allele_count
  df$alt_allele_count = alt_allele_count
  
  #calculating the log odds score
  df$score = log((alt_allele_count/1000)/0.25) - log((ref_allele_count/1000)/0.25) 
  
  #merging
  mutation_annotate_ib = rbind(mutation_annotate_ib, df)
  
  #@print check
  print(paste("(Part 1) For IB population variants - Gone over row : ",row, " of total rows : ",dim(mutation_output_ib)[1],sep = ""))
  
}

#cleaning up
rm(counts_pos, df, var, alt_allele, alt_allele_count, counts_pos_reverse, ref_allele, ref_allele_count, row)

#removing NA columns -- please check later what is wrong here
mutation_annotate_ib = mutation_annotate_ib[!is.na(mutation_annotate_ib$ref_allele_count),]

### important part - checking if my calculated score is similar to the score predicted by TEMPLE - if not entirely similar atleast that they have a positive correlation ###
cor.test(mutation_annotate_ib$EffectOnNucleotide, mutation_annotate_ib$score)

#calculating ratio score for all the listed variants
mutation_annotate_ib$ratio_score = by(mutation_annotate_ib, seq_len(nrow(mutation_annotate_ib)), function(x) {abs(x$ref_allele_count - x$alt_allele_count)/max(x$ref_allele_count,x$alt_allele_count)})

### second pop is SWE
#creating an empty dataframe to fill in the annotations
mutation_annotate_ns = data.frame()

#going over every variant sequentially
for(row in 1:nrow(mutation_output_ns)){
  
  #@print check
  print("Going over all variants in the NS population ... ")
  
  #selecting a single variant
  var = mutation_output_ns[row,]
  
  #### from TEMPLE's manual - the PosInTFBS column specifically states the position of the variant with reference of the positive strand
  #### post confirming if the reported positions of the variants on the reverse strand actually contain the polarized alleles I found that they indeed were
  #### to summarise - TEMPLE's reporting of variants is in reference to the forward strand even if the variants occur on the reverse strand
  
  #extracting the ref and alt allele
  ref_allele = as.character(var$Al1)
  alt_allele = as.character(var$Al2)
  
  #### reverse-stranded alleles vs positive stranded alleles - in the output all the alleles that are listed under reverse = true are actually occuring on the REVERSE strand
  #### alleles that are listed under reverse = FALSE are actually occuring on the POSITIVE/SENSE strand
  #### for alleles that are occurring on the positive strand the logic of obtaining the ref and alt allele counts is straightforward - for reverse stranded alleles it is a bit COMPLICATED
  #### for reverse stranded variants - I will have to obtain the reverse complemented position and the allele counts will be the allele counts of the reverse complemented alleles on the forward strand
  #### so for variant on position 9 on the reverse strand with AL1 -> A; and AL2 -> G - I obtain the reverse complemented position on the forward strand (in a PWM of length 12 this will be 9 --> 4)
  #### then on the reverse complemented position the allele counts of reverse complemented alleles will be given; so count of allele A will be the count of allele T and count of allele G will be count of allele C
  if(var$Reverse == "true"){ #if you are inside it means you are on the reverse strand
    
    #extracting the coutns info the specific pwm of interest and changing the row names
    counts_pos = tfbs_counts_info[tfbs_counts_info$name == var$PWM,]
    rownames(counts_pos) = NULL
    
    #since the allele is on the reverse strand - extracting the reverse complemented position
    counts_pos_reverse = (dim(counts_pos)[1] - var$PosInTFBS) + 1
    counts_pos = counts_pos[counts_pos_reverse,]
    
    ## extracting the count of ref and alt allele
    ref_allele_count = as.numeric(counts_pos[colnames(counts_pos) == as.character(reverseComplement(DNAString(ref_allele)))]) + 1
    alt_allele_count = as.numeric(counts_pos[colnames(counts_pos) == as.character(reverseComplement(DNAString(alt_allele)))]) + 1
    
  }else{ #if you inside it means you are on the positive strand
    
    #extracting the coutns info the specific pwm of interest and changing the row names
    counts_pos = tfbs_counts_info[tfbs_counts_info$name == var$PWM,]
    rownames(counts_pos) = NULL
    
    #extracting the position within the TFBS
    counts_pos = counts_pos[var$PosInTFBS,]
    
    #extracting the count of ref and alt alleles
    ref_allele_count = as.numeric(counts_pos[colnames(counts_pos) == ref_allele]) + 1
    alt_allele_count = as.numeric(counts_pos[colnames(counts_pos) == alt_allele]) + 1
    
  }
  
  #writing the information into a dataframe and merging
  df = var
  
  #adding the counts of ref and alt alleles
  df$ref_allele_count = ref_allele_count
  df$alt_allele_count = alt_allele_count
  
  #calculating the log odds score
  df$score = log((alt_allele_count/1000)/0.25) - log((ref_allele_count/1000)/0.25) 
  
  #merging
  mutation_annotate_ns = rbind(mutation_annotate_ns, df)
  
}

#cleaning up
rm(counts_pos, df, var, alt_allele, alt_allele_count, counts_pos_reverse, ref_allele, ref_allele_count, row)

#removing NA columns -- please check later what is wrong here
mutation_annotate_ns = mutation_annotate_ns[!is.na(mutation_annotate_ns$ref_allele_count),]

### important part - checking if my calculated score is similar to the score predicted by TEMPLE - if not entirely similar atleast that they have a positive correlation ###
cor.test(mutation_annotate_ns$EffectOnNucleotide, mutation_annotate_ns$score)

#calculating shannon's entropy for all the listed variants
mutation_annotate_ns$ratio_score = by(mutation_annotate_ns, seq_len(nrow(mutation_annotate_ns)), function(x) {abs(x$ref_allele_count - x$alt_allele_count)/max(x$ref_allele_count,x$alt_allele_count)})

#writing all info to a file now
#first all variant specific entropy info
write.table(x = mutation_annotate_ib, file = "ratios_info/mutation_annotate_info_ib.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = mutation_annotate_ns, file = "ratios_info/mutation_annotate_info_ns.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


####### PWM lengths for IB ########

#declaring empty dataframes to fill in allele-specific information - please note the alleles here are in no particular order of preference just that pairs are maintained
pwm_allele1_df_ib = data.frame()
pwm_allele2_df_ib = data.frame()
pwm_allele3_df_ib = data.frame()

#going over every PWM and their respective positions within the alignments
for(row in 1:nrow(mutation_annotate_ib)){
  
  #extracting the var information
  var = mutation_annotate_ib[row,]
  
  ##### Step 1 - Extracting the CRM sequence info and the TFBS sequence in question within this CRM #####
  
  #extracting the coordinates of the CRM
  crm_chr = var$chr
  crm_start = var$start_crm
  crm_end = var$stop_crm
  
  #what the crm sequence filename will be
  crm_filename = paste("crm_",crm_chr,"_",crm_start,"_",crm_end,".fa", sep = "")
  
  #giving a call to the function to extract the outgroup sequence directly from the sequence files
  out_seq = extract_outgroup_sequence(crm_filename)
  
  #extracting the sequence where the TFBS is predicted to be binding to now
  ## start of the TFBS ==> (var$Position - Var$PosInTFBS + 1) ------ var$Position - var$PosInTFBS tells me where the first TFBS position is, adding 1 converts bed0 to bed1
  ## end of TFBS ==> start of TFBS + length of the PWM
  pwm_seq = substr(x = out_seq, start = var$Position - var$PosInTFBS + 1, stop = (var$Position - var$PosInTFBS + 1) + dim(tfbs_counts_info[tfbs_counts_info$name == var$PWM,])[1] - 1)
  
  ### the problem of reverse complementing DNA sequences 
  ### if variants are found on the reverse strand (which is recognized by "reverse == TRUE") then the whole PWM is by default recognized on the reverse strand
  ### in such cases I will also have to reverse complement the outgroup seqeuence to obtain the ancestral allele for the reverse complemented positions
  if(var$Reverse == "true"){
    
    #if you are inside it means that the pwm is identified on the reverse complemented strand - so reverse complementing the ancestral sequence
    pwm_seq = as.character(reverseComplement(DNAString(pwm_seq)))
    
  } #if Reverse == FALSE then continue without deroute
  
  #extracting the tfbs for the pwm
  pwm_info = tfbs_counts_info[tfbs_counts_info$name == var$PWM, c("A", "C", "G", "T", "name")]
  rownames(pwm_info) = NULL
  
  ###### In some weird cases for drosophila (possibly for other species too) the outgroup sequences is terminated in between for certain PWMs
  ###### controlling for that by checking if the length of the outgroup PWM sequence is similar to the pwm length
  if(dim(pwm_info)[1] == nchar(pwm_seq)){
    
    #adding the ancestral state information here
    pwm_info$anc_allele = unlist(strsplit(pwm_seq, ""))
    
    ##### Step 2 - Calculating the pairwise entropies per position per PWM - every position will have 3 entropies #####
    
    #first checking and removing any pwm positions that have a gap in the ancestral allele
    pwm_info = pwm_info[!pwm_info$anc_allele == "-",]
    
    #if there is no information left then bypassing this pwm
    if(dim(pwm_info)[1] > 0){
      
      #going over every position in a for loop
      for(pos in 1:nrow(pwm_info)){
        
        #extracting the pos information
        pos = pwm_info[pos,]
        
        #listing all the other alleles that are not ancestral alleles
        var_alleles = colnames(pos)[1:4][!colnames(pos)[1:4] == pos$anc_allele]
        
        #storing counts of all alleles and adding 1 to them
        anc_count = as.numeric(pos[colnames(pos) == pos$anc_allele]) + 1
        allele1_count = as.numeric(pos[colnames(pos) == var_alleles[1]]) + 1
        allele2_count = as.numeric(pos[colnames(pos) == var_alleles[2]]) + 1
        allele3_count = as.numeric(pos[colnames(pos) == var_alleles[3]]) + 1
        
        #calculating the pairwise entropy - here allele numbering is in order of A,C,G and T - not that it matters
        allele1_ratio_score = abs(anc_count - allele1_count)/max(anc_count, allele1_count)
        allele2_ratio_score = abs(anc_count - allele2_count)/max(anc_count, allele2_count)
        allele3_ratio_score = abs(anc_count - allele3_count)/max(anc_count, allele3_count)
        
        #merging information for individual dataframes
        df1 = data.frame(anc_allele_count = anc_count, allele_count = allele1_count, ratio_score = allele1_ratio_score, name = pos$name)
        pwm_allele1_df_ib = rbind(pwm_allele1_df_ib, df1)
        
        df2 = data.frame(anc_allele_count = anc_count, allele_count = allele2_count, ratio_score = allele2_ratio_score, name = pos$name)
        pwm_allele2_df_ib = rbind(pwm_allele2_df_ib, df2)
        
        df3 = data.frame(anc_allele_count = anc_count, allele_count = allele3_count, ratio_score = allele3_ratio_score, name = pos$name)
        pwm_allele3_df_ib = rbind(pwm_allele3_df_ib, df3)
        
      }
    }
  }
  
  #@print check
  print(paste("(Part 3) For IB population - Gone over row : ",row, " of total rows : ",dim(mutation_annotate_ib)[1],sep = ""))
  
}

#binding all the allele outputs in one table
allele_lengths_pwm_ib = rbind(pwm_allele1_df_ib, pwm_allele2_df_ib, pwm_allele3_df_ib)


####### PWM lengths for NS ########

#declaring empty dataframes to fill in allele-specific information - please note the alleles here are in no particular order of preference just that pairs are maintained
pwm_allele1_df_ns = data.frame()
pwm_allele2_df_ns = data.frame()
pwm_allele3_df_ns = data.frame()

#going over every PWM and their respective positions within the alignments
for(row in 1:nrow(mutation_annotate_ns)){
  
  #extracting the var information
  var = mutation_annotate_ns[row,]
  
  ##### Step 1 - Extracting the CRM sequence info and the TFBS sequence in question within this CRM #####
  
  #extracting the coordinates of the CRM
  crm_chr = var$chr
  crm_start = var$start_crm
  crm_end = var$stop_crm
  
  #what the crm sequence filename will be
  crm_filename = paste("crm_",crm_chr,"_",crm_start,"_",crm_end,".fa", sep = "")
  
  #giving a call to the function to extract the outgroup sequence directly from the sequence files
  out_seq = extract_outgroup_sequence(crm_filename)
  
  #extracting the sequence where the TFBS is predicted to be binding to now
  ## start of the TFBS ==> (var$Position - Var$PosInTFBS + 1) ------ var$Position - var$PosInTFBS tells me where the first TFBS position is, adding 1 converts bed0 to bed1
  ## end of TFBS ==> start of TFBS + length of the PWM
  pwm_seq = substr(x = out_seq, start = var$Position - var$PosInTFBS + 1, stop = (var$Position - var$PosInTFBS + 1) + dim(tfbs_counts_info[tfbs_counts_info$name == var$PWM,])[1] - 1)
  
  ### the problem of reverse complementing DNA sequences 
  ### if variants are found on the reverse strand (which is recognized by "reverse == TRUE") then the whole PWM is by default recognized on the reverse strand
  ### in such cases I will also have to reverse complement the outgroup seqeuence to obtain the ancestral allele for the reverse complemented positions
  if(var$Reverse == "true"){
    
    #if you are inside it means that the pwm is identified on the reverse complemented strand - so reverse complementing the ancestral sequence
    pwm_seq = as.character(reverseComplement(DNAString(pwm_seq)))
    
  } #if Reverse == FALSE then continue without deroute
  
  #extracting the tfbs for the pwm
  pwm_info = tfbs_counts_info[tfbs_counts_info$name == var$PWM, c("A", "C", "G", "T", "name")]
  rownames(pwm_info) = NULL
  
  ###### In some weird cases for drosophila (possibly for other species too) the outgroup sequences is terminated in between for certain PWMs
  ###### controlling for that by checking if the length of the outgroup PWM sequence is similar to the pwm length
  if(dim(pwm_info)[1] == nchar(pwm_seq)){
    
    #adding the ancestral state information here
    pwm_info$anc_allele = unlist(strsplit(pwm_seq, ""))
    
    ##### Step 2 - Calculating the pairwise entropies per position per PWM - every position will have 3 entropies #####
    
    #first checking and removing any pwm positions that have a gap in the ancestral allele
    pwm_info = pwm_info[!pwm_info$anc_allele == "-",]
    
    #if there is no information left then bypassing this pwm
    if(dim(pwm_info)[1] > 0){
      
      #going over every position in a for loop
      for(pos in 1:nrow(pwm_info)){
        
        #extracting the pos information
        pos = pwm_info[pos,]
        
        #listing all the other alleles that are not ancestral alleles
        var_alleles = colnames(pos)[1:4][!colnames(pos)[1:4] == pos$anc_allele]
        
        #storing counts of all alleles and adding 1 to them
        anc_count = as.numeric(pos[colnames(pos) == pos$anc_allele]) + 1
        allele1_count = as.numeric(pos[colnames(pos) == var_alleles[1]]) + 1
        allele2_count = as.numeric(pos[colnames(pos) == var_alleles[2]]) + 1
        allele3_count = as.numeric(pos[colnames(pos) == var_alleles[3]]) + 1
        
        #calculating the pairwise entropy - here allele numbering is in order of A,C,G and T - not that it matters
        allele1_ratio_score = abs(anc_count - allele1_count)/max(anc_count, allele1_count)
        allele2_ratio_score = abs(anc_count - allele2_count)/max(anc_count, allele2_count)
        allele3_ratio_score = abs(anc_count - allele3_count)/max(anc_count, allele3_count)
        
        #merging information for individual dataframes
        df1 = data.frame(anc_allele_count = anc_count, allele_count = allele1_count, ratio_score = allele1_ratio_score, name = pos$name)
        pwm_allele1_df_ns = rbind(pwm_allele1_df_ns, df1)
        
        df2 = data.frame(anc_allele_count = anc_count, allele_count = allele2_count, ratio_score = allele2_ratio_score, name = pos$name)
        pwm_allele2_df_ns = rbind(pwm_allele2_df_ns, df2)
        
        df3 = data.frame(anc_allele_count = anc_count, allele_count = allele3_count, ratio_score = allele3_ratio_score, name = pos$name)
        pwm_allele3_df_ns = rbind(pwm_allele3_df_ns, df3)
        

      }
    }
  }
  #@print check
  print(paste("(Part 4) For NS population PWM - Gone over row : ",row, " of total rows : ",dim(mutation_annotate_ns)[1],sep = ""))
  
}

#binding all the allele outputs in one table
allele_lengths_pwm_ns = rbind(pwm_allele1_df_ns, pwm_allele2_df_ns, pwm_allele3_df_ns)

#second all the positions within the PWM and their entropies
write.table(x = allele_lengths_pwm_ib, file = "ratios_info/pwm_position_ratioscore_info_aratha_ib.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = allele_lengths_pwm_ns, file = "ratios_info/pwm_position_ratioscore_info_aratha_ns.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
