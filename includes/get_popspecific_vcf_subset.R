#function that will take in the pops-specific accessions and crm coordinates to constuct filtered pop-specific vcf file
get_popspecific_vcf_subset = function(ib_pop, ns_pop, crm_coord){

  ## for A. thaliana the vcf file is directly accessed from the 1001 genomes ##
  
  # #constructing the url for REST API
  # #the main url
  # rest_command_main = "https://tools.1001genomes.org/api/v1/vcfsubset/"
  # 
  # #strains per population
  # #IB
  # rest_command_pop_ib = paste("strains/", paste(ib_pop$V1, collapse = ","), "/", sep = "")
  # #NS
  # rest_command_pop_ns = paste("strains/", paste(ns_pop$V1, collapse = ","), "/", sep = "")
  # 
  # #region of interest
  # rest_command_region = paste("regions/Chr", crm_coord$seqnames, ":", crm_coord$start, "..", crm_coord$end, "/", sep = "")
  # 
  # #meta information
  # rest_command_meta = "type/fullgenome/format/vcf"
  # 
  
  #### Maria has 1001G vcf file downloaded on netscratch so using it locally ####
  #declaring constants in the command
  bcftools_main = paste("bcftools view /netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project_new/aratha/1001g_vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz ")
  bcftools_region = paste(" -r ",crm_coord$seqname, ":", crm_coord$start, "-", crm_coord$end, sep = "")
  
  #strains per population
  #IB
  bcftools_accessions_ib = paste(" -s ", paste(ib_pop$V1, collapse = ","), sep = "")
  #NS
  bcftools_accessions_ns = paste(" -s ", paste(ns_pop$V1, collapse = ","), sep = "")
  
  ### first for IB
  #patching together the URL for IB and issuing the command
  request_ib = system(paste(bcftools_main, bcftools_accessions_ib, bcftools_region, sep = ""), intern = TRUE)
  
  #checking if 1001 genomes has given me any output - status code should be 200
  if(length(request_ib) > 18){
    
    #post-processing the output
    #bcftools_command_ceu = bcftools_command_ceu[-grep(pattern = "^##", x = bcftools_command_ceu)]
    # bcftools_command_ceu[1] = str_remove(string = bcftools_command_ceu[1], pattern = "#")
    # bcftools_command_ceu = read.table(text = bcftools_command_ceu, sep="\t", header = TRUE)
    # 
    # #changing the first column name
    # colnames(bcftools_command_ceu)[1] = "#CHROM"
    
    #declaring the name of the file to write the raw vcf into
    path_to_raw_vcf_ib = paste("backup_files/vcf_files/ib_", crm_coord$seqnames,"_",crm_coord$start,"_",crm_coord$end,"_raw.vcf", sep = "")

    #writing to a file now
    writeLines(request_ib, path_to_raw_vcf_ib)
    
    # #writing to the file now
    # cat(content(request_ib), file = path_to_raw_vcf_ib)
    
    ##filtering the vcf file now##
    #delcaring the name of the file to write the filtered info
    path_to_filtered_vcf_ib = paste("backup_files/vcf_files/ib_", crm_coord$seqnames,"_",crm_coord$start,"_",crm_coord$end,"_filtered.vcf", sep = "")
    
    #filtering vcf file and writing the filtered file
    filtering_command = paste("vcftools --vcf ",path_to_raw_vcf_ib, " --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf ",maf," --recode --stdout > ",path_to_filtered_vcf_ib, sep="")
    
    #issuing command
    system(filtering_command)
    
    # ## re-adjusting the file to insert the vcf format - required for vcfkit to work ##
    # #reading the vcf filtered file here
    # filtered_vcf_ceu = readLines(path_to_filtered_vcf_ceu)
    # 
    # #reading the file here as a list - adding the file format in the first position without disturbing the numbering -so 1st element will be 2nd, 2nd is 3rd ..
    # filtered_vcf_ceu = append(filtered_vcf_ceu, values = "##fileformat=VCFv4.1", after = 0) #after 0 to insert this as the first element
    # 
    # #writing to a the same vcf file
    # writeLines(text = filtered_vcf_ceu, path_to_filtered_vcf_ceu)
  }
  
  ### second for NS
  #patching together the URL for IB and issuing the command
  request_ns = system(paste(bcftools_main, bcftools_accessions_ns, bcftools_region, sep = ""), intern = TRUE)
  
  #checking if 1001 genomes has given me any output - status code should be 200
  if(length(request_ns) > 18){
    
    #post-processing the output
    #bcftools_command_ceu = bcftools_command_ceu[-grep(pattern = "^##", x = bcftools_command_ceu)]
    # bcftools_command_ceu[1] = str_remove(string = bcftools_command_ceu[1], pattern = "#")
    # bcftools_command_ceu = read.table(text = bcftools_command_ceu, sep="\t", header = TRUE)
    # 
    # #changing the first column name
    # colnames(bcftools_command_ceu)[1] = "#CHROM"
    
    # #declaring the name of the file to write the raw vcf into
    path_to_raw_vcf_ns = paste("backup_files/vcf_files/ns_", crm_coord$seqnames,"_",crm_coord$start,"_",crm_coord$end,"_raw.vcf", sep = "")

    #writing to a file now
    writeLines(request_ns, path_to_raw_vcf_ns)
    
    # #writing to the file now
    # cat(content(request_ns), file = path_to_raw_vcf_ns)

    ##filtering the vcf file now##
    #delcaring the name of the file to write the filtered info
    path_to_filtered_vcf_ns = paste("backup_files/vcf_files/ns_", crm_coord$seqnames,"_",crm_coord$start,"_",crm_coord$end,"_filtered.vcf", sep = "")
    
    #filtering vcf file and writing the filtered file
    filtering_command = paste("vcftools --vcf ",path_to_raw_vcf_ns, " --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf ",maf," --recode --stdout > ",path_to_filtered_vcf_ns, sep="")
    
    #issuing command
    system(filtering_command)
    
    # ## re-adjusting the file to insert the vcf format - required for vcfkit to work ##
    # #reading the vcf filtered file here
    # filtered_vcf_ceu = readLines(path_to_filtered_vcf_ceu)
    # 
    # #reading the file here as a list - adding the file format in the first position without disturbing the numbering -so 1st element will be 2nd, 2nd is 3rd ..
    # filtered_vcf_ceu = append(filtered_vcf_ceu, values = "##fileformat=VCFv4.1", after = 0) #after 0 to insert this as the first element
    # 
    # #writing to a the same vcf file
    # writeLines(text = filtered_vcf_ceu, path_to_filtered_vcf_ceu)
  }
  
}