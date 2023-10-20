## templeRun - this is a wrapper around TEMPLE - here we execute TEMPLE across multiple noncoding regions of interest by constructing the input files internally  

###################################################
#! ----- declaring paths here----- !#
#loading populations - both pops have 45 individuals
ib_pop = read.delim("pop_file/IB_pop_45inds.txt", header = FALSE)
ns_pop = read.delim("pop_file/NS_pop_45inds.txt", header = FALSE)

#listing all the CRMs for which I have to extract the sequence
crm_info = read.delim("input_files/crm_30_2kb_upstream_noCDS_wPWMinfo.txt", header = TRUE, sep = "\t")

#pwm file that will be used for TEMPLE - this pwm file comes with a p-value calculation from TEMPLE
pwm_file = "input_files/processed_PWMforTEMPLE_withCRMinfo.txt"

#path to the ingroup reference sequence
ingroup_ref = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project/misc/temple_analysis_perind/aratha/input_files/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"

###################################################

#importing relevant libraries
library(httr)
library(jsonlite)
library(xml2)
library(Biostrings)
library(stringr)
library(rlist)

#sourcing relevant functions
source("includes/misc.R")
source("includes/get_popspecific_vcf_subset.R")
source("includes/make_sequences_for_temple.R")
source("includes/get_outgroup_sequence.R")
source("includes/make_regionfile_for_temple.R")

#declaring constants
maf = 0.0001

#reordering the positions
crm_info = crm_info[order(crm_info$seqnames, crm_info$start),]
rownames(crm_info) = NULL

#declaring an empty dataframe that will collect information on the outlier crms
crm_outlier_df = data.frame()

#going over all the crms sequentially and do the following
#1) making population and region specific vcf files and filtering
#2) constructing individual specific and outgroup sequences with variant information
#3) constructing region file
#4) executing TEMPLE
#5) storing information on outlier CRMs - if there are any -- outlier crms are regions that have >0.8 gaps when aligned with outgroup species

#using for-loop to go over every file
for(crm in 1:nrow(crm_info)){
  
  #single crm region extraction
  crm_coord = crm_info[crm,]

  #### Step 1 - making population and region specific vcf files, and filtering ####
  #giving call to a pre-set function
  get_popspecific_vcf_subset(ib_pop = ib_pop, ns_pop = ns_pop, crm_coord = crm_coord)
  
  #### Step 2 - constructing sequences of the variant positions - if there are any ####
  #important deroute - I have to give paths for the filtered vcf files for both pops to the next function for creating individual specific sequences
  path_to_filtered_vcf_ib = paste("backup_files/vcf_files/ib_", crm_coord$seqnames,"_",crm_coord$start,"_",crm_coord$end,"_filtered.vcf", sep = "")
  path_to_filtered_vcf_ns = paste("backup_files/vcf_files/ns_", crm_coord$seqnames,"_",crm_coord$start,"_",crm_coord$end,"_filtered.vcf", sep = "")
  
  #giving the call for generating sequences (individual-specific and outgroup sequence)
  #the return to the following function is only triggered by outlier crm coordinates (more than 0.8 gaps in alignment with outgroup)
  #if alignments are okay, the return is an empty dataframe - the sequence construction is performed by-default and the return does not
  #serve as an indicator for sequence construction
  crm_outlier = make_sequences_for_temple(path_to_filtered_vcf_ib = path_to_filtered_vcf_ib, path_to_filtered_vcf_ns = path_to_filtered_vcf_ns, crm_coord = crm_coord)
        
  ##### Step 3 - constructing region file #####
  make_regionfile_for_temple(crm_coord = crm_coord)
    
  ##### Step 4 - Running TEMPLE now #####
  #delcaring the location for the sequence file
  seq_file = paste("backup_files/sequence_files/crm_", crm_coord$seqnames,"_",crm_coord$start, "_", crm_coord$end, ".fa", sep = "")
  
  #declaring the location for the region file
  reg_file = paste("backup_files/region_files/crm_",crm_coord$seqnames,"_",crm_coord$start,"_",crm_coord$end,"_regionfile.txt", sep = "")
  
  ## important point in TEMPLE - checking if both sequence and region files are constructed ##
  if(file.exists(seq_file) & file.exists(reg_file)){
    
    #adding stuff together
    #1 1 2 --> number of threads, threshold for missing data, analysis mode ---- for two population analysis the mode is 2
    temple_command_construct = paste("java -jar gTEMPLE_v1.0.jar ", pwm_file, " ", seq_file, " ", reg_file, " 1 1 2", sep = "")
    
    #launching command
    system(temple_command_construct)
    
    #everytime temple makes the file with the same name, once the file is created I will have to rename the file to make space for the next file
    file.rename(from = "pop_IB_pop_NS/", to = paste("crm_output_",crm_coord$seqnames,"_",crm_coord$start,"_",crm_coord$end, "/", sep = ""))
  
  }
  
  ##### Step 5 - checking for outlier crms from Step 2 - if current crm is outlier then collecting info #####
  #checking if crm is an outlier by checking if the dataframe is empty from Step 2
  if(dim(crm_outlier)[1] > 0){
    
    #merging with the main dataframe
    crm_outlier_df = rbind(crm_outlier_df, crm_outlier)
    
  }

  #@print check
  print(paste("Gone over CRM : ", crm_coord$seqnames, "_", crm_coord$start, "_", crm_coord$end, sep = ""))

}#gone over a single CRMs

#writing the outlier CRMs to a file
write.table(x = crm_outlier_df, file = "outlier_crms_coordinates.txt", row.names = FALSE, col.names = TRUE, sep = "\t")