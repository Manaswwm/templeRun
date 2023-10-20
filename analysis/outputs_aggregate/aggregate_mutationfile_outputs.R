####### In case of a large number of output files - please consider running on a high performance cluster ######

#######################################################
#! ----- declaring paths ----- !#
output_directory = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project_new/aratha/temple_batches_nopval/batch_1/"

#######################################################  
  
#script that aggregates output files from TEMPLE and creates a combined file that could furthur processed
#importing relevant libraries
library(stringr)

#listing all folders and entering their subfolders for all batches
output_filenames = list.files(path = output_directory, pattern = "crm_*/*", full.names = TRUE)
output_filenames = list.files(path = paste(output_filenames, "/", sep = ""), full.names = TRUE)

### done listing all output folder names - turning to specific files and listing them
mutation_filenames = list.files(path = paste(output_filenames,"/", sep = ""), pattern = "*/*_MUTATION.csv", full.names = TRUE)

#going over all mutation files in a sequential manner and storing the data in a dataframe - using for loop

#delcaring the empty dataframe
mutation_output = data.frame()

#inserting for loop here
for(item in mutation_filenames){
  
  #opening a single mutation file - have to skip 57 lines as they contain metadata
  #have to modify the following code because of the absence of the last column
  df = read.delim(file = item, skip = 103, sep = ",", header = FALSE)
  
  #removing last column - this is NULL as we do not have the difference in p value scores
  df = df[,1:(ncol(df)-1)]
  
  #giving back column names
  colnames(df) = df[1,]
  
  #removing first entry as these are the column names themselves
  df = df[-1,]
  rownames(df) = NULL
  
  #checking if the mutation output file is non-empty
  if(dim(df)[1] > 0){
    
    #giving back the chromosome and start and stop position of the CRM
    df$chr = strsplit(str_match(item, "crm_output_\\s*(.*?)\\s*//")[2], "_")[[1]][1]
    df$start_crm = strsplit(str_match(item, "crm_output_\\s*(.*?)\\s*//")[2], "_")[[1]][2]
    df$stop_crm = strsplit(str_match(item, "crm_output_\\s*(.*?)\\s*//")[2], "_")[[1]][3]
    
    #merging with the original dataframe
    mutation_output = rbind(mutation_output, df)
    
    #@print check
    print(paste("Gone over file numer : ", match(item, mutation_filenames), " of ",length(mutation_filenames), sep = ""))
  }
  
}

#writing to a file
write.table(x = mutation_output, file = paste("mutation_output_aratha_",length(batch_mutation_filenames),"_crms.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
