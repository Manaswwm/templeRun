#this script will take in the information on the ratio scores and calculate the syn and nonsyn ratios for populations
#importing relevant packages
library(ggplot2)
library(dplyr)

#importing the pop-specific ratios score and the PWM ratios score
pos_ratios_ib = read.delim(file = "ratios_info/pwm_position_ratioscore_info_aratha_ib.txt", header = TRUE, sep = "\t")
ib_pop_ratios = read.delim(file = "ratios_info/mutation_annotate_info_ib.txt", header = TRUE, sep = "\t")

pos_ratios_ns = read.delim(file = "ratios_info/pwm_position_ratioscore_info_aratha_ns.txt", header = TRUE, sep = "\t")
ns_pop_ratios = read.delim(file = "ratios_info/mutation_annotate_info_ns.txt", header = TRUE, sep = "\t")

#for every mutation in both pops - assigning a unique position id
ib_pop_ratios$pos = paste(ib_pop_ratios$chr, "_",ib_pop_ratios$start_crm + ib_pop_ratios$Position - 1, sep = "")
ns_pop_ratios$pos = paste(ns_pop_ratios$chr, "_",ns_pop_ratios$start_crm + ns_pop_ratios$Position - 1, sep = "")

#keeping only unique point mutations in both populations
ib_pop_ratios = ib_pop_ratios[!duplicated(ib_pop_ratios$pos),]
ns_pop_ratios = ns_pop_ratios[!duplicated(ns_pop_ratios$pos),]

##### assinging nonsyn and syn for both populations and the backgounrd pwm positions #####
#### population IB ####

##taking coding stats which will be used later
### Part 1 - listing output files for both regions
cds_stats_files_ib = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_aratha_IB/", pattern = "poly_div_gene_*", full.names = TRUE)

### Part 2 - aggregating all batch files
poly_div_stats_cds_ib = lapply(cds_stats_files_ib, function(x){read.delim(x, sep = "\t")})
poly_div_stats_cds_ib = do.call(rbind, poly_div_stats_cds_ib)
poly_div_stats_cds_ib = unique(poly_div_stats_cds_ib) ## -- taking unique due to the nature of DNABD naming file - scope for improvement

####### taking syn length within the noncoding is inflating the estimates of pi_n/pi_s because overall there is an excess of syn sites in TFBS as compared to the coding regions in terms of proportion to nonsyn
## an alternative approach is to just use the pi_s from the coding regions

#taking the max count in the pair of counts
ib_pop_ratios$max_count = by(ib_pop_ratios, seq_len(nrow(ib_pop_ratios)), function(x) {max(x$ref_allele_count, x$alt_allele_count)})
ib_pop_ratios$max_count = as.numeric(ib_pop_ratios$max_count)

#seperating syn and nonsyn
ib_nonsyn = ib_pop_ratios[ib_pop_ratios$ratio_score >= 0.6 & ib_pop_ratios$max_count >= 400,]
ib_syn = ib_pop_ratios[!ib_pop_ratios$pos %in% ib_nonsyn$pos,]

#next - going to the background pwm positions and identifying the nonsyn and syn lengths
#first identifying the max count alleles and listing them
pos_ratios_ib$max_count = by(pos_ratios_ib, seq_len(nrow(pos_ratios_ib)), function(x){max(x$anc_allele_count, x$allele_count)})
pos_ratios_ib$max_count = as.numeric(pos_ratios_ib$max_count)

#checking the dimensions for the nonsyn criteria
nonsyn_length_ib = dim(pos_ratios_ib[pos_ratios_ib$ratio_score >= 0.6 & pos_ratios_ib$max_count >= 400,])[1]
syn_length_ib = dim(pos_ratios_ib)[1] - nonsyn_length_ib

##making ratios now
#constructing pi_n and pi_s for noncoding
pi_n_ib = sum(ib_nonsyn$Pi)/nonsyn_length_ib
pi_s_ib = sum(ib_syn$Pi)/syn_length_ib

#pi_n/pi_s with noncoding pi_s
pi_n_ib/pi_s_ib

#pi_n/pi_s with coding pi_s -- coding pi_s stays constant
pi_n_ib/pi_syn_coding

#### preparing for asymptotic MK by taking nonsyn from noncoding and syn from coding - the sfs will only consist of 10 classes to normalize ####

## importing data for the coding regions
#listing backup files
ib_backup_files = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_aratha_IB/", pattern = "backup_*", full.names = TRUE)

#declaring an empty dataframe to fill in the information
freq_cds_ib = data.frame()

#going over all batch files sequentially
for(file in ib_backup_files){
  
  #reading a single batch file
  load(paste(file,"/frequencies_cds", sep = ""))
  
  #merging
  freq_cds_ib = rbind(freq_cds_ib, frequencies_cds)
  
  
}

#cleaning up
rm(file, ib_backup_files, frequencies_cds)

#taking unique
freq_cds_ib = unique(freq_cds_ib)

#merging
syn_freq_ib = data.frame(freq = freq_cds_ib$freq_der[freq_cds_ib$type == "syn"])
syn_freq_ib$rec_level = cut(syn_freq_ib$freq, breaks = seq(0,1,by=0.05))
syn_freq_ib = syn_freq_ib %>%
  rowwise() %>%
  group_by(rec_level) %>%
  tally(name = "n") %>%
  mutate(prop = n/sum(n))


## here we use Al2 as this is the derived allele and Al1 is the ancestral allele
nonsyn_freq_ib = data.frame(freq = ib_nonsyn$Al2_pop_IB)
nonsyn_freq_ib$rec_level = cut(nonsyn_freq_ib$freq, breaks = seq(0,1,by=0.05))
nonsyn_freq_ib = nonsyn_freq_ib %>%
  rowwise() %>%
  group_by(rec_level) %>%
  tally(name = "n") %>%
  mutate(prop = n/sum(n))

#merging information together
syn_nonsyn_freq_ib = data.frame(X = c(seq(0.05,0.95, by = 0.05), 0.99) , pN = nonsyn_freq_ib$n, pS = syn_freq_ib$n)

#writing this to a file now
write.table(x = syn_nonsyn_freq_ib, file = "asympMK_analysis/ib_freq_noncoding.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#cleaning up
rm(syn_freq_ib, nonsyn_freq_ib, syn_nonsyn_freq_ib, ib_pop_ratios, ib_syn)

#### calculating class specific alpha locally for noncoding nonsyn and coding syn #####

##declaring constants first - some of these will be taken directly 
div_syn_coding_ib = 0.13728
div_nonsyn_noncoding_ib = 0.02950
nonsyn_length_ib = nonsyn_length_ib
syn_length_ib = sum(poly_div_stats_cds_ib$syn_length)

#declaring an empty dataframe to collect alpha values
noncoding_alpha_ib = data.frame()

#subsetting
freq_cds_syn_ib = freq_cds_ib[freq_cds_ib$type == "syn",]

#going over interval classes and calculating alpha
for(int in seq(0, 1, by = 0.05)){
  
  #extracting frequencies belonging to these classes
  syn_freq = freq_cds_ib$freq_der[freq_cds_ib$freq_der >= int & freq_cds_ib$freq_der < (int+0.05)]
  nonsyn_freq = ib_nonsyn$Al2_pop_IB[ib_nonsyn$Al2_pop_IB >= int & ib_nonsyn$Al2_pop_IB < (int+0.05)]
  
  #calculating pi
  pi_s = sum(2*syn_freq*(1-syn_freq))/syn_length_ib
  pi_n = sum(2*nonsyn_freq*(1-nonsyn_freq))/nonsyn_length_ib
  
  #calculating alpha
  alpha = 1-((div_syn_coding_ib/div_nonsyn_noncoding_ib)) * (pi_n/pi_s)
  
  #binding together
  df = data.frame(alpha = alpha, int = int)
  
  #merging
  noncoding_alpha_ib = rbind(noncoding_alpha_ib, df)
  
}

#removing last class and renaming first class
noncoding_alpha_ib = noncoding_alpha_ib[1:(dim(noncoding_alpha_ib)[1]-1),]
noncoding_alpha_ib$int[1] = 0.01

#plotting estimates
ggplot(data = noncoding_alpha_ib, aes(x = int, y = alpha)) + geom_point(size = 2) + theme_bw() + xlab("SFS classes") + scale_y_continuous(breaks = seq(-2,2, by = 0.2)) + ggtitle("A. thaliana - IB CRM alpha estimates per SFS class")
ggsave(filename = "asympMK_analysis/class_estimates/ib_alpha_sfs_perclass.png", device = "png", dpi = 300)

#### population NS ####

####### taking syn length within the noncoding is inflating the estimates of pi_n/pi_s because overall there is an excess of syn sites in TFBS as compared to the coding regions in terms of proportion to nonsyn
## an alternative approach is to just use the pi_s from the coding regions

#taking the max count in the pair of counts
ns_pop_ratios$max_count = by(ns_pop_ratios, seq_len(nrow(ns_pop_ratios)), function(x) {max(x$ref_allele_count, x$alt_allele_count)})
ns_pop_ratios$max_count = as.numeric(ns_pop_ratios$max_count)

#seperating syn and nonsyn
ns_nonsyn = ns_pop_ratios[ns_pop_ratios$ratio_score >= 0.6 & ns_pop_ratios$max_count >= 400,]
ns_syn = ns_pop_ratios[!ns_pop_ratios$pos %in% ns_nonsyn$pos,]

#next - going to the background pwm positions and identifying the nonsyn and syn lengths
#first identifying the max count alleles and listing them
pos_ratios_ns$max_count = by(pos_ratios_ns, seq_len(nrow(pos_ratios_ns)), function(x){max(x$anc_allele_count, x$allele_count)})
pos_ratios_ns$max_count = as.numeric(pos_ratios_ns$max_count)

#checking the dimensions for the nonsyn criteria
nonsyn_length_ns = dim(pos_ratios_ns[pos_ratios_ns$ratio_score >= 0.6 & pos_ratios_ns$max_count >= 400,])[1]
syn_length_ns = dim(pos_ratios_ns)[1] - nonsyn_length_ns

##making ratios now
#constructing pi_n and pi_s for noncoding
pi_n_ns = sum(ns_nonsyn$Pi)/nonsyn_length_ns
pi_s_ns = sum(ns_syn$Pi)/syn_length_ns

#pi_n/pi_s with noncoding pi_s
pi_n_ns/pi_s_ns

#pi_n/pi_s with coding pi_s -- coding pi_s stays constant
pi_n_ns/pi_syn_coding

#### the hypothesis for pi_n/pi_s being larger than 1 is that syn and nonsyn length difference is not as high as it is usually in the coding region
#### for coding regions nonsyn_length/syn_length  ~ 3; however so far in the coding regions this is : nonsyn_length/syn_length ~1.2
#checking with PWMs are contributing to higher syn length and then confirming if this signal is true
length_check = aratha_pos_ratios
length_check$type = "syn"
length_check$type[length_check$ratio_score > 0.6 & length_check$max_count > 400] = "nonsyn"

#checking the top 10 contributors to the syn_length
sort(table(length_check$name[length_check$type == "syn"]), decreasing = TRUE)[1:10]


##### making the sfs ratios #####
ib_sfs = data.frame(freq_der = c(ib_nonsyn$Al2_pop_IB, ib_syn$Al2_pop_IB), type = rep(c("nonsyn", "syn"), c(781, 760)))

#processing to make it compatible with sfs construction
ib_sfs$rec_level = cut(ib_sfs$freq_der, breaks = seq(0,1,by=0.1))
ib_sfs = ib_sfs %>%
  rowwise() %>%
  group_by(type, rec_level) %>%
  tally(name = "n") %>%
  mutate(prop = n/sum(n))

#adding another column to denote frequencies
ib_sfs$freq = rep(seq(0.1,1, by = 0.1), 2)

#making the sfs plots -- unique variants
ggplot(ib_sfs, aes(x = freq, y = prop, fill = type)) + geom_bar(stat="identity", position = position_dodge(), colour = "black", alpha = 0.6) + theme_classic()+
  xlab("Frequencies") + ylab("Proportion") + labs(fill = "Population") + scale_x_continuous(breaks = seq(0,1,by=0.1))+
  ggtitle("Iberian population SFS - nonsyn and syn variants") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "output_plots/sfs_NS_uniquevariants_5485CRMs.png", device = "png", width = 8, height = 9)