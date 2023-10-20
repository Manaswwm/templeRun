#collecting outputs from TEMPLE analysis (aggregated) and splitting the variants for the two populations

##############################################
#! -------- declaring paths -------- !#
aggregated_output_path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBS_project_new/aratha/aggregate_outputs/mutation_output_aratha_8508_crms.txt"

##############################################

#importing relevant libraries
library(stringr)
library(ggplot2)
library(dplyr)

#importing mutation output file directly here
mutation_output_subset = read.delim(file = aggregated_output_path, header = TRUE, sep = "\t")

### In most cases I have a single variant affecting the binding property of multiple TFs - so single position variant in a CRM affecting multiple TFBSs ###
#first - creating a new single column that will give me the positional id of CRM that is unique to every CRM as it will contain chromosome and start and stop positons
mutation_output_subset$crm_id = paste(mutation_output_subset$chr,":",mutation_output_subset$start_crm,"-",mutation_output_subset$stop_crm, sep = "")

#changing data-type of effect on nucleotide scores to subset unique mutations
mutation_output_subset$EffectOnNucleotide = as.numeric(mutation_output_subset$EffectOnNucleotide)

#sanity check - removing any variants that have NA for effect on nucleotide score
mutation_output_subset = mutation_output_subset[!is.na(mutation_output_subset$EffectOnNucleotide),]

#### trial where I take all point mutations and not only one representative ones that has the highest effect on nucleotide score
#### TEMPLE reports point mutations but a single point mutation could affect multiple TFBS with different intensity
#### by only keeping mutations that have the highest effect on nucleotide score (absolute score) I lose the potentially neutral (syn) variants by a lot and I lose a lot of variants in general

#seperating dataframes for both populations to contain no missing data
ib_mutation_subset = mutation_output_subset[, c("PWM", "Position", "PosInTFBS", "Reverse", "Polarization", "Al1", "Al2", "Al1_pop_IB", "Al2_pop_IB", "NumbMisLines_pop_IB", "Pi", "EffectOnNucleotide", "chr", "start_crm", "stop_crm")]
ib_mutation_subset = ib_mutation_subset[ib_mutation_subset$NumbMisLines_pop_IB == 0,]
ib_mutation_subset = ib_mutation_subset[ib_mutation_subset$Polarization == 1,]

ns_mutation_subset = mutation_output_subset[, c("PWM", "Position", "PosInTFBS", "Reverse", "Polarization", "Al1", "Al2", "Al1_pop_NS", "Al2_pop_NS", "NumbMisLines_pop_NS", "Pi", "EffectOnNucleotide", "chr", "start_crm", "stop_crm")]
ns_mutation_subset = ns_mutation_subset[ns_mutation_subset$NumbMisLines_pop_NS == 0,]
ns_mutation_subset = ns_mutation_subset[ns_mutation_subset$Polarization == 1,]

#removing fixed or lost alleles
ib_mutation_subset = ib_mutation_subset[!ib_mutation_subset$Al1_pop_IB %in% c(0,45),]

ns_mutation_subset = ns_mutation_subset[!ns_mutation_subset$Al1_pop_NS %in% c(0,45),]

#constructing frequencies of alleles for sfs
ib_mutation_subset$Al1_pop_IB = as.numeric(ib_mutation_subset$Al1_pop_IB)
ib_mutation_subset$Al2_pop_IB = as.numeric(ib_mutation_subset$Al2_pop_IB)
ns_mutation_subset$Al1_pop_NS = as.numeric(ns_mutation_subset$Al1_pop_NS)
ns_mutation_subset$Al2_pop_NS = as.numeric(ns_mutation_subset$Al2_pop_NS)

ib_mutation_subset$Al1_pop_IB = ib_mutation_subset$Al1_pop_IB/45
ib_mutation_subset$Al2_pop_IB = ib_mutation_subset$Al2_pop_IB/45
ns_mutation_subset$Al1_pop_NS = ns_mutation_subset$Al1_pop_NS/45
ns_mutation_subset$Al2_pop_NS = ns_mutation_subset$Al2_pop_NS/45

### writing the pop-specific mutations to a file now for calculating lengths and then pin/pis ratios
write.table(x = ib_mutation_subset, file = "pop_specific_variants/ib_pop_mutation_subset.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = ns_mutation_subset, file = "pop_specific_variants/ns_pop_mutation_subset.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

#plotting the p-value distribution
ggplot(data = ib_mutation_subset) + geom_histogram(aes(x = EffectOnNucleotide), color = "black", fill = "white", binwidth = 1) + theme_classic() +
  xlab("Difference in score p-value")
ggsave("pvalue_distribution_8414crms.png", device = "png", height = 5, width = 8)




