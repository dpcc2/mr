#### PQTL #### 

## Load packages ##### 

library(dplyr)
library(TwoSampleMR)
library(readxl)


###### Import datasets ####

pqtl_exposure <- read_excel("Documents/eqtl.pqtl.ocac.project/pqtl/ncomms14357-s6.xlsx", sheet = "Supplemental Data 5 ")

ocac_outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr22.txt") # Change chromosome number


###### Manipulate pqtl dataset ###### 

names(pqtl_exposure)[names(pqtl_exposure) == "TargetFullName"] <- "protein"
names(pqtl_exposure)[names(pqtl_exposure) == "Beta"] <- "beta.exposure"
names(pqtl_exposure)[names(pqtl_exposure) == "P-value"] <- "pval.exposure"
names(pqtl_exposure)[names(pqtl_exposure) == "Chr"] <- "chr"
names(pqtl_exposure)[names(pqtl_exposure) == "Position"] <- "pos"
names(pqtl_exposure)[names(pqtl_exposure) == "Stat"] <- "stat.exposure"
names(pqtl_exposure)[names(pqtl_exposure) == "Allele"] <- "effect_allele.exposure"
pqtl_exposure$se.exposure <- pqtl_exposure$beta.exposure/pqtl_exposure$stat.exposure
pqtl_exposure <- pqtl_exposure[ , c(1:6, 11, 7, 8, 10)]


###### Manipulate OCAC dataset ###### 

outcome <- ocac_outcome[,c(2:6, 49:52)] # Change selected columns for different subtypes
names(outcome)[names(outcome) == "X1000G_SNPname"] <- "SNP"
names(outcome)[names(outcome) == "Chromosome"] <- "chr"
names(outcome)[names(outcome) == "Position"] <- "pos"
names(outcome)[names(outcome) == "Effect"] <- "effect_allele.outcome"
names(outcome)[names(outcome) == "Baseline"] <- "other_allele.outcome"
names(outcome)[6] <- "beta.outcome"
names(outcome)[7] <- "se.outcome"
names(outcome)[8] <- "stat.outcome"
names(outcome)[9] <- "pval.outcome"
outcome$SNP <- sub(":\\S*", "", outcome$SNP)
outcome <- outcome[outcome$SNP != 1 , ]


outcome$effect_allele.outcome <- as.character(outcome$effect_allele.outcome)
outcome$other_allele.outcome <- as.character(outcome$other_allele.outcome)


outcome <- outcome %>% 
  mutate(effect_allele.outcome = case_when( nchar(effect_allele.outcome) < nchar(other_allele.outcome) ~ "-",
                                            TRUE ~ effect_allele.outcome),
         other_allele.outcome = case_when( nchar(effect_allele.outcome) > nchar(other_allele.outcome) ~ "-",
                                           TRUE ~ other_allele.outcome))


###### Merging datasets and matching effects ###### 

# For chromosome 12, merge only by "chr" and "pos", 
# outcome <- outcome[, -1]
pqtl_ocac <- merge(pqtl_exposure, outcome, by = c("chr", "pos", "SNP"))

pqtl_ocac <- pqtl_ocac %>%
  mutate(beta.exposure = case_when( effect_allele.exposure != effect_allele.outcome ~ beta.exposure*-1,
                                    effect_allele.exposure == effect_allele.outcome ~ beta.exposure*1))

write.table(pqtl_ocac, file = "~/Documents/eqtl.pqtl.ocac.project/pqtl/endometrioid/pqtl_chr22") # Change chromosome number

#### Merge all chromosomes ####

all_pqtl <- NULL
dat <- list.files("~/Documents/eqtl.pqtl.ocac.project/pqtl/clearcell/")
for(i in dat){
  x <- paste("~/Documents/eqtl.pqtl.ocac.project/pqtl/clearcell/", i, sep = "")
  y <- read.csv(x, sep = " ")
  all_pqtl <- rbind(all_pqtl, y)
}
rm(i, dat, x, y)

write.table(all_pqtl, file = "~/Documents/eqtl.pqtl.ocac.project/pqtl/results/all_clearcell")

###### Looping MR analysis through proteins ###### 

names(all_pqtl)[names(all_pqtl) == "protein"] <- "exposure"
names(all_pqtl)[names(all_pqtl) == "N"] <- "samplesize.exposure"
all_pqtl$samplesize.outcome <- "OCAC Study"
all_pqtl$id.exposure <- all_pqtl$exposure
all_pqtl$id.outcome <- "Ovarian Cancer"
all_pqtl$mr_keep <- TRUE
all_pqtl$outcome <- all_pqtl$id.outcome
pqtl_results <- mr_singlesnp(all_pqtl, single_method = "mr_wald_ratio", all_method = c("mr_ivw"))
pqtl_results <- pqtl_results[complete.cases(pqtl_results), ]

write.table(pqtl_results, file = "~/Documents/eqtl.pqtl.ocac.project/pqtl/results/mr_clearcell")

####

#### Extracting SNPs for BEDtools analysis ####

all <- read.csv("~/Documents/eqtl.pqtl.ocac.project/pqtl/pqtl_ocac_merge/all_pqtl_ocac", sep="")
proteins <- all[,"exposure"]
mr_loci <- all[all$exposure %in% proteins,]
mr_loci <- mr_loci [ mr_loci$SNP != "All - Inverse variance weighted", ]
y <- cis.exposure[ , c("chr", "SNP")]
all_loci <- all
all_loci <- merge(all_loci, y)
all_loci <- unique(all_loci) 
all_loci <- all_loci[, c("chr", "SNP", "exposure")]
all_loci$SNP <- as.numeric(levels(all_loci$SNP))[all_loci$SNP]
all_loci$start <- all_loci$SNP - 1
all_loci$stop <- all_loci$SNP
all_loci <- all_loci[, c("chr", "start", "stop", "exposure")]
all_loci$chr <- paste( "chr", all_loci$chr, sep = "")
write.table(all_loci, file = "~/Documents/eqtl.pqtl.ocac.project/eqtl/all_loci", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")





mr_loci <- merge( mr_loci, y)
mr_loci <- unique(mr_loci) 
mr_loci <- mr_loci[, c("chr", "SNP", "exposure")]
mr_loci$SNP <- as.numeric(levels(mr_loci$SNP))[mr_loci$SNP]
mr_loci$start <- mr_loci$SNP - 1
mr_loci$stop <- mr_loci$SNP
mr_loci <- mr_loci[, c("chr", "start", "stop", "exposure")]
mr_loci$chr <- paste( "chr", mr_loci$chr, sep = "")

write.table(mr_loci, file = "~/Documents/eqtl.pqtl.ocac.project/eqtl/mr_loci", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")