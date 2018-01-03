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

#### Looping MR analysis through proteins #### 

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

#### Extracting SNPs for BEDtools analysis ####

pqtl_loci <- pqtl_exposure[, c(2,3,10)]
pqtl_loci$start <- pqtl_loci$pos -1
pqtl_loci$stop <- pqtl_loci$pos
pqtl_loci <- pqtl_loci[, c("chr", "start", "stop", "protein")]
pqtl_loci$chr <- paste( "chr", pqtl_loci$chr, sep = "")
write.table(pqtl_loci, file = "~/Documents/eqtl.pqtl.ocac.project/pqtl/pqtl_loci", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#### Finding the most likely causal protein within each genome bin ####

#### Load files ####
pqtls <- pqtl_exposure
ld_all <- read.delim("~/Documents/eqtl.pqtl.ocac.project/eqtl/EUR/fourier_ls-chr2.bed")
oc_chr <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr2.txt")

#### Load LD blocks ####
n <- paste(1:length(ld_all$chr))
ld_all <- cbind(ld_all, n)
rm(n)

#### Organise OC dataframe ####
oc_chr <- oc_chr[ oc_chr$X1000G_SNPname != "",  ]
oc_chr <- oc_chr[, c(2, 3, 4, 9:12, 17:20, 29:32, 37:40, 49:52, 53:56)]
names(oc_chr)[names(oc_chr) == "X1000G_SNPname"] <- "rsid"
names(oc_chr)[names(oc_chr) == "Position"] <- "pos"
names(oc_chr)[names(oc_chr) == "Chromosome"] <- "chr"
oc_chr$rsid <- sub(":\\S*", "", oc_chr$rsid)


#### merge oc, ld, eqtl files ####
f = function(x) min(which(x >= ld_all$start & x <= ld_all$stop))
f = Vectorize(f)
oc_chr$x = f(oc_chr$pos)
ld_snp <- merge(oc_chr, pqtls, by = c("chr", "pos"))

#### Find top SNP per LD block ####

ld_snp <- ld_snp %>% 
  group_by(x, protein) %>%
  top_n(-1, pval.exposure)

write.table(ld_snp, file = "~/Documents/eqtl.pqtl.ocac.project/pqtl/results/ld_blocks/chr2", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

