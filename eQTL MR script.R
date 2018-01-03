### eQTL data management

## Load dependancies
library(TwoSampleMR)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

## load in exposure data sets

cis.exposure <- read.delim("~/Documents/eqtl.pqtl.ocac.project/eqtl/cis_eQTL_table_conditional_ALL")

## Drop non-informative data
cis.exposure <- cis.exposure[ cis.exposure$Gene != ".", ]

## Rename columns 
names(cis.exposure)[names(cis.exposure) == "rs.SNP"] <- "rsid"
names(cis.exposure)[names(cis.exposure) == "SNPs.used.for.conditioning"] <- "Conditioning"
names(cis.exposure)[names(cis.exposure) == "Beta"] <- "beta.exposure"
names(cis.exposure)[names(cis.exposure) == "Gene"] <- "gene.exposure"
names(cis.exposure)[names(cis.exposure) == "P"] <- "pval.exposure"


## Sort datasets by gene, then conditioning, then P-value. 
cis.exposure <- cis.exposure[order(cis.exposure[,5], cis.exposure[,1], cis.exposure[,7]),]

## Take the top snp from each gene/conditioning group

cis.exposure <- cis.exposure %>% 
  group_by(gene.exposure, Conditioning) %>% 
  top_n(-1, pval.exposure) %>%
  filter(row_number() == n())


## Calculate test statistics and standard error from p-values and attach to dataframe 
cis.exposure$beta.exposure / ( qnorm( cis.exposure$pval.exposure /2)) -> se.exposure
se.exposure <- sqrt(se.exposure^2)
cis.exposure <- cbind.data.frame( cis.exposure, se.exposure)
rm( se.exposure)


## Re-order the columns 
cis.exposure <- subset(cis.exposure, select = c( "Conditioning", "rsid", "Probeset","gene.exposure" , "SNP", "beta.exposure",
                                                 "se.exposure", "pval.exposure", "FDR"))

## Change FDR to a numeric
cis.exposure$FDR <- sub("<", "", cis.exposure$FDR)
cis.exposure$FDR <- as.numeric(as.character(cis.exposure$FDR))


## Extract the effect alleles/chr/pos

effect_allele.exposure <- str_sub(cis.exposure$SNP, -1)
cis.exposure <- cbind.data.frame(cis.exposure, effect_allele.exposure)
cis.exposure <- subset(cis.exposure, select = c( "Conditioning", "rsid","gene.exposure", "SNP", "effect_allele.exposure",
                                                 "beta.exposure", "se.exposure", "pval.exposure", "FDR" ))

rm(effect_allele.exposure)

cis.exposure$chr <- sub(":.*", "", cis.exposure$SNP)
cis.exposure$pos <- sub("(.*?:)", "", cis.exposure$SNP)
cis.exposure$pos <- sub("[:punct:].*", "", cis.exposure$pos )
cis.exposure$SNP <- sub("_.*", "", cis.exposure$pos )
cis.exposure$pos <- NULL 
cis.exposure <- subset(cis.exposure, select = c( "Conditioning", "rsid","gene.exposure", "chr", "SNP", "effect_allele.exposure",
                                                 "beta.exposure", "se.exposure", "pval.exposure", "FDR" ))


### Separate statistically significant in cis

stringent <- subset(cis.exposure, pval.exposure <= 1e-5)
## indels <- subset(stringent, effect_allele.exposure != "R" & effect_allele.exposure == "D" | effect_allele.exposure == "I")
## r <- subset(stringent, effect_allele.exposure == "R")
stringent.no.r <- subset(stringent, effect_allele.exposure != "R")


### Load in OCAC summary data

ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr1.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr2.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr3.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr4.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr5.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr6.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr7.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr8.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr9.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr10.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr11.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr12.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr13.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr14.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr15.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr16.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr17.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr18.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr19.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr20.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr21.txt")
ocac.outcome <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr22.txt")


### Extracting overall summary data from OCAC

outcome <- ocac.outcome[, c(2,3,4,5,6,9,10,12)]
outcome <- outcome[ ocac.outcome$X1000G_SNPname != "", ]
names(outcome)[names(outcome) == "X1000G_SNPname"] <- "rsid"
names(outcome)[names(outcome) == "Position"] <- "SNP"
names(outcome)[names(outcome) == "Chromosome"] <- "chr"
names(outcome)[names(outcome) == "Effect"] <- "effect_allele.outcome"
names(outcome)[names(outcome) == "Baseline"] <- "other_allele.outcome"
names(outcome)[6] <- "beta.outcome"
names(outcome)[7] <- "se.outcome"
names(outcome)[8] <- "pval.outcome"


## Extract rsID from outcome data

outcome$rsid <- sub(":\\S*", "", outcome$rsid)

### Preparing outcome data for analysis 

outcome$effect_allele.outcome <- as.character(outcome$effect_allele.outcome)
outcome$other_allele.outcome <- as.character(outcome$other_allele.outcome)

outcome <- outcome %>% 
  mutate(effect_allele.outcome = case_when( nchar(effect_allele.outcome) < nchar(other_allele.outcome) ~ "D",
                                            nchar(effect_allele.outcome) > nchar(other_allele.outcome) ~ "I",
                                            nchar(effect_allele.outcome) > 1                           ~ substr(effect_allele.outcome, 1, 1),
                                            TRUE ~ effect_allele.outcome))


### Merge exposure and outcome datasets and prepare for analysis
stringent.no.r <- subset(stringent, effect_allele.exposure != "R")
stringent.no.r <- merge(stringent.no.r, outcome, by = c("chr", "SNP"))
stringent.no.r <- stringent.no.r[stringent.no.r$effect_allele.exposure == stringent.no.r$effect_allele.outcome | 
                                   stringent.no.r$effect_allele.exposure == stringent.no.r$other_allele.outcome, ]
stringent.no.r$id.exposure <- "eQTL"
stringent.no.r$id.outcome <- "Ovarian Cancer"
stringent.no.r$mr_keep <- TRUE
stringent.no.r$outcome <- "Ovarian Cancer"
stringent.no.r$samplesize.exposure <- 4896
stringent.no.r$samplesize.outcome <- "OCAC Study"

## Dealing with problematic SNPs

stringent.no.r <- stringent.no.r[stringent.no.r$effect_allele.outcome != "I" | stringent.no.r$gene.exposure != "BNC2", ]
stringent.no.r <- stringent.no.r[stringent.no.r$Conditioning != 1 | stringent.no.r$gene.exposure != "B2M", ]


## Automate/loop gene anlysis by chromosome

stringent.no.r[ "3" == stringent.no.r[, "chr"], "gene.exposure"] -> chr3 ## CHANGE CHR NUMBERS
unique(chr3) -> chr3 ## CHANGE CHR NUMBERS
store.results <- cbind("id.exposure", "id.outcome",  "outcome", "exposure" ,"method",
                               "nsnp", "b", "se", "pval", deparse.level = 2 )
results.single.snp <- store.results[0, ]
results.lenient.sig <- store.results[0, ]

for (gene in chr3){ ## CHANGE CHR
  print(gene)
  stringent.no.r$exposure <- gene
  subset(stringent.no.r, gene.exposure == gene) -> gene
  gene$effect_allele.outcome <- factor(gene$effect_allele.outcome, levels = levels(gene$effect_allele.exposure))
  gene[gene$effect_allele.exposure != gene$effect_allele.outcome, "beta.exposure"]* -1 -> 
    gene[gene$effect_allele.exposure != gene$effect_allele.outcome, "beta.exposure"]
  
  
##  mr(gene, method_list = c("mr_wald_ratio", "mr_ivw" )) -> results
  mr_singlesnp(gene, single_method = "mr_wald_ratio", all_method = c("mr_ivw")) -> results
  results <- results[complete.cases(results), ]
  results.single.snp <- rbind(results.single.snp, results)

}

#########################################################################################
## Manual gene analysis

CARD16 <- stringent.no.r[stringent.no.r$gene.exposure == "CARD16", ]
CARD16$exposure <- "CARD16"
CARD16$effect_allele.outcome <- factor(CARD16$effect_allele.outcome, levels = levels(CARD16$effect_allele.exposure))
CARD16[CARD16$effect_allele.exposure != CARD16$effect_allele.outcome, "beta.exposure"]* -1 -> 
  CARD16[CARD16$effect_allele.exposure != CARD16$effect_allele.outcome, "beta.exposure"]
x <- mr(CARD16) 
mr_scatter_plot(x, CARD16)


#########################################################################################
## Exporting files
length(unique(stringent.no.r$gene.exposure))
length(unique(results.single.snp$exposure))
write.table(results.single.snp, file = "~/Documents/eqtl.pqtl.ocac.project/eqtl/endometrioid_results/chr22")



## Pull in all overall results

dat <- list.files("~/Documents/eqtl.pqtl.ocac.project/eqtl/endometrioid_results/")

for (i in dat){
  y <- paste("~/Documents/eqtl.pqtl.ocac.project/eqtl/endometrioid_results/", i, sep = "")
  assign(i, read.csv(y, sep = " "))
}

single_snp <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22) 

write.table(single_snp, file = "~/Documents/eqtl.pqtl.ocac.project/eqtl/endometrioid_results/endometrioid_results")

##########################################################################################
## Importing results files all together

overall_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/overall_results/overall_results", sep="")
serous_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/serous_results/serous_results", sep="")
serous_hg_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/serous_hg_results/serous_hg_results", sep="")
ser_lg_lmp_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/ser_lg_lmp_results/ser_lg_lmp_results", sep="")
serous_lmp_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/serous_lmp_results/serous_lmp_results", sep="")
serouslowgrade_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/serouslowgrade_results/serouslowgrade_results", sep="")
mucinous_all_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/mucinous_all_results/mucinous_all_results", sep="")
mucinous_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/mucinous_results/mucinous_results", sep="")
mucinous_lmp_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/mucinous_lmp_results/mucinous_lmp_results", sep="")
lmp_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/lmp_results/lmp_results", sep="")
clearcell_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/clearcell_results/clearcell_results", sep="")
endometrioid_results <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/endometrioid_results/endometrioid_results", sep="")
non_overlapping <- read.delim("~/Documents/eqtl.pqtl.ocac.project/eqtl/non_overlapping", header=FALSE)
overlapping <- read.delim("~/Documents/eqtl.pqtl.ocac.project/eqtl/overlapping", header=FALSE)


##########################################################################################
## Extracting SNPs for BEDtools analysis 

all <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/overall_results/overall_results", sep="")
sig <- read.csv("~/Documents/eqtl.pqtl.ocac.project/eqtl/overall_results_lenient/sig.lenient", sep="")
genes <- sig[sig$pval < 0.05, "exposure"]
mr_loci <- all[all$exposure %in% genes,]
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

##########################################################################################
## Finding the most likely causal gene within each genome bin 

##Load files
eqtls <- read.delim("~/Documents/eqtl.pqtl.ocac.project/eqtl/cis_eQTL_table_conditional_ALL")
ld_all <- read.delim("~/Documents/eqtl.pqtl.ocac.project/eqtl/EUR/fourier_ls-chr22.bed")
oc_chr <- read.csv("~/Documents/eqtl.pqtl.ocac.project/ocac.summary/Summary_chr22.txt")


##sort eqtl file. 

eqtls <- eqtls[ eqtls$Gene != ".", ]
names(eqtls)[names(eqtls) == "rs.SNP"] <- "rsid"
names(eqtls)[names(eqtls) == "SNPs.used.for.conditioning"] <- "Conditioning"
names(eqtls)[names(eqtls) == "Beta"] <- "beta.exposure"
names(eqtls)[names(eqtls) == "Gene"] <- "gene.exposure"
names(eqtls)[names(eqtls) == "P"] <- "pval.exposure"
effect_allele.exposure <- str_sub(eqtls$SNP, -1)
eqtls <- cbind.data.frame(eqtls, effect_allele.exposure)
rm(effect_allele.exposure)
eqtls$chr <- sub(":.*", "", eqtls$SNP)
eqtls$pos <- sub("(.*?:)", "", eqtls$SNP)
eqtls$pos <- sub("[:punct:].*", "", eqtls$pos )
eqtls$SNP <- sub("_.*", "", eqtls$pos )
eqtls$pos <- NULL 
eqtls$Probeset <- NULL
eqtls$FDR <- NULL
eqtls$effect_allele.exposure <- NULL
eqtls$rsid <- NULL 

## Load LD blocks
n <- paste(1:length(ld_all$chr))
ld_all <- cbind(ld_all, n)
rm(n)

## Organise OC dataframe
oc_chr <- oc_chr[oc_chr$endometrioid_pvalue != -99 & oc_chr$X1000G_SNPname != "",  ]
oc_chr <- oc_chr[, c(2, 3, 4, 53, 54, 56)]
names(oc_chr)[names(oc_chr) == "X1000G_SNPname"] <- "rsid"
names(oc_chr)[names(oc_chr) == "Position"] <- "SNP"
names(oc_chr)[names(oc_chr) == "Chromosome"] <- "chr"
oc_chr$rsid <- sub(":\\S*", "", oc_chr$rsid)


## merge oc, ld, eqtl files
f = function(x) min(which(x >= ld_all$start & x <= ld_all$stop))
f = Vectorize(f)
oc_chr$x = f(oc_chr$SNP)
ld_snp <- merge(oc_chr, eqtls, by = c("chr", "SNP"))

## Find top SNP per LD block

ld_snp <- ld_snp %>% 
  group_by(x) %>% 
  top_n(-1, endometrioid_pvalue) %>%
  group_by(gene.exposure) %>%
  top_n(-1, pval.exposure)

write.table(ld_snp, file = "~/Documents/eqtl.pqtl.ocac.project/eqtl/ld_snps/endometrioid/ld_chr22", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
