library(ieugwasr)
library(here)
library(dplyr)
library(ggplot2)
library(data.table)

generate_vid <- function(d, ea = "ea", nea = "nea", eaf = "eaf", beta = "beta", rsid = "rsid", chr = "chr", position = "position") {
  toflip <- d[[ea]] > d[[nea]]
  d[[eaf]][toflip] <- 1 - d[[eaf]][toflip]
  d[[beta]][toflip] <- d[[beta]][toflip] * -1
  temp <- d[[nea]][toflip]
  d[[nea]][toflip] <- d[[ea]][toflip]
  d[[ea]][toflip] <- temp
  d[["rsido"]] <- d[[rsid]]
  d[[rsid]] <- paste0(d[[chr]], ":", d[[position]], "_", strtrim(d[[ea]], 5), "_", strtrim(d[[nea]], 5))
  d
}

read <- function(dati, main_hits) {
    a <- fread(here("data", dati$filename))
    a$A1 <- toupper(a$A1)
    a$A2 <- toupper(a$A2)
    a <- generate_vid(a, "A1", "A2", "MAF", "BETA", "SNP", "chr", "pos")
    head(a)
    b <- subset(a, SNP %in% main_hits$rsid) %>%
        mutate(trait = dati$trait, age = dati$age)
    return(b)
}


dat <- tribble(
    ~trait, ~age, ~filename,
    "ADHD", "Adult", "ADHD_M1_adults_forFUMA.txt.gz", 
    "ADHD", "Child", "ADHD_M1_child_forFUMA.txt.gz", 
    "ADHD", "All", "ADHD_M1_forFUMA.txt.gz", 
    "ASD", "Adult", "ASD_M1_adults_forFUMA.txt.gz", 
    "ASD", "Child", "ASD_M1_child_forFUMA.txt.gz", 
    "ASD", "All", "ASD_M1_forFUMA.txt.gz", 
    "Anxiety", "Adult", "Anxiety_M1_adults_forFUMA.txt.gz", 
    "Anxiety", "Child", "Anxiety_M1_child_forFUMA.txt.gz", 
    "Anxiety", "All", "Anxiety_M1_forFUMA.txt.gz", 
    "Depression", "Adult", "Depressionk_M1_adults_forFUMA.txt.gz", 
    "Depression", "Child", "Depressionk_M1_child_forFUMA.txt.gz", 
    "Depression", "All", "Depressionk_M1_forFUMA.txt.gz", 
    "Neuroticism", "All", "Neuro_M1_forFUMA.txt.gz", 
    "PLE", "All", "PLE_M1_forFUMA.txt.gz", 
    "Wellbeing", "All", "Wellbeing_M1_forFUMA.txt.gz"
)


vardat <- lapply(1:nrow(dat), function(i) {
    message(i)
    read(dat[i,], main_hits)
}) %>% bind_rows()


neur <- tophits("ebi-a-GCST90029028") %>% generate_vid() %>%
    mutate(trait = "Neuroticism")
ple <- tophits("ieu-b-5099") %>% generate_vid() %>%
    mutate(trait = "PLE")
depsnps <- scan(here("data", "depsnps.txt"), what="character")
dep <- associations(depsnps, "ebi-a-GCST005902") %>% generate_vid() %>%
    mutate(trait = "Depression")

well <- read.csv(here("data", "inst", "wellbeing.csv")) %>%
    select(rsid = SNP.ID, beta, se, chr=Chr., position=Position..bp., ea=EA, nea=NEA) %>%
    mutate(trait="Wellbeing", chr=as.character(chr)) %>%
    generate_vid()
well

asd <- read.csv(here("data", "inst", "asd.csv")) %>%
    filter(!is.na(se)) %>%
    select(rsid = Index.var, chr=CHR, position=BP, beta=beta, se, ea=A1, nea=A2) %>%
    mutate(trait="ASD", chr=as.character(chr)) %>%
    generate_vid()
asd

adhd <- read.csv(here("data", "inst", "adhd.csv")) %>%
    select(rsid = "rs.ID", chr=Chr., position=Position..bp., beta=OR, se=s.e., ea=A1, nea=A2) %>%
    mutate(trait="ADHD", beta=log(beta), chr=as.character(chr)) %>%
    generate_vid()
adhd

main_hits <- bind_rows(neur, ple, dep, well, asd, adhd)

save(main_hits, vardat, file=here("data", "polygenic.rda"))

# Trait Nsnp Reference
# Wellbeing 3 https://www.nature.com/articles/ng.3552
# ADHD 27 https://www.nature.com/articles/s41588-022-01285-8
# ASD 12 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6454898/
# Neuroticism 112 https://pubmed.ncbi.nlm.nih.gov/29892013/
# PLE 217 https://pubmed.ncbi.nlm.nih.gov/35396580/
# Depression 101 https://pubmed.ncbi.nlm.nih.gov/29662059/
# Educational attainment 74 https://pubmed.ncbi.nlm.nih.gov/27225129/
# BMI 458 https://www.medrxiv.org/content/10.1101/2021.06.28.21259622v1


## Analysis

analysis <- function(tr) {
    temp1 <- inner_join(main_hits %>% filter(trait == tr), vardat %>% filter(trait == tr), by=c("rsid"="SNP", "trait"))
    if(nrow(temp1) == 0) return(NULL)
    p <- ggplot(temp1, aes(x=beta, y=BETA)) +
    geom_errorbar(aes(ymin=BETA-1.96*SE, ymax=BETA+1.96*SE), colour="grey") +
    geom_errorbarh(aes(xmin=beta-1.96*se, xmax=beta+1.96*se), colour="grey") +
    geom_point() +
    geom_smooth(method="lm") +
    labs(x="Population-based main effect estimate", y="MZ difference effect estimate") +
    facet_grid(age ~ .)
    ggsave(p, file=here("images", paste0("polygenic_", tr, ".png")), height=length(unique(temp1$age))*3, width=5)

    o <- group_by(temp1, trait, age) %>%
        do({
            o <- summary(lm(BETA ~ beta, data=., weight=1/.$SE^2))$coef
            o <- as.data.frame(o)
            names(o) <- c("beta", "se", "tval", "pval")
            o$nsnp <- nrow(.)
            o[2,]
        })
    return(list(o=o, p=p))
}

res <- lapply(unique(dat$trait), analysis)

res2 <- lapply(res, \(x) x$o) %>% bind_rows()
res2
write.csv(res2, file=here("images", "polygenic.csv"))
