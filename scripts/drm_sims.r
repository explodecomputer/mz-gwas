# to setup data:
# cd ../data
# wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
# tar xzvf 1kg.v3.tgz
# plink \
# --bfile EUR \
# --chr 22 \
# --recode A \
# --out eur22 \
# --from-bp 0 --to-bp 17000000

library(glue)
library(data.table)
library(dplyr)
library(ggplot2)
library(here)
library(parallel)
options(mc.cores=30)

fast_assoc <- function(y, x)
{
    index <- is.finite(y) & is.finite(x)
    n <- sum(index)
    y <- y[index]
    x <- x[index]
    vx <- var(x)
    bhat <- stats::cov(y, x)/vx
    ahat <- mean(y) - bhat * mean(x)
    rsq <- (bhat * vx)^2/(vx * var(y))
    fval <- rsq * (n - 2)/(1 - rsq)
    tval <- sqrt(fval)
    se <- abs(bhat/tval)
    p <- stats::pf(fval, 1, n - 2, lower.tail = FALSE)
    return(list(ahat = ahat, bhat = bhat, se = se, fval = fval, pval = p, n = n))
}

test_drm <- function(g, y)
{
  y.i <- tapply(y, g, median, na.rm=T)
  z.ij <- abs(y - y.i[g+1])
  fast_assoc(z.ij, g) %>% as_tibble()
}

sim_full <- function(rsq, geno, n, cormat)
{
  i <- sample(1:ncol(geno), 1)
  y <- as.numeric(scale(geno[,i])) * sqrt(rsq) + rnorm(nrow(geno), 0, sqrt(1-rsq))
  res <- lapply(1:ncol(geno), function(i)
  {
      test_drm(geno[1:n,i], y[1:n])
  }) %>% bind_rows() %>% mutate(snp=1:n(), ldrsq=cormat[i,]^2)
  res %>% 
    mutate(fdr = p.adjust(pval, "fdr")) %>%
    summarise(
        vqtl = which.min(pval),
        minp = min(pval, na.rm=T),
        nfdr = sum(fdr < 0.05, na.rm=T),
        vqtl_ldrsq = ldrsq[vqtl]
        ) %>%
    mutate(rsq = rsq, n = n, qtl=i, af=sum(geno[,i])/(2*nrow(geno))) %>%
    ungroup()
}

geno <- fread(here("data", "eur22.raw"))
fam <- geno[,1:6]
geno <- as.matrix(geno[,-c(1:6)])
geno[1:10,1:10]
dim(geno)

cormat <- cor(geno)
dim(cormat)

## Determining spurious effects

param <- expand.grid(
    rsq=seq(0, 0.5, by=0.01),
    n=c(250,500),
    nrep=1:100
)

dim(param)
set.seed(1234)

res <- mclapply(1:nrow(param), function(i)
{
    message(i)
    sim_full(param$rsq[i], geno, param$n[i], cormat)
}) %>% bind_rows

save(res, file="drm_sims.rdata")

o <- sim_full(0.5, geno, 500, cormat)
o %>% ggplot(., aes(x=ldrsq, y = -log10(pval))) +
geom_point() +
geom_smooth()

set.seed(1234)
o <- sim_full(0.5, geno, 500, cormat)
o %>% ggplot(., aes(x=ldrsq, y = -log10(p))) +
geom_point() +
geom_smooth()


ggplot(res, aes(x=as.factor(rsq), y=nfdr)) +
geom_boxplot(aes(fill=as.factor(n)))

ggplot(res, aes(x=as.factor(rsq), y=-log10(minp))) +
geom_boxplot(aes(fill=as.factor(n)))

plot(I(-log10(minp)) ~ vqtl_ldrsq, res)
