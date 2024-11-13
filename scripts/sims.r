library(dplyr)
library(ggplot2)
library(tidyr)
library(simulateGP)
library(here)
library(parallel)
library(furrr)
mc.cores <- 60
plan("multicore", workers = 60)

sim_pop <- function(n, beta1, beta2, af, h2)
{
  g <- rbinom(n, 2, af)
  prs <- g * beta1
  vg <- rnorm(n, 0, h2)
  v <- rnorm(n, 0, beta2 * g)
  ve <- rnorm(n, 0, sqrt(1 - var(vg) - var(v) - var(prs)))
  y <- prs + v + vg + ve
  return(tibble(
    g, y
  ))
}

sim_mz <- function(n, beta1, beta2, af, h2)
{
  g <- rbinom(n, 2, af)
  prs <- g * beta1
  vg <- rnorm(n, 0, h2)
  v1 <- rnorm(n, 0, beta2 * g)
  ve1 <- rnorm(n, 0, sqrt(1 - var(vg) - var(v1) - var(prs)))
  y1 <- prs + v1 + vg + ve1
  v2 <- rnorm(n, 0, beta2 * g)
  ve2 <- rnorm(n, 0, sqrt(1 - var(vg) - var(v2) - var(prs)))
  y2 <- prs + v2 + vg + ve2
  return(tibble(
    g, y1, y2
  ))
}

test_drm <- function(g, y)
{
  y.i <- tapply(y, g, median, na.rm=T)  
  z.ij <- abs(y - y.i[g+1])
  fast_assoc(z.ij, g) %>%
    as_tibble() %>%
    mutate(method="drm")
}

test_mz <- function(g, y1, y2)
{
  yd1 <- abs(y1-y2)
  r1 <- fast_assoc(yd1, g) %>%
    as_tibble() %>%
    mutate(method="mzdiff")
  r1
}

param <- expand.grid(
    beta1 = 0,
    beta2 = 0.115,
    h2 = seq(0.05, 0.95, by=0.05),
    af = 0.3,
    rep=1:10000
)
dim(param)

res <- future_map(1:nrow(param), function(i)
  {
  a1 <- do.call(sim_mz, param[i,] %>% mutate(n=10000) %>% select(-c(rep)))
  a2 <- do.call(sim_mz, param[i,] %>% mutate(n=500000) %>% select(-c(rep)))
  a3 <- do.call(sim_mz, param[i,] %>% mutate(n=20000) %>% select(-c(rep)))
  if(any(is.na(a1$y1)) | any(is.na(a1$y2)) | any(is.na(a2$y1)) | any(is.na(a2$y2)))
  {
    return(NULL)
  }
  bind_rows(
    test_mz(a1$g, a1$y1, a1$y2),
    test_drm(a2$g, a2$y1),
    test_drm(a3$g, a3$y1)
  ) %>%
    bind_cols(., param[i,])
}, .progress=TRUE) %>%
  bind_rows()

save(res, file=here("scripts/mzpowersim.rdata"))

##

load(here("scripts/mzpowersim.rdata"))

p2 <- bind_rows(res) %>%
group_by(h2, method, n) %>%
summarise(pow = sum(pval < 5e-8)/n(), fval=mean(fval)) %>%
mutate(method = case_when(method == "drm" ~ paste0("Population, n = ", n/1000, "k"), TRUE ~ paste0("MZ pairs, n = ", n/1000, "k"))) %>%
ggplot(., aes(x=h2, y=fval)) +
  geom_line(aes(colour=method)) +
  geom_point(aes(colour=method)) +
  labs(y="vQTL F-statistic (log10 scale)", x="Narrow sense heritability", colour="Design") +
  scale_colour_brewer(type="qual") +
  scale_y_log10()
p2
ggsave(p2, file=here("images/pow_diffn.pdf"), height=6, width=11)

p2 <- bind_rows(res2, res2a) %>%
filter(n != 500000) %>%
mutate(n2 = case_when(n == 10000 ~ 20000, TRUE ~ 20000)) %>%
mutate(method = case_when(method == "drm" ~ paste0("Population, n = ", n/1000, "k"), TRUE ~ paste0("MZ pairs, n = ", n/1000, "k"))) %>%
filter(h2 < 0.9) %>%
ggplot(., aes(x=h2, y=fval)) +
  geom_smooth(aes(colour=method)) +
  labs(y="vQTL discovery power", x="Narrow sense heritability", colour="Design") +
  scale_colour_brewer(type="qual")
p2
ggsave(p2, file=here("images/pow_diffn.pdf"), height=6, width=11)


sim_mz2 <- function(g, beta1, beta2, h2)
{
  n <- length(g)
  prs <- g * beta1
  vg <- rnorm(n, 0, h2)
  v1 <- rnorm(n, 0, beta2 * g)
  ve1 <- rnorm(n, 0, sqrt(1 - var(vg) - var(v1) - var(prs)))
  y1 <- prs + v1 + vg + ve1
  v2 <- rnorm(n, 0, beta2 * g)
  ve2 <- rnorm(n, 0, sqrt(1 - var(vg) - var(v2) - var(prs)))
  y2 <- prs + v2 + vg + ve2
  return(tibble(
    g, y1, y2
  ))
}

test_mz <- function(g, y1, y2)
{
  yd1 <- abs(y1-y2)
  r1 <- summary(lm(yd1 ~ g))$coef %>%
    as_tibble() %>%
    slice(2) %>%
    mutate(method="mzdiff")
  r1
}

gendatp <- function(n, p1, p2, p3, r1)
{
	# dat <- simulateGP:::simulate_geno(n, r1, p1, p2) %>% as_tibble
	dat <- correlated_binomial(n, p1, p2, r1) %>% as_tibble()
	names(dat) <- c("g1", "g2")
	dat$g3 <- rbinom(n, 1, p3)
	return(dat)
}

run_simp_mz <- function(param, i)
{
	set.seed(i*10)
	dat <- gendatp(param$n[i], param$p1[i], param$p2[i], param$p3[i], param$r1[i])
	mzdat <- sim_mz2(dat$g1, param$beta1[i], param$beta2[i], param$h2)
	#x <- dat$y1 + rnorm(nrow(dat), sd=sd(dat$y1)/2)
	o1 <- test_drm(dat$g1, mzdat$y1)
	o2 <- test_drm(dat$g2, mzdat$y1)
	o3 <- test_drm(dat$g3, mzdat$y1)
	m1 <- test_mz(dat$g1, mzdat$y1, mzdat$y2)
	m2 <- test_mz(dat$g2, mzdat$y1, mzdat$y2)
	m3 <- test_mz(dat$g3, mzdat$y1, mzdat$y2)
	param$drm1[i] <- o1$`Pr(>|t|)`
	param$drm2[i] <- o2$`Pr(>|t|)`
	param$drm3[i] <- o3$`Pr(>|t|)`
	param$mz1[i] <- m1$`Pr(>|t|)`
	param$mz2[i] <- m2$`Pr(>|t|)`
	param$mz3[i] <- m3$`Pr(>|t|)`
	return(param[i,])
}

param <- expand.grid(
	p1=0.1,
	p2=0.1,
	p3=0.5,
	p4=0.1,
	n=1000,
	r1=seq(0, 1, by=0.2),
	beta1=1,
	beta2=0,
	h2=0.5,
	sim=1:500,
	r2=NA,
	F=NA,
	drm1=NA,
	drm2=NA,
	drm3=NA,
	mz1=NA,
	mz2=NA,
	mz3=NA
)

resmz <- lapply(1:nrow(param), function(i) run_simp_mz(param, i)) %>% bind_rows()
resmz %>% str

p3 <- resmz %>% 
  dplyr::select(r1, MZ=mz2, pop=drm2) %>% gather(., "key", "value", MZ, pop) %>%
  mutate(key = case_when(key == "pop" ~ "Population", TRUE ~ "MZ")) %>%
  ggplot(., aes(x=r1, y=-log10(value))) +
  geom_boxplot(aes(fill=as.factor(r1))) +
  scale_fill_brewer(type="seq") +
  facet_grid(. ~ key) +
  labs(y="vQTL -log10 p under the null", x="LD between tagging\nvariant and causal variant", fill="LD")
p3
ggsave(p3, file=here("images/inflation.pdf"), height=6, width=11)

save(resmz, res1, res2, res2a, file=here("images/sim.rdata"))
