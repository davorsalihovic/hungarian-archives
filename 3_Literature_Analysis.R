#  Quantitative Analysis of the Pre-Mohacs Collection at the Hungarian National Archives ------------
#    Davor Salihovic, University of Antwerp                                                 
#                      davor.salihovic@uantwerpen.be

library(dplyr); library(tidyr)
library(ggplot2); library(patchwork); library(ggpubr)
library(rethinking); library(car); library(performance)


# Clean the R environment
rm(list = ls())

## LOAD AND PREPARE DATA ----
load("data/data.RData")

szazadok <- read.csv("data/szazadok.csv") %>%
  mutate(journal = "sz")
tortsz <- read.csv("data/tort_sz.csv") %>%
  mutate(journal = "t")
literature <- rbind(tortsz, szazadok)

d <- documents

d$year <- substr(documents$Date, 1, 4)
d$year <- as.numeric(documents$year)
d <- d[d$year < 2000,]

general_theme <- 
  theme(text = element_text(family = "serif"),
        axis.title.x = element_text(family = "serif"),
        axis.title.y = element_text(family = "serif"),
        axis.text.x = element_text(family = "serif"),
        axis.text.y = element_text(family = "serif"),
        plot.title = element_text(family = "serif"))


# 3. LITERATURE ANALYSIS ----
arch.counts <- d %>% group_by(year) %>%
  summarise(n = n()) %>%
  mutate(
    decade = year %/% 10 * 10) %>%
  group_by(decade) %>%
  summarise(n = sum(n))

lit.counts <- literature %>% group_by(y_studied) %>%
  summarise(n = n()) %>%
  mutate(
    decade = y_studied %/% 10 * 10) %>%
  group_by(decade) %>%
  summarise(n = sum(n)) %>%
  filter(!(decade %in% c(1000, 1100, 1200, 1300, 1400, 1500)))

counts <- arch.counts %>%
  left_join(lit.counts, by = "decade", suffix = c("", "_lit"))

counts <- counts[complete.cases(counts), ]

p8 <- ggplot(data = counts, aes(x = n, y = n_lit)) +
  geom_jitter(aes(x = n, y = n_lit), color = "blue", alpha = .5, width = 0.5, height = 0.5) +
  geom_text(aes(label = decade), hjust = 0, vjust = -0.5, family = "serif", size = 3) +
  theme_classic() +
  general_theme +
  ylab("References to years in scholarship") +
  xlab("Archival prevalence")

p9 <- ggplot(data = counts, aes(x = n_lit, y = n)) +
  geom_jitter(aes(x = n_lit, y = n), color = "blue", alpha = .5, width = 0.5, height = 0.5) +
  geom_text(aes(label = decade), hjust = 0, vjust = -0.5, family = "serif", size = 3) +
  theme_classic() +
  general_theme +
  xlab("References to years in scholarship") +
  ylab("Archival prevalence")

p10 <- ggarrange(p8, p9, nrow = 2)
p10


## 3.1 MODEL FOR PREVALENCE IN SCHOLARSHIP ----
cor.test(counts$n, log(counts$n_lit))
cor.test(log(counts$n), log(counts$n_lit))

counts <- counts %>%
  mutate(log_n = log(n),
         log_n_lit = log(n_lit))

set.seed(123)
mPrev <- ulam(
  alist(
    log_n_lit ~ dnorm(mu, sigma),
    mu <- a + b*log_n,
    a ~ dnorm(1, 10),
    b ~ dnorm(2, 10),
    sigma ~ dunif(0, 100)
  ),
  data = counts,
  chains = 4, 
  iter = 4000, 
  warmup = 1000
)

set.seed(123)
mPrev1 <- ulam(
  alist(
    n_lit ~ dnorm(mu, sigma),
    mu <- a + n^b,
    a ~ dnorm(0, 10),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 100)
  ),
  data = counts,
  chains = 4, 
  iter = 4000, 
  warmup = 1000
)

pred <- data.frame(n = seq(min(counts$n), max(counts$n), length.out = length(counts$n)))
lin <- link(mPrev1, data = pred)
linmu <- apply(lin, 2, mean)
linpi <- apply(lin, 2, PI)

mplot <- cbind(counts$n, counts$n_lit, pred, linmu)
names(mplot) <- c("n", "n_lit", "pred_n", "mu")
mplot$u <- linpi[2, ]
mplot$l <- linpi[1, ]

p11 <- ggplot(data = mplot, aes(x = n, y = n_lit)) + 
  geom_point(col = "blue", alpha = .5) +
  geom_line(aes(x = pred_n, y = mu)) +
  geom_line(aes(x = pred_n, y = u), lty = 2) +
  geom_line(aes(x = pred_n, y = l), lty = 2) +
  theme_classic() +
  general_theme + 
  annotate("text", x = 1e4, y = 10, label = expression(y[i] == 1.46 + x[i]^0.34), size = 4.5, 
           family = "serif") +
  xlab("Archival prevalence") +
  ylab("References in scholarship") 
p11

## 3.2 DENSITY DIFFERENCE ----
literatureA <- literature[literature$y_studied < 1500 & 
                            !(literature$y_studied %in% c(1000, 1100, 1200, 1300, 1400, 1500)),]
densL <- density(literatureA$y_studied)
densLf <- approxfun(densL$x, densL$y, rule = 2)
integrate(densLf, lower = 800, upper = 1200)
integrate(densLf, lower = 1200, upper = 1300)
integrate(densLf, lower = 1400, upper = 1500)

densA <- density(d$year[!is.na(d$year)])

densdiff <- data.frame(year = densA$x, diff = (densL$y - densA$y))
p12 <- ggplot(data = densdiff, aes(x = year, y = diff)) +
  geom_bar(stat = "identity", fill = "blue", alpha = .5) +
  theme_classic() +
  general_theme +
  scale_x_continuous(breaks = seq(800, 1500, 50)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1398.2420, lty = 2) +
  ylab("Density discrepancy") +
  xlab("Time")
p12

p13 <- ggplot(data = literatureA, aes(x = y_studied)) +
  geom_histogram(aes(y = ..density..), colour = "blue", fill = "white", alpha = .5, col = NA, bins = 50) +
  geom_density(fill = "red", alpha = .4, col = NA, bw = 29.17) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(700, 1526), breaks = seq(700, 1526, 150)) +
  theme_classic() +
  general_theme +
  ylab("Density") +
  xlab("Years")

p14 <- ggplot(data = d, aes(x = year)) +
  geom_histogram(aes(y = ..density..), colour = "blue", fill = "white", alpha = .5, col = NA, bins = 50) +
  geom_density(fill = "red", alpha = .4, col = NA, bw = 4.8) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(700, 1526), breaks = seq(700, 1526, 150)) +
  theme_classic() +
  general_theme +
  ylab("Density") +
  xlab("Years")

p15 <- ggarrange(p13, p14, nrow = 2, labels = LETTERS[1:2])
p16 <- p12 + inset_element(p15, .01, .01, .6, .6)
p16


## 3.3 MODELS FOR SHIFT IN TEMPORAL FOCUS ----
litA <- literatureA %>%
  mutate(ystd = (year - mean(year)) / sd(year),
         journal_n = ifelse(journal == "t", 1, 0))

lm1 <- lm(y_studied ~ ystd, data = litA)
lm2 <- lm(y_studied ~ ystd + journal_n, data = litA)
lm3 <- lm(y_studied ~ ystd*journal_n, data = litA)

mTemp0 <- ulam(
  alist(
    y_studied ~ dnorm(mu, sigma),
    mu <- a + bY*ystd,
    a ~ dnorm(0, 100),
    bY ~ dnorm(0, 100),
    sigma ~ dunif(0, 100)
  ),
  data = litA,
  chains = 4, 
  iter = 4000, 
  warmup = 1000)

mTemp1 <- ulam(
  alist(
    y_studied ~ dnorm(mu, sigma),
    mu <- a + bY*ystd + bJ*journal_n + bYJ*ystd*journal_n,
    a ~ dnorm(0, 100),
    c(bY, bJ, bYJ) ~ dnorm(0, 100),
    sigma ~ dunif(0, 100)
  ),
  data = litA,
  chains = 4, 
  iter = 4000, 
  warmup = 1000
)

traceplot(mTemp0)
traceplot(mTemp1)

check_autocorrelation(lm3)

litPred1 <- predict.lm(lm3, newdata = data.frame(ystd = litA$ystd[litA$journal_n == 1], journal_n = 1), interval = "confidence")
litPred0 <- predict.lm(lm3, newdata = data.frame(ystd = litA$ystd, journal_n = 0), interval = "confidence")

litPred1 <- cbind(litPred1, x = litA$ystd[litA$journal_n == 1])
litPred0 <- cbind(litPred0, x = litA$ystd, journal_n = 0)

p17 <- ggplot(data = litA, aes(x = ystd, y = y_studied)) +
  geom_point(aes(col = journal_n, alpha = .5)) +
  geom_line(data = litPred1, aes(x = x, y = fit), lty = 1, col = "purple") +
  geom_line(data = litPred0, aes(x = x, y = fit), lty = 1, col = "grey40") +
  geom_ribbon(data = litPred0, aes(ymin = lwr, ymax = upr, x = x), inherit.aes = F, 
              fill = "grey40", alpha = .2) +
  geom_ribbon(data = litPred1, aes(ymin = lwr, ymax = upr, x = x), inherit.aes = F, 
              fill = "purple", alpha = .2) +
  theme_classic() +
  general_theme +
  theme(legend.position = "none") +
  labs(x = "Issue year (std)", y = "") +
  scale_x_continuous(breaks = seq(min(litA$ystd), max(litA$ystd), .5),
                     labels = round(seq(min(litA$ystd), max(litA$ystd), .5) * sd(litA$year) + mean(litA$year)))

ggsave("p10.tiff", p10, h = 7, w = 7)
ggsave("p11.tiff", p11, h = 7, w = 7)
ggsave("p16.tiff", p16, h = 7, w = 7)
ggsave("p17.tiff", p17, h = 6, w = 6)
