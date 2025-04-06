#  Quantitative Analysis of the Pre-Mohacs Collection at the Hungarian National Archives ------------
#    Davor Salihovic, University of Antwerp                                                 
#                      davor.salihovic@uantwerpen.be

library(dplyr); library(tidyr)
library(ggplot2); library(patchwork); library(ggpubr)
library(rethinking)

# Clean the R environment
rm(list = ls())

## LOAD AND PREPARE DATA ----
load("data/data.RData")

d <- documents

d$year <- substr(documents$Date, 1, 4)
d$year <- as.numeric(documents$year)
d <- d[d$year < 2000,]

# 2. DOCUMENT ANALYSIS ----

## 2.2 ANALYSIS ----
### 2.2.1 Inflection points and plots ----
hst <- hist(d$year[d$year < 1526], breaks = 436)
minima <- hst$mids[which(diff(sign(diff(hst$counts))) == 2) + 1]
maxima <- hst$mids[which(diff(sign(diff(hst$counts))) == -2) + 1]

counts <- d %>% group_by(year) %>%
  summarise(n = n()) %>%
  filter(year <= 1526)
counts <- counts[!is.na(counts$year), ]

period0 <- counts[counts$year < 1280 & counts$year > 1200,]
period1 <- counts[counts$year < 1380 & counts$year > 1280,]
period2 <- counts[counts$year < 1430 & counts$year > 1380,]
period3 <- counts[counts$year < 1460 & counts$year > 1430,]
period4 <- counts[counts$year < 1490 & counts$year > 1460,]

minima <- rbind(
  period0[which.min(period0$n), ],
  period1[which.min(period1$n), ],
  period2[which.min(period2$n), ],
  period3[which.min(period3$n), ],
  period4[which.min(period4$n), ]
)

maxima <- rbind(
  period0[which.max(period0$n), ],
  period1[which.max(period1$n), ],
  period2[which.max(period2$n), ],
  period3[which.max(period3$n), ],
  period4[which.max(period4$n), ]
)

smoothed <- smooth.spline(hst$counts, nknots = 30) # This seems fine
ys <- smoothed$y
infl <- c(FALSE, diff(diff(ys) > 0) != 0)

inflpoints <- data.frame(x = hst$mids[infl], y = smoothed$y[infl])
inflpoints <- inflpoints[inflpoints$x > 1200,]

inf1 <- inflpoints[c(1, 3, 5, 7, 9),]
inf2 <- inflpoints[c(2, 4, 6, 8),]

smoothdat <- data.frame(x = hst$mids, y = smoothed$y)
smoothdat <- smoothdat[smoothdat$x > 1200,]

p1 <- ggplot(data = counts, aes(x = year, y = n)) +
  geom_line(col = "blue", alpha = .7) +
  geom_point(data = minima, aes(x = year, y = n), col = "grey20", fill = "grey20", alpha = .7) +
  geom_point(data = maxima, aes(x = year, y = n), col = "#8f0000", fill = "#8f0000", alpha = .7) +
  geom_text(data = minima, aes(x = year, y = n, label = year), hjust = -.2, family = "serif") +
  geom_text(data = maxima, aes(x = year, y = n, label = year), hjust = .2, vjust = -1, family = "serif") +
  theme_classic() +
  theme(text = element_text(family = "serif"),
        axis.title.x = element_text(family = "serif"),
        axis.title.y = element_text(family = "serif"),
        axis.text.x = element_text(family = "serif"),
        axis.text.y = element_text(family = "serif"),
        legend.text = element_text(family = "serif")) +
  labs(x = "Year", y = "Number of documents")

p2 <- ggplot(data = smoothdat, aes(x = x, y = y)) +
  geom_line(col = "#670091", alpha = .5) +
  geom_point(data = inf1, aes(x = x, y = y), alpha = .5) +
  geom_point(data = inf2, aes(x = x, y = y), alpha = .5) +
  geom_text(data = inf1, aes(x = x, y = y, label = x), family = "serif", hjust = .2, vjust = -1, size = 3) +
  geom_text(data = inf2, aes(x = x, y = y, label = x), family = "serif", hjust = .2, vjust = 2, size = 3) +
  theme_light() +
  theme(text = element_text(family = "serif"),
        axis.title.x = element_text(family = "serif"),
        axis.title.y = element_text(family = "serif"),
        axis.text.x = element_text(family = "serif"),
        axis.text.y = element_text(family = "serif"),
        legend.text = element_text(family = "serif")) +
  labs(x = "Year", y = "Number of documents")

p3 <- p1 + inset_element(p2, 0.02, .3, 0.6, .99)
p3

## 2.3 MODELS OF DOCUMENT GROWTH ----
set.seed(123)
model.data <- counts %>%
  mutate(yearstd = (year - mean(year)) / sd(year)) 

m1 <- ulam(
  alist(
    n ~ dpois(lambda),
    log(lambda) <- a + b*yearstd,
    a ~ dnorm(0, 10),
    b ~ dnorm(0, 10)
  ),
  data = model.data,
  cores = 4,
  chains = 4, 
  iter = 4000,
  warmup = 1000
)

traceplot(m1)
postcheck(m1)

saveRDS(m1, "document_growth.Rds")

pred <- seq(min(model.data$yearstd), max(model.data$yearstd), length.out = length(model.data$yearstd))
lambda <- link(m1, data = data.frame(yearstd = pred))
lambda.mu <- apply(lambda, 2, mean)
lambda.PI <- apply(lambda, 2, PI)

modelplotdata <- data.frame(
  n = model.data$n, ystd = model.data$yearstd, year = model.data$year,
  mu = lambda.mu, y = pred,
  u = lambda.PI[2, ], l = lambda.PI[1, ]
)

p4 <- ggplot(data = modelplotdata, aes(x = ystd, y = n)) +
  geom_point(col = "blue", alpha = .5) +
  geom_line(aes(x = pred, y = mu), lwd = .7) +
  geom_line(aes(x = pred, y = u), lty = 2) +
  geom_line(aes(x = pred, y = l), lty = 2) +
  theme_classic() +
  theme(text = element_text(family = "serif"),
        axis.title.x = element_text(family = "serif", size = 16),
        axis.title.y = element_text(family = "serif", size = 16),
        axis.text.x = element_text(family = "serif", size = 16),
        axis.text.y = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 16),
        plot.title = element_text(family = "serif", size = 16)) +
  annotate("text", x = -1, y = 500, label = expression(y[i] == e^{5.79 + x[i] * 1.42}), size = 6, 
           family = "serif") +
  labs(x = "Year (std)", y = "Number of documents")
p4

### 2.3.1 Rates of change ----
# Rate of change in intervals
exp(5.79+1.42*-1) - exp(5.79+1.42*-2) 
exp(5.79+1.42*0) - exp(5.79+1.42*-2)
exp(5.79+1.42*1) - exp(5.79+1.42*-1)  

# Average rate of change:
(exp(5.79+1.42*-2) - exp(5.79+1.42*-4)) / 2 
(exp(5.79+1.42*1) - exp(5.79+1.42*-2)) / (1 - (-2))

# Inst. rate of change, at points
exp(5.79+1.42*-2)*1.42
exp(5.79+1.42*0)*1.42
exp(5.79+1.42*1)*1.42

  # Tangent plot
a <- 5.79  # Intercept
b <- 1.42  # Coefficient

poisson_function <- function(x) {
  exp(a + b * x)
}

x_points <- c(-2, 0, 1)
slopes <- poisson_function(x_points) * b 
values <- poisson_function(x_points)          

tangent_lines <- function(x0, slope, y0, x) {
  y0 + slope * (x - x0)
}

tangent_data <- data.frame()
for (i in seq_along(x_points)) {
  x_range <- seq(x_points[i] - .3, x_points[i] + .3, length.out = 100)
  y_tangent <- tangent_lines(x_points[i], slopes[i], values[i], x_range)
  tangent_data <- rbind(tangent_data,
                        data.frame(x = x_range, y = y_tangent, group = as.factor(i)))
}

p5 <- ggplot(data = model.data, aes(x = yearstd, y = n)) +
  geom_point(color = rangi2, alpha = .25) +
  geom_line(data = tangent_data, aes(x = x, y = y, group = group),
            lwd = 1, color = "black") +
  theme_classic() +
  theme(text = element_text(family = "serif"),
        axis.title.x = element_text(family = "serif"),
        axis.title.y = element_text(family = "serif"),
        axis.text.x = element_text(family = "serif"),
        axis.text.y = element_text(family = "serif"),
        plot.title = element_text(family = "serif")) +
  scale_x_continuous(breaks = seq(min(model.data$yearstd), max(model.data$yearstd), .75),
                     labels = round(seq(min(model.data$yearstd), max(model.data$yearstd), .75) * sd(model.data$year) + mean(model.data$year))) +
  annotate("text", x = -2, y = 100, label = 27.13, size = 6, 
           family = "serif") +
  annotate("text", x = -0.2, y = 400, label = 464.36, size = 6, 
           family = "serif") +
  annotate("text", x = .6, y = 1400, label = 1921.1, size = 6, 
           family = "serif") +
  labs(x = "Year", y = "Number of documents")
p5

# Integrate over periods
dens <- density(d$year[!is.na(d$year)])
densf <- approxfun(dens$x, dens$y, rule = 2)
integrate(densf, lower = 800, upper = 1200)
integrate(densf, lower = 1200, upper = 1300)
integrate(densf, lower = 1400, upper = 1500)


## 2.4 METADATA ANALYSIS ----
dm <- d
names(dm) <- c("date", "dl_df", "issuer", "survival", "old_ref", "id", "abstract", "alt_date",
                "index", "language", "subject", "place", "year")

index <- dm[!is.na(dm$index),]
abstract <- dm[!is.na(dm$abstract),]
language <- dm[!is.na(dm$language),]
place <- dm[!is.na(dm$place),]
subject <- dm[!is.na(dm$subject),]
issuer <- dm[!is.na(dm$issuer),]

general_theme <- 
  theme(text = element_text(family = "serif"),
        axis.title.x = element_text(family = "serif"),
        axis.title.y = element_text(family = "serif"),
        axis.text.x = element_text(family = "serif"),
        axis.text.y = element_text(family = "serif"),
        plot.title = element_text(family = "serif"))

limits <- c(1100, 1526)

plots <- list()
metadata <- list(index, abstract, language, place, subject, issuer)

for (i in seq_along(metadata)) {
  data <- metadata[[i]]
  plots[[i]] <- ggplot(data = data, aes(year)) +
    geom_histogram(bins = length(unique(issuer$year)), fill = "blue", alpha = .6) +
    geom_vline(xintercept = 1433, color = "red", lwd = .7) +
    theme_classic() +
    general_theme +
    labs(x = "", y = "Frequency") +
    scale_x_continuous(limits = limits)
}

p6 <- ggarrange(plots[[1]], plots[[2]], plots[[3]],
                plots[[4]], plots[[5]], plots[[6]], nrow = 3, ncol = 2,
                align = "h", labels = LETTERS[1:6],
                font.label = list(size = 11))

# Original documents and later copies
dm$original <- ifelse(grepl("eredet", dm$survival, ignore.case = T), 1, 0)
propdm <- dm %>%
  group_by(year) %>%
  summarise(N_original = sum(original),
            total = n(),
            N_copies = total - N_original,
            Proportion_original = N_original / total)
propdm$Proportion_copies <- 1 - propdm$Proportion_original

propdm <- propdm %>%
  mutate(
    Log_original = log(N_original + 1),
    Log_copies = log(N_copies + 1),
    Sq_original = sqrt(N_original),
    Sq_copies = sqrt(N_copies)
  )

y_vals <- c("Log_original", "Log_copies", "Proportion_original", "Proportion_copies")
dplots <- list()

for (i in seq_along(y_vals)) {
  col_name <- y_vals[i]
  dplots[[i]] <- ggplot(data = propdm, aes(x = year, y = .data[[col_name]])) +
    geom_bar(stat = "identity", fill = "grey20", alpha = .4) +
    theme_classic() +
    general_theme +
    labs(x = "Year", y = col_name) +
    scale_x_continuous(breaks = seq(800, 1500, 100))
}

p7 <- ggarrange(dplots[[1]], dplots[[2]], dplots[[3]],
                dplots[[4]], nrow = 2, ncol = 2,
                align = "h", labels = LETTERS[1:4],
                font.label = list(size = 11))

ggsave("p1.tiff", p1, h = 8, w = 8)
ggsave("p2.tiff", p2, h = 8, w = 8)
ggsave("p3.tiff", p3, h = 8, w = 8)
ggsave("p4.tiff", p4, h = 8, w = 8)
ggsave("p5.tiff", p5, h = 8, w = 8)
ggsave("p6.tiff", p6, h = 8, w = 8)
ggsave("p7.tiff", p7, h = 8, w = 8)


