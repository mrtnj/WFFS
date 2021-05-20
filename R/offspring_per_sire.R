
## Look at SWB statistics to create a reasonable distribution
## for simulating the number of offspring per sire.

library(countreg)
library(dplyr)
library(MASS)
library(readr)
library(purrr)


foal <- map_dfr(system("ls foal_statistics/*.txt", intern = TRUE),
                read_tsv,
                col_types = "cn",
                .id = "file")

counts <- summarise(group_by(foal, name),
                    count = sum(offspring))

ex <- na.exclude(counts$count)



fit <- MASS::fitdistr(ex, "negative binomial")

fit


fakedata <- rnbinom(300, size = fit$estimate[1], mu = fit$estimate[2])

fakedata <- rztnbinom(300, mu = fit$estimate[2], size = fit$estimate[1])


get_offspring_distribution_qztnb <- function(n_sires,
                                             n_offspring) {

    quantiles <- qztnbinom(mu = 30,
                           size = 0.3,
                           p = seq(from = 0, to = 1 - 1/300, by = 1/300))
    
    o <- round(quantiles/sum(quantiles) * n_offspring)
 
    o   
}

o_total <- get_offspring_distribution_qztnb(300, 6000)
o_goal1 <- get_offspring_distribution_qztnb(300 * 0.6, 6000 * 0.6)
o_goal2 <- get_offspring_distribution_qztnb(300 * 0.4, 6000 * 0.4)




## Based on deciles

nonzero <- ex[ex > 0]

sorted <- sort(nonzero)

n <- length(sorted)


## Divide the sires into ten windows with 10% of sires each (deciles)

window_start <- round(seq(from = 1, to = n, by = n/10))
window_end <- c(window_start[-1] - 1, n)

offspring_in_window <- numeric(10)

for (window_ix in 1:10) {
    offspring_in_window[window_ix] <- sum(sorted[window_start[window_ix]:window_end[window_ix]])
}
    

## Find the fraction offspring in each decile of sires

fraction_in_window <- offspring_in_window/sum(nonzero)

get_offspring_distribution <- function(n_sires,
                                       n_offspring) {
    
    offspring_per_sire_window <- round(fraction_in_window * n_offspring / (n_sires / 10))
    
    rep(offspring_per_sire_window, each = n_sires / 10)

}

offspring_tot <- get_offspring_distribution(300,
                                            6000)

offspring_per_sire_window <- round(fraction_in_window * 6000 / (300 / 10))

## Breeding goal 1: 40% of sires

offspring_g1 <- get_offspring_distribution(300 * 0.4,
                                           6000 * 0.4)

## Breeding goal 2: 60% of sires

offspring_g2 <- get_offspring_distribution(300 * 0.6,
                                           6000 * 0.6)


