
## Simple simulation to test h2 and intensity

library(AlphaSimR)
library(dplyr)
library(ggplot2)


source("R/simulation_functions.R")


results <- vector(mode = "list", length = 3 * 3)
k <- 1


for (h2 in c(0.01, 0.1, 0.3)) {
    
    for (sires in c(100, 300, 600)) {
        
        cat(h2, " ", sires, "\n")
        
        setup <- make_simulation(n_ind = 6000,
                                 n_chr = 10,
                                 n_qtl = 10,
                                 n_snp = 400,
                                 h2 = h2)
        
        founders <- setup$founderpop
        simparam <- setup$simparam
        
        chosen_lethal <- pick_lethal(founders,
                                     lethal_is = "snp",
                                     simparam)
        lethal_ix <- chosen_lethal$lethal_ix
        
        
        n_gen <- 40
        
        generations <- list(n_gen)
        
        generations[[1]] <- founders
        
        for (gen_ix in 2:n_gen) {
            
            generations[[gen_ix]] <-
                breed_avoiding_carrier_x_carrier(generations[[gen_ix - 1]],
                                                 lethal_ix,
                                                 lethal_is = "snp",
                                                 n_sires = sires,
                                                 simparam)
        }
        
        stats <- get_stats(generations)
        stats$h2 <- h2
        stats$sires <- sires
        
        results[[k]] <- stats
        k <- k + 1
    }
}


saveRDS(results,
         file = "simulations/simple_simulations/h2_intensity_results.Rds")
        


combined <- Reduce(rbind, results)

plot_gain <- qplot(y = mean_g, x = generation, colour = paste(h2, sires), data = combined,
                   geom = "line")


## Average gain per generation the first ten generations

sd_per_year <- lapply(results, get_sd_per_year)

lapply(sd_per_year, function(x) mean(x[1:10]))
