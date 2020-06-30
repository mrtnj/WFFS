
## Simple simulation with these features:
## * Discrete generations
## * 3000 females, 200 males
## * Selection only on the male side
## * One fully recessive lethal ~ 5% starting frequency
## * Assuming the lethal is unknown, but no matings on affected individuals


library(AlphaSimR)


source("R/simulation_functions.R")


## Read command line argument -- will be used as suffix for output files

args <- commandArgs(trailingOnly = TRUE)

rep_ix <- args[1]

    
print(rep_ix)

## Set up simulation

setup <- make_simulation(n_ind = 6000,
                         n_chr = 10,
                         h2 = 0.4)

founders <- setup$founderpop
simparam <- setup$simparam


## Pick a lethal variant

chosen_lethal <- pick_lethal(founders,
                             lethal_is = "snp",
                             simparam)
lethal_ix <- chosen_lethal$lethal_ix
other_snp_ix <- chosen_lethal$other_snp_ix


## Breeding simulation

n_gen <- 40

generations <- list(n_gen)

generations[[1]] <- founders

for (gen_ix in 2:n_gen) {
    
    generations[[gen_ix]] <-
        breed_unknown_lethal(generations[[gen_ix - 1]],
                             lethal_ix,
                             lethal_is = "snp",
                             simparam)
}


## Create simulation object

carrier_status <- lapply(generations,
                         carrier_test,
                         lethal_ix = lethal_ix,
                         lethal_is = "snp",
                         simparam)

carriers <- get_carriers(carrier_status)

stats <- get_stats(generations)

simulation_results <- list(simparam = simparam,
                           carrier_status = carrier_status,
                           carriers = carriers,
                           stats = stats,
                           lethal_ix = lethal_ix,
                           other_snp_ix = other_snp_ix)

saveRDS(generations,
        file = paste("simulations/simple_simulations_unknown/populations_",
                     rep_ix,
                     ".Rds",
                     sep = ""))

saveRDS(simulation_results, 
        file = paste("simulations/simple_simulations_unknown/results_",
                     rep_ix,
                     ".Rds",
                     sep = ""))



