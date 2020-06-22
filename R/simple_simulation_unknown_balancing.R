
## Simple simulation with these features:
## * Discrete generations
## * 3000 females, 200 males
## * Selection only on the male side
## * One fully recessive lethal ~ 5% starting frequency
## * Assuming the lethal is unknown, but no matings on affected individuals
## * Pleiotropy between lethal and breeding goal


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
                             lethal_is = "qtl",
                             simparam)
lethal_ix <- chosen_lethal$lethal_ix
other_snp_ix <- chosen_lethal$other_snp_ix
founder_lethal_frequency <- chosen_lethal$founder_lethal_frequency


## Breeding simulation

n_gen <- 40

generations <- list(n_gen)

generations[[1]] <- founders

for (gen_ix in 2:n_gen) {
    
    generations[[gen_ix]] <-
        breed_unknwon_lethal(generations[[gen_ix - 1]],
                             lethal_ix,
                             lethal_is = "qtl",
                             simparam)
}


## Create simulation object

carrier_status <- lapply(generations,
                         carrier_test,
                         lethal_ix = lethal_ix,
                         lethal_is = "qtl",
                         simparam)

carriers <- get_carriers(carrier_status)

stats <- get_stats(generations)

simulation_results <- list(simparam = simparam,
                           carrier_status = carrier_status,
                           carriers = carriers,
                           stats = stats,
                           lethal_ix = lethal_ix,
                           other_snp_ix = other_snp_ix,
                           lethal_qtl_effect = simparam$traits[[1]]@addEff[lethal_ix],
                           other_qtl_effects = simparam$traits[[1]]@addEff[-lethal_ix],
                           founder_lethal_frequency = founder_lethal_frequency,
                           founder_h2 = varG(founders)/varP(founders))

simulation_results$founder_variance_explained_by_lethal <-
    2 * simulation_results$lethal_qtl_effect^2 *
    simulation_results$founder_lethal_frequency *
    (1 - simulation_results$founder_lethal_frequency) *
    simulation_results$founder_h2

saveRDS(generations,
        file = paste("simulations/simple_simulations_balancing/populations_",
                     rep_ix,
                     ".Rds",
                     sep = ""))

saveRDS(simulation_results, 
        file = paste("simulations/simple_simulations_balancing/results_",
                     rep_ix,
                     ".Rds",
                     sep = ""))


    
    