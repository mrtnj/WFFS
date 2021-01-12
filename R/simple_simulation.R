
## Simple simulation with these features:
## * Discrete generations
## * 3000 females, 300 males
## * Selection only on the male side
## * One fully recessive lethal ~ 5% starting frequency
## * No carrier x carrier matings


library(AlphaSimR)


source("R/simulation_functions.R")


## Read command line arguments
## * selection_rule -- tells us what mating rule to use
## * lethal_is -- tells us whether to use balancing selection or not
## * outfile for populations
## * outfile for results


args <- commandArgs(trailingOnly = TRUE)

selection_rule <- args[1]
lethal_is <- args[2]
n_top_exempt <- as.numeric(args[3])
outfile_populations <- args[4]
outfile_results <- args[5]

print(selection_rule)
print(lethal_is)
print(n_top_exempt)
print(outfile_populations)
print(outfile_results)


## Set up simulation

setup <- make_simulation(n_ind = 6000,
                         n_chr = 10,
                         n_qtl = 10,
                         n_snp = 400,
                         h2 = 0.3)

founders <- setup$founderpop
simparam <- setup$simparam


## Pick a lethal variant

chosen_lethal <- pick_lethal(founders,
                             lethal_is = lethal_is,
                             simparam)
lethal_ix <- chosen_lethal$lethal_ix
other_snp_ix <- chosen_lethal$other_snp_ix
founder_lethal_frequency <- chosen_lethal$founder_lethal_frequency


## Breeding simulation

n_gen <- 20

generations <- list(n_gen)

generations[[1]] <- founders

for (gen_ix in 2:n_gen) {
    
    if (selection_rule == "avoid_carrier_x_carrier") {
    
        generations[[gen_ix]] <-
            breed_avoiding_carrier_x_carrier(generations[[gen_ix - 1]],
                                             lethal_ix,
                                             lethal_is = lethal_is,
                                             n_sires = 300,
                                             simparam = simparam)
        
    } else if (selection_rule == "unknown_lethal") {
        
        generations[[gen_ix]] <-
            breed_unknown_lethal(generations[[gen_ix - 1]],
                                 lethal_ix,
                                 lethal_is = lethal_is,
                                 n_sires = 300,
                                 simparam = simparam)
           
    } else if (selection_rule == "selection_against") {
        
        generations[[gen_ix]] <-
            breed_against_lethal(generations[[gen_ix - 1]],
                                 lethal_ix,
                                 lethal_is = lethal_is,
                                 n_sires = 300,
                                 n_top_exempt = n_top_exempt,
                                 simparam = simparam)
    }
}


## Create simulation object

carrier_status <- lapply(generations,
                         carrier_test,
                         lethal_ix = lethal_ix,
                         lethal_is = lethal_is,
                         simparam)

carriers <- get_carriers(carrier_status)

stats <- get_stats(generations)

if (lethal_is == "snp") {

    simulation_results <- list(simparam = simparam,
                               carrier_status = carrier_status,
                               carriers = carriers,
                               stats = stats,
                               lethal_ix = lethal_ix,
                               other_snp_ix = other_snp_ix)

} else if (lethal_is == "qtl") {

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
}

saveRDS(generations,
        file = outfile_populations)

saveRDS(simulation_results, 
        file = outfile_results)



