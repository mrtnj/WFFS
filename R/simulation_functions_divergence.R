## Create a founder pouplation and a simulation parameters object
## with two correlated breeding goal traits

make_simulation_divergence <- function(n_ind,
                                       n_chr,
                                       h2,
                                       corA) {
    
    founderpop <- runMacs(nInd = n_ind,
                          nChr = n_chr,
                          segSites = 500)
    
    simparam <- SimParam$new(founderpop)
    simparam$restrSegSites(minQtlPerChr = 100,
                           minSnpPerChr = 400,
                           overlap = FALSE)
    simparam$setGender("yes_sys")
    
    simparam$addTraitA(nQtlPerChr = 10,
                       mean = c(0, 0),
                       var = c(1, 1),
                       corA = matrix(c(1, corA,
                                       corA, 1),
                                     byrow = TRUE,
                                     ncol = 2,
                                     nrow = 2))
    simparam$setVarE(h2 = c(h2, h2))
    
    simparam$addSnpChip(nSnpPerChr = 400)
    
    founders <- newPop(founderpop,
                       simParam = simparam)
    
    list(founderpop = founders,
         simparam = simparam)   
}


## Simulate breeding with an unknown lethal. Individuals homozygous for the
## lethal are not used in breeding.
##
## Parameters:
## * parent_generation -- population object with the parents
## * lethal_ix -- number of lethal in SNP or QTL genotypes
## * lethal_is -- indicator for pleiotropy, either "snp" (neutral variant) or 
##                "qtl" for pleiotropy
## * n_sires -- number of sires used each generation
## * prop_goal2 -- proportion of sires selected based on breeding goal trait 2
## * simparam -- simulation parameter object

breed_unknown_lethal_divergence <- function(parent_generation,
                                            lethal_ix,
                                            lethal_is,
                                            n_sires,
                                            prop_goal2,
                                            simparam) {
    
    ## Exclude dams who are affected
    dams <- parent_generation[parent_generation@gender == "F"]
    
    dam_carrier_status <- carrier_test(dams,
                                       lethal_ix,
                                       lethal_is,
                                       simparam)
    
    nonaffected_dams <- dams[dam_carrier_status < 2]
    
    ## Exclude sires who are affected
    potential_sires <- parent_generation[parent_generation@gender == "M"]
    
    sire_carrier_status <- carrier_test(potential_sires,
                                        lethal_ix,
                                        lethal_is,
                                        simparam)
    
    nonaffected_sires <- potential_sires[sire_carrier_status < 2]
    
    
    sires_goal1 <- selectInd(pop = nonaffected_sires,
                       trait = 1,
                       use = "pheno",
                       gender = "M",
                       nInd = n_sires * (1 - prop_goal2),
                       simParam = simparam)
    
    sires_goal2 <- selectInd(pop = nonaffected_sires,
                             trait = 2,
                             use = "pheno",
                             gender = "M",
                             nInd = n_sires * prop_goal2,
                             simParam = simparam)

    sires <- c(sires_goal1,
               sires_goal2)    
    
    ## Create matings
    
    offspring <- randCross2(females = nonaffected_dams,
                            males = sires,
                            nProgeny = 1,
                            nCrosses = 6000,
                            simParam = simparam)
    
    offspring
}