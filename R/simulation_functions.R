
## Helper functions


## Performs carrier test for a population, returning genotypes at the
## lethal locus
##
## Parameters:
## * pop -- population object
## * lethal_ix -- number of lethal in SNP or QTL genotypes
## * lethal_is -- indicator for pleiotropy, either "snp" (neutral variant) or 
##                "qtl" for pleiotropy
## * simparam -- simulation parameter object

carrier_test <- function(pop,
                         lethal_ix,
                         lethal_is,
                         simparam) {
    
    if (lethal_is == "snp") {
        result <- pullSnpGeno(pop,
                              simParam = simparam)[, lethal_ix]
    } else if (lethal_is == "qtl") {
        result <- pullQtlGeno(pop,
                              simParam = simparam,)[, lethal_ix]
    }
    
    result
}


## Picks a lethal variant from the founder population
##
## Parameters:
## * founderpop -- founder population object
## * lethal_is -- indicator for pleiotropy, either "snp" (neutral variant) or 
##                "qtl" for pleiotropy
## * simparam -- simulation parameter object

pick_lethal <- function(founderpop,
                        lethal_is,
                        simparam) {

    if (lethal_is == "snp") {
        founder_geno <- pullSnpGeno(founderpop,
                                    simParam = simparam)
    } else if (lethal_is == "qtl") {
        founder_geno <- pullQtlGeno(founderpop,
                                    simParam = simparam)
        effects <- simparam$traits[[1]]@addEff
    }
    
    f <- colSums(founder_geno)/(2 * nrow(founder_geno))
    candidates <- which(f > 0.04 & f < 0.06)
    
    if (lethal_is == "qtl") {
        candidate_effects <- effects[candidates]
        candidates <- candidates[candidate_effects > 0]
    }
    
    if (length(candidates) < 1) {
        stop("Found no candidate loci with appropriate frequency.")
    } else if (length(candidates) == 1) {
        lethal_ix <- candidates
    
    } else {
        lethal_ix <- sample(candidates, 1)
    }
    
    other_snp_ix <- setdiff(candidates, lethal_ix)
 
    list(lethal_ix = lethal_ix,
         other_snp_ix = other_snp_ix,
         founder_lethal_frequency = f[lethal_ix])
}

## Create a founder pouplation and a simulation parameters object

make_simulation <- function(n_ind,
                            n_chr,
                            h2) {
    
    founderpop <- runMacs(nInd = n_ind,
                          nChr = n_chr,
                          segSites = 500)
    
    simparam <- SimParam$new(founderpop)
    simparam$restrSegSites(minQtlPerChr = 100,
                           minSnpPerChr = 400,
                           overlap = FALSE)
    simparam$setGender("yes_sys")
    
    simparam$addTraitA(nQtlPerChr = 10)
    simparam$setVarE(h2 = h2)
    
    simparam$addSnpChip(nSnpPerChr = 400)
    
    founders <- newPop(founderpop,
                       simParam = simparam)
 
    list(founderpop = founders,
         simparam = simparam)   
}

## Simulate breeding that avoids carrier--carrier matings (but allows
## carrier--noncarrier matings on both sides)
##
## Parameters:
## * parent_generation -- population object with the parents
## * lethal_ix -- number of lethal in SNP or QTL genotypes
## * lethal_is -- indicator for pleiotropy, either "snp" (neutral variant) or 
##                "qtl" for pleiotropy
## * simparam -- simulation parameter object

breed_avoiding_carrier_x_carrier <- function(parent_generation,
                                             lethal_ix,
                                             lethal_is,
                                             simparam) {

    ## Take all females as dams and split them in carrier/noncarrier
    dams <- parent_generation[parent_generation@gender == "F"]

    dam_carrier_status <- carrier_test(dams,
                                       lethal_ix,
                                       lethal_is,
                                       simparam)

    carrier_dams <- dams[dam_carrier_status > 0]
    noncarrier_dams <- dams[dam_carrier_status == 0]

    n_carrier_dams <- sum(dam_carrier_status > 0)

    ## Select the top sires and split out noncarriers
    sires <- selectInd(pop = parent_generation,
                       trait = 1,
                       use = "pheno",
                       gender = "M",
                       nInd = 200,
                       simParam = simparam)

    sire_carrier_status <- carrier_test(sires,
                                        lethal_ix,
                                        lethal_is,
                                        simparam)

    noncarrier_sires <- sires[sire_carrier_status == 0]

    ## Create matings for noncarrier/noncarrier and carrier/noncarrier separately
    ## Splitting matings proportionally between the groups
    
    noncarrier_dam_x_any_sire <- randCross2(females = noncarrier_dams,
                                            males = sires,
                                            nProgeny = 1,
                                            nCrosses = 2 * noncarrier_dams@nInd,
                                            simParam = simparam)

    if (n_carrier_dams > 0) {
        carrier_dam_x_noncarrier_sire <- randCross2(females = carrier_dams,
                                                    males = noncarrier_sires,
                                                    nProgeny = 1,
                                                    nCrosses = 2 * carrier_dams@nInd,
                                                    simParam = simparam)
        
    offspring <- c(carrier_dam_x_noncarrier_sire,
                   noncarrier_dam_x_any_sire)
    } else {
        offspring <- noncarrier_dam_x_any_sire
    }
    
    offspring
}



###########################

## Functions for analysing simulation results


## Get genetic trends

get_stats <- function(generations) data.frame(generation = 1:length(generations),
                                          mean_g = unlist(lapply(generations, meanG)),
                                          var_g = unlist(lapply(generations, varG)),
                                          stringsAsFactors = FALSE)

## Get carrier numbesr for lethal

get_carriers <- function(carrier_status) {
    carriers <- unlist(lapply(carrier_status,
                              function (x) sum(x > 0)))
    data.frame(generation = 1:length(carrier_status),
               carriers = carriers)
}