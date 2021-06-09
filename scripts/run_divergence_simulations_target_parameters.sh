#!/bin/bash

if [ ! -d simulations/divergence_simulations/target_parameters ]; then
    mkdir simulations/divergence_simulations/target_parameters
fi

for REP in {1..100}; do

    Rscript R/divergence_simulation_target_parameters.R \
        avoid_carrier_x_carrier \
        qtl \
        0 \
        simulations/divergence_simulations/target_parameters/avoid_carrier_x_carrier_qtl_populations_$REP.Rds \
        simulations/divergence_simulations/target_parameters/avoid_carrier_x_carrier_qtl_results_$REP.Rds 
                
done

