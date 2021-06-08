#!/bin/bash

if [ ! -d simulations/divergence_simulations/target_parameters ]; then
    mkdir simulations/divergence_simulations/target_parameters
fi

for REP in {1..100}; do

    Rscript R/divergence_simulation.R \
        avoid_carrier_x_carrier \
        qtl \
        0 \
        simulations/divergence_simulations/target_parameters/avoid_carrier_x_carrier_qtl/populations_$REP.Rds \
        simulations/divergence_simulations/target_parameters/avoid_carrier_x_carrier_qtl/results_$REP.Rds 
                
done

