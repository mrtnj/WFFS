#!/bin/bash

if [ ! -d simulations ]; then
    mkdir simulations
fi

if [ ! -d simulations/simple_simulations ]; then
    mkdir simulations/simple_simulations
fi

for REP in 1..100; do

    Rscript R/simple_simulation.R $REP
    
done