#!/bin/bash

# Emission of the source
gnuplot plot_emission.gplt

# Backflux
#plot_absorption_distr.gplt
gnuplot plot_exit.gplt

# Deposition
gnuplot plot_depos_r_phi.gplt
gnuplot plot_depos_z_phi.gplt
