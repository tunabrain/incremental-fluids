7) Variable Density
========================================

One of the things we glossed over in the last version is that temperature and fluid concentration not only influence the force exerted by buoyancy, but also have an effect on the fluid density. The fluid density has so far been a global constant, but to achieve higher visual fidelity, we now want to introduce spatially varying density. This allows for more interesting effects, such as heavy fluid displacing lighter fluid with ease.

To do this, we first have to compute the fluid density at each grid cell. Since we only need the density at the U/V velocity locations, we introduce two additional grids - one for the U grid cells and one for the V grid cells. The fluid density is then computed as a function of temperature relative to ambient temperature and fluid concentration, following Bridson's book[1]. These fluid densities can then be used for the construction of the pressure matrix and the pressure update, replacing the constant density from before. See Bridson's book[1] for details.

### References and recommended reading:

  1. Robert Bridson - Fluid Simulation for Computer Graphics
       - Chapter 5.3
       - Chapter 4.5