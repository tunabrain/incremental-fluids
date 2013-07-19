8) Fluid Implicit Particle (FLIP)
========================================

In the final installment of this series, we will implement a method to almost completely eliminate numerical dissipation from our solver. Despite our efforts with higher-order interpolation, we still lose some of the small-scale detail in the simulation to dissipation. This can be remedied using the particle-based FLIP method introduced by Brackbill et al. [2]

One of the fundamental ideas of particle-in-cell methods is to combine particles and grids. Grids exceed at providing spatial derivatives, required for e.g. the pressure computation, whereas particles can naturally handle advection. Therefore, it makes sense to combine the strengths of both by computing advection using particles and everything else on the grid.

The way this works in FLIP is that we have a large number of particles distributed on the grid, carrying the values of all the fluid properties in the simulation. At the beginning of an iteration, the particles transfer their quantities onto the grid, after which the pressure computations etc. are carried out as usual. At the end of the iteration, the _change_ of the values on the grid is transferred to the particles. Then, the particles are advected in the fluid velocity field using standard numerical integration.

Why this almost completely eliminates dissipation is better explained in Bridson's book[1]. What we're more concerned with is the implementation of the method - the outline of FLIP is very simple, but in practice there are a few pitfalls one should look out for.

Initializing the particles is fairly simple, as we simply spawn a certain number (4 in this implementation) of particles inside each grid cell and randomly jitter them a bit from the grid position. Then their values are initialized by interpolating the values from the grid (using the existing interpolation routines).

Conceptually, for transferring the particle values onto the grid we sum the weighted particle contributions in each grid cell and divide by the sum of weights to get the final result. The particle weights stem from some filter function - for example, one could use a simple box filter, where for each cell we sum the values of all the particles in it and divide by the number of particles in the cell. While this is simple, it unfortunately leads to objectionable results, so we opt to use a linear hat filter instead; in essence, particles contribute less to grid cells that are further away than to ones that are closer. Either way, we limit the radius of influence of each particle to at most one grid cell - that is, it only ever influences the four closest grid cells around it.

To transfer the change in fluid quantities back to the particles, we keep a copy of each fluid quantity after the particle-to-grid transfer. After the rest of the fluid updates (pressure etc.), we can then subtract this copy from the freshly computed values to get the change in quantity. Then we can simply use the existing interpolation methods to get the change in quantity at each of the particle positions. Following Bridson's book[1], we also mix in a little bit of the standard PIC update to avoid noise.

Finally, advection can be performed using the existing third order Runge-Kutta implementation - only this time forward in time instead of backwards. Similar to before, if the newly computed position ends up inside a solid, we move it out of the solid to the nearest point in the fluid.

There is unfortunately one big problem one encounters when implementing this method. Namely, as the particles are moved around in the flow, they are typically no longer uniformly distributed. In some areas, particles are pushed together, resulting in overcrowded cells, whereas in other areas gaps may open up, leaving some cells with no particles near them. This can result in cells with undefined values after the particle-to-grid transfer.

Battling the first problem is as easy as deleting particles from a cell if it contains more than a certain number of particles (12 in this implementation). The latter problem is unfortunately not as easy to deal with, as in addition to spawning new particles in empty cells, we also have to somehow to replace undefined values with quantities that make sense.

Similar to our treatment of solid bodies, we can handle this as an extrapolation problem, and add some additional code to our extrapolation routine to also compute the values in empty cells from the closest fluid cells. After all the values are defined, we can then spawn new particles where necessary and initialize their values by interpolating from the grid. 

### References and recommended reading:

  1. Robert Bridson - Fluid Simulation for Computer Graphics
       - Chapter 10
  2. Brackbill et al: FLIP: A low-dissipation, particle-in-cell method for fluid flow
       - http://robotics.upenn.edu/~alla/courses/cis563/Flip.pdf