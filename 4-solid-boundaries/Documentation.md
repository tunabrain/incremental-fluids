4) Solid Boundaries
========================================

So far, the only boundary conditions taken into account in the solver are the walls of the simulation domain. While empty boxes are nice to play with, they're not particularly interesting - smooth walls tend not to generate complex fluid flow.

One way to remedy this is to introduce additional solid objects into the domain. For our purposes, these are stationary, unmovable objects, and we'll concentrate on simple shapes, such as boxes and spheres for now.

To account for additional boundaries, the solver has to be modified in a handful of places. First of all, the solver should only update cells which are not occupied by a solid - this includes the pressure update and advection. Additionally, the systems of equations has to be modified to include solid boundaries, which is performed closely following Bridson's book [1]. Finally, we have to deal with advection and interpolation.

One of the things that might happen during advection is that, after tracing back through the velocity field, the lookup point ends up inside a solid. This obviously doesn't make much sense, as there is no flow from the solid into the fluid. Since this is obviously an error, most likely caused by numerical inaccuracy or too large timesteps, we somehow have to fix this gracefully. One method of solving this problem is to simply move the lookup point from inside the solid onto the nearest point in the fluid domain, which is implemented in the accompanying code.

The next problem is that during interpolation, the interpolation method might require samples from inside the solid if the interpolation point is close to the surface of the solid. This is a problem, since fluid quantities are not really defined inside solids. While we can simply discard sample points that land inside solids, there is a more robust way of handling this problem. The solution implemented here is to extrapolate the closest fluid quantities into the solid so that interpolation can use all samples and returns consistent results.

For extrapolation, we deviate slightly from Bridson's book and use a fast marching method to solve the PDE on page 95 instead of using a fast sweeping method. The performance impact is negligible, but the results are more accurate in our case. 

### References and recommended reading:

  1. Robert Bridson - Fluid Simulation for Computer Graphics
       - Chapter 4
       - Chapter 6