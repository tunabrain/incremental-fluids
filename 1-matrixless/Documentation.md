1) Matrixless Gauss-Seidel Solver
========================================

This implementation is about as simple as we will get in this project. It already provides a basic framework for fluids on a MAC grid, which we will build on for the next few solvers. It has most of the structure we need in place; only the method implementations themselves are what we are going to modify in the future.

The selection of methods used in this implementation is very similar to the well known Fluid Dynamics For Games [3], just on a MAC grid in this case; I recommend having at least read the paper or, even better, implemented the method.

The layout of our solver and, in fact, most Eulerian solvers, is as follows:

- Add sources (inflows, body forces)
- Force incompressibility
  - Build pressure right hand side
  - Build pressure matrix
  - Solve for pressure
  - Apply pressure
- Advect

Some of these points, such as inflows, are fairly straightforward and will stay the same for all of our implementations. Others, however, can be implemented in many different ways - some more accurate than others. This is one of the great things about Eulerian fluid simulation: The basic workings of the solver always stay the same - we just have individual components which we can mix and match to get the accuracy and speed we want.

For this implementation, less accurate methods have been chosen for simplicity's sake - but we will slowly improve on them in the next versions.

#### Adding Sources

We only add inflows for this version; no body forces are applied. Inflows are as simple as "set fluid quantity y to value x in specified region" and don't require much explanation.

#### Forcing Incompressibility

The pressure right hand side is quite standard in this case; it is simply the negative divergence of the velocity field, which can be computed easily from the MAC grid.

As you can probably guess from the title, we skip explicitly building the matrix altogether, since it has a very simple structure - namely, it is the five-point stencil or laplacian matrix. Each row of the matrix only has 5 non-zero entries, which refer to the current fluid cell itself and its four neighbours.

Instead, we will simply use the matrix implicitly in the pressure solve, which keeps in line with most of the "simple fluid" literature out there [3][4].

For the pressure solve itself, we will use a Gauss-Seidel solver [5], which is a very simple iterative solver for linear systems of equations. It doesn't show excellent convergence, but it is simple to implement and requires little memory.

Finally, applying the pressure is again quite standard, as we simply subtract the pressure gradient from the velocity field, which is easily done on the MAC grid.

#### Advection

A popular method of handling advection is semi-Lagrangian advection, in which we trace backwards in time for each fluid cell to determine the location of the new fluid value and interpolate from the grid at that point.

There are two components which we can fiddle with here:

  - Tracing back in time
  - Interpolating from the grid

The former falls in the category of numerical integration in a velocity field. We will pick the simplest available method, forward Euler, for this implementation.

The latter is a simple interpolation problem. There are many interpolation methods we can pick from, but to keep things simple, we will use linear interpolation for now.

Both Euler and linear interpolation are not very good, since they quickly cause loss of detail in the simulation - we will improve on that later.


### References and recommended reading:

  1. Robert Bridson - Fluid Simulation for Computer Graphics
     Pages 3 - 50
  2. David Cline et al. - Fluid Flow for the Rest of Us
     http://people.sc.fsu.edu/~jburkardt/pdf/fluid_flow_for_the_rest_of_us.pdf
  3. Jos Stam - Real-Time Fluid Dynamics for Games
     http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
  4. Mark J. Harris - Fast Fluid Dynamics Simulation on the GPU
     http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html
  5. http://en.wikipedia.org/wiki/Gauss-Seidel