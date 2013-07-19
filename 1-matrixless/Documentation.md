1) Matrixless Gauss-Seidel Solver
========================================

In the first installment of this series, we will have a look at a minimalistic implementation of a fluid solver on a MAC grid. It already provides a basic framework we're going to build on in the future; most of the structure we need is in place, and it is only the method implementations we're going to modify in future versions.

To give a small overview, the layout of one iteration of our solver (and most Eulerian solvers) is as follows:

- Add sources (inflows, body forces)
- Force incompressibility
  - Build pressure right-hand side
  - Build pressure matrix
  - Solve for pressure
  - Apply pressure
- Advect

Some of the items in the above list, such as inflows, are fairly straightforward and will stay almost the same for all of our implementations. Others, however, can be implemented in many different ways - some more accurate than others. This gives great flexibility in implementing the solver, as we have individual components which we can mix and match to get the speed and accuracy we want.

For this implementation, less accurate methods have been chosen for simplicity's sake, but we will slowly improve on them in the next versions of the solver.

#### Adding Sources

We only add inflows for this implementation; no body forces are applied. For inflows, we simply set the fluid quantities inside a specified rectangle to a fixed value.

#### Forcing Incompressibility

The pressure right hand side is quite standard in this case; it is simply the negative divergence of the velocity field, which can be computed easily from the MAC grid.

As you can probably guess from the title, we skip explicitly building the matrix for this implementation, since it has a very simple structure. Instead, we will simply construct the matrix implicitly in the pressure solve, which keeps in line with most of the "simple fluid" literature out there[3][4].

For the pressure solve itself, we will use a Gauss-Seidel solver [5], which is a very simple iterative solver for linear systems of equations. It doesn't show excellent convergence, but it is simple to implement and requires little memory.

Finally, applying the pressure is performed by subtracting the pressure gradient from the velocity field, which can be easily performed on the MAC grid.

#### Advection

A popular method of handling advection is the Semi-Lagrangian scheme. The idea of this approach is, given a point inside a grid cell, to determine the location this point would have been one timestep ago - that is, we trace the point backwards on the velocity field. Then we interpolate the fluid quantities at the determined location and set the result as the new value of the grid cell.

The Semi-Lagrangian scheme has two components we can fiddle with; tracing the velocity field and interpolating from the grid. The former falls in the category of numerical integration in a velocity field. We will pick the simplest available method, forward Euler, for this implementation. The latter is a simple interpolation problem. There are many interpolation methods we can pick from, but to keep things simple, we will use linear interpolation for now.

Both Euler and linear interpolation are not optimal, since they quickly cause loss of detail in the simulation - we will improve on that later.


### References and recommended reading:

  1. Robert Bridson - Fluid Simulation for Computer Graphics
       - Pages 3 - 50
  2. David Cline et al. - Fluid Flow for the Rest of Us
       - http://people.sc.fsu.edu/~jburkardt/pdf/fluid_flow_for_the_rest_of_us.pdf
  3. Jos Stam - Real-Time Fluid Dynamics for Games
       - http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
  4. Mark J. Harris - Fast Fluid Dynamics Simulation on the GPU
       - http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html
  5. http://en.wikipedia.org/wiki/Gauss-Seidel