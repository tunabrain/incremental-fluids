2) Better advection
========================================

Normally, when we compute a numerical solution to an analytical problem, the computed solution will not be accurate - it contains errors.

These can be caused by the limited precision of the machine the code runs on, by the discretization (when we compute a finite number of grid cells, not a continuous domain) or by the numerical algorithms used  (such as the approximate Gauss-Seidel solver).

These errors are objectionable enough when we compute a single solution based on some accurate input data; however, they can become disastrous if we have to compute a series of solutions, with the computation of each (inaccurate) solution using the (also inaccurate) solution before it as input. Errors accumulate with each successive computation and can destroy the accuracy of the solutions. Unfortunately, this is exactly what our fluid solver is doing - the result of each iteration is based on the results of the last iteration, and errors will accumulate over time.

If the accumulated errors become unbounded, the computation becomes unstable and the solution worthless - this is commonly referred to as 'blow-up'. For the last version of our fluid solver, this fortunately cannot happen unless we get too extreme with the timestep, as most numerical methods we used are unconditionally stable - that is, they will never suffer from blow-up.

However, we still accumulate error. It won't destroy the simulation, but it will show up as inaccuracies in the solution - in this case, the errors will become visible as _numerical dissipation_. Essentially, numerical dissipation is a side-effect that causes the simulation to "smooth" out interesting features in the fluid flow, not unlike running a blur operator on all fluid quantities after each iteration.

This causes the density distribution to become blurry over time, destroying important visual cues; but what is even worse is that the velocities will be smoothed out as well, destroying the nice small-scale vortices and swirly flow that make fluids interesting.

The main culprit here is our simplistic advection routine - linear interpolation and forward Euler are not particularly good at preventing dissipation. Luckily, interpolation and numerical integration are both well explored topics and there are a lot of alternatives for us to choose from.

To improve interpolation, we will replace the linear interpolation with a cubic Catmull-Rom interpolation spline[2]. This is a bit of an arbitrary choice - any other interpolation spline of any order would work too - but cubic splines already make a large difference in accuracy without being overly complex to implement.

To improve interpolation, forward Euler is replaced with a third-order Runge-Kutta method. There is actually a large family of Runge-Kutta integration methods[3] of which the implemented one is only one of many members, but we won't delve too much into that. Again, it greatly improves integration accuracy without being overly complex.

One of the things that should be noted is that the Catmull-Rom spline can lead to oscillations. One of the main conditions of the method is that the function which is being interpolated is smooth. When this is not the case, such as at borders of the fluid domain or near sharp details in the fluid flow, the interpolation spline can over- or undershoot and lead to curious stripes in the fluid. For now, this is something we have to deal with, but we will be able to improve this situation once we introduce FLIP. One of the consequences of this is that we will have to slightly tweak our inflows - instead of setting fluid quantities to a specific value inside a rectangle, we place smooth "blobs" instead to ensure the fluid quantities near the inflow stay smooth to avoid oscillation. 

### References and recommended reading:

  1. Robert Bridson - Fluid Simulation for Computer Graphics
       - Chapter 3
       - Appendix A.2
  2. http://en.wikipedia.org/wiki/Catmull-Rom#Interpolation_on_the_unit_interval_without_exact_derivatives
  3. http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods