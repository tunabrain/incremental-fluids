6) Heat
========================================

So far, the only source of velocity in our simulation was due to inflows injecting fluid at a certain initial velocity. One way to make this more interesting is to include heat into the mix. Since hot air rises and cold air sinks, this will allow us to also introduce fluid motion due to changes in temperature.

To achieve this, we first have to introduce temperature as a fluid quantity. It is initialized and advected like all the other quantities and is centered within the fluid cells, just like the fluid concentration. However, unlike any of the other quantities, one additional step we have take into account is heat diffusion. Heat does not stay in place; hot areas exchange energy with cold areas, seeking to bring temperature to an equilibrium.

Heat diffusion is described by the heat equation[2]. There are several possible ways of solving it; for example, one could use forward Euler to directly have cells exchange heat with their neighbours. However, this approach is not very stable for large timesteps or diffusion rates, so we will use backward Euler instead - rather than directly computing the next solution from the current one, we solve for a future solution, such that if we went back in time, we would arrive at the current one. This can be formulated as a system of equations that shows the same matrix properties as the pressure matrix and can be solved using the exact same preconditioner and solver. This method is stable and not too computationally expensive, as it normally doesn't need many iterations to converge compared to the pressure solve.

One we have the heat distribution, we can finally introduce velocity due to buoyancy. For this purpose, we define an ambient temperature and add vertical force depending on the difference of fluid temperature to ambient temperature. Additionally, we also add in downward force due to fluid concentration - dense fluid is heavier than light fluid - which also allows for some interesting effects.

### References and recommended reading:

  1. Robert Bridson - Fluid Simulation for Computer Graphics
       - Chapter 5
  2. http://en.wikipedia.org/wiki/Heat_equation
  3. http://en.wikipedia.org/wiki/Backward_Euler_method