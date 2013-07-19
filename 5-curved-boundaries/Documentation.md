5) Curved Boundaries
========================================

One of the flaws of the previous solid boundary implementation is that we treat cells as either completely fluid or completely solid. Consequently, the solver only sees vertical or horizontal walls. This is fine when the actual solid surface is horizontal or vertical as well (an unrotated box, for example), but in any other case the fluid flow won't accurately follow the real boundary.

To remedy this, we implement the variational pressure solve introduced by Batty et al[2]. Instead of classifying cells as completely solid or completely fluid, we compute the fraction of the cell occupied by fluid. We then only have to modify the pressure matrix and right-hand side as outlined in Bridson's book[1].

The book proposes using supersampling to compute the volume fractions, but we opt to use a more accurate analytic method instead. See code for details.

### References and recommended reading:

  1. Robert Bridson - Fluid Simulation for Computer Graphics
       - Chapter 4.5
  2. Batty et al., A Fast Variational Framework for Accurate Solid-Fluid Coupling
       - http://hal.archives-ouvertes.fr/docs/00/38/47/25/PDF/variationalFluids.pdf