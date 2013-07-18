3) Conjugate Gradient Method
========================================

The next issue we are going to improve on is the method we have in place for solving the system of linear equations arising from the pressure equations - Gauss-Seidel.

While there are many iterative solvers for systems of equations out there, for our particular problem - large, sparse, symmetric positive definite matrices - there is almost no method as established as the conjugate gradient method (CG). CG will give us more control over the error in the solution and will also converge quicker, especially for larger systems of equations (if the grid resolution is increased). 

The underlying workings of the method - and _why_ it works - are not quite easy to understand, but if you're interested, there exists excellent literature about it[2]. Luckily though, implementing CG is quite easy regardless whether one understands it or not, simply by following an algorithm outline[3] - only basic vector and matrix operations are required.

There is one important tuning factor in the method which has a large impact on how quickly the solver converges, namely the preconditioner. Without going into much detail, there are certain properties about a matrix that influence how quickly CG converges; applying a preconditioner modifies the matrix (and subsequently the system of equations) such that it still describes the same solution, but tweaks the matrix properties to hopefully increase the convergence rate. For this implementation, we choose the modified incomplete Cholesky preconditioner (MIC) described in Bridson's book[1]. 

It is at this point that it makes sense to also replace our implicit matrix construction by an explicit one, and build the matrix in memory before the pressure solve. Since the pressure matrix is very sparse and symmetric, we can use a particular memory saving storage scheme outlined in Bridson's book[1]. This will also make later revisions easier, when we modify the systems of equations to include more complex effects.

### References and recommended reading:

  1. Robert Bridson - Fluid Simulation for Computer Graphics
       - Chapter 4
  2. Jonathan Shewchuck - An Introduction to the Conjugate Gradient Method Without the Agonizing Pain
       - http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
  3. http://en.wikipedia.org/wiki/Conjugate_gradients#The_resulting_algorithm