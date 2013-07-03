Incremental fluids
==================

The purpose of this project is to provide simple, easy to understand fluid solver implementations in C++, together with code documentation, algorithm explanation and recommended reading. It is meant for people with beginner to intermediate knowledge of computational fluid dynamics looking for working reference implementations to run and study.

The solvers come in a large variety, ranging from inaccurate but simple to complex but fully featured. All solvers are Eulerian in nature and run on a staggered Marker-and-Cell grid.

The different solvers are sorted into subfolders marked with difficulty levels. Each folder contains a small markdown file explaining the basic ideas behind the code and provides a list of recommended literature to read.

Codes with higher difficulty build on codes with lower difficulty, either adding on features or replacing methods with better ones. The basic classes and concepts, however, always stay the same - ideally, you start with the simplest solver and work your way through to understand and visually confirm the difference in simulation.

All solvers are single-file and require no external libraries apart from lodepng to save individual frames. Compilation should be straightforward on all platforms.

It should be noted that this project will not explain the Navier-Stokes equations or provide in-depth descriptions of numerical methods - many other resources provide this kind of information. Ideally, you have a copy of "Fluid Simulation for Computer Graphics" by Robert Bridson sitting on your shelf, which will make it easy to follow along with the code.