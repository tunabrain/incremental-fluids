![Fluid](https://raw.github.com/tunabrain/incremental-fluids/master/Header.png)

Incremental fluids
==================

The purpose of this project is to provide simple, easy to understand fluid solver implementations in C++, together with code documentation, algorithm explanation and recommended reading. It is meant for people with beginner to intermediate knowledge of computational fluid dynamics looking for working reference implementations to run and study.

This project closely follows Robert Bridson's book, "Fluid Simulation for Computer Graphics", and implements a selection of the methods explained in the book. Ideally, you have a copy of the book sitting on your shelf, which will make it a lot easier to follow along with the code.

The solvers in this project come in a large variety, ranging from minimalistic to complex. All solvers are Eulerian in nature and run on a staggered Marker-and-Cell grid.

The different solvers are sorted into subfolders marked with a number and a short description. Each folder contains a small markdown file explaining the basic ideas behind the code and provides a list of recommended literature to read.

The number of the solver defines a progression - codes with higher number build on codes with lower number, either adding on features or replacing methods with better ones. The basic classes and concepts, however, always stay the same - ideally, you just start with the simplest solver and work your way through to understand and visually confirm the difference in simulation. Code is only explained once when it is introduced and not in any of the successing solvers to avoid clutter.

All solvers are single-file and require no external libraries apart from lodepng to save individual frames. Compilation should be straightforward on all platforms.

If you want, you can also check out a couple of videos rendered with code from this project (just ported to the GPU):

 - http://www.youtube.com/watch?v=D5uqoB13UBM
 - http://www.youtube.com/watch?v=SzqYnjIR4n0
 - http://www.youtube.com/watch?v=e_S-a5VNXbg
 - http://www.youtube.com/watch?v=vMB8elqhum0
 