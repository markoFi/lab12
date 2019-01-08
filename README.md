# Lab 12
### Linear Schrödinger Equation with quadratic potential
----
Discretization of the linear Schrödinger equation
<p align="center">
<img src="stuffy_stuff/f1.png" width="400">
</p>
a la Crank-Nicolson results in a system of linear equations of the form
<p align="center">
<img src="stuffy_stuff/f2.png" width="120">
</p>
where
<p align="center">
<img src="stuffy_stuff/f3.png" width="320">
</p>
with
<p align="center">
<img src="stuffy_stuff/f4.png" width="250">
</p>

---
Implement the scheme above in order to solve the Schrödinger equation with  <img src="stuffy_stuff/f5.png" width="170">.

The initial condition is
<p align="center">
<img src="stuffy_stuff/f6.png" width="250">
</p>
where <img src="stuffy_stuff/f7.png" width="170">. Use <img src="stuffy_stuff/f8.png" width="120">. Place N=300 grid points in the interval x=-40 to x=40.

The analytical solution is
<p align="center">
<img src="stuffy_stuff/f9.png" width="550">
</p>
where <img src="stuffy_stuff/f10.png" width="120">. The period of the solution is <img src="stuffy_stuff/f11.png" width="60">.

---
This problem is taken from *Computational Methods for Physicists* by *S. Sirca* and *M. Horvat*.
