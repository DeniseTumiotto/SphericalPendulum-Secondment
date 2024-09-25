# Contractivity on Riemaniann Manifold - A test case
## Project collaboration MLU - NTNU
### Theory
The theoretical background for the study can be found in the [paper](https://www.aimsciences.org//article/doi/10.3934/jcd.2024002)

### Code

`a_main.m`

In the main script, the solution of the spherical pendulum is evaluated. Where the system is  
q'= ω × q  
ω'= ge<sub>3</sub> × q − cω  
where (q, ω) is a vector of length six, g is the gravity acceleration and c is a damping parameter.  
  
No input is requires and the output is saved in the folder `out/` as a txt file.
The solution is saved as a matrix: in each column there is the solution in the current time, such that the total number of columns is equal to the number of time steps done during the simulation.  
  
`readAll.m` or `readSol.m`  
  
In order to read the solution .txt file, the two reading functions are created. The first one accept up to one input specifying how many solutions to read (starting from the oldest). The second one accept up to one input, it will read one specified file (0 --> most recent, -1 --> second most recent, and so on...)  
The output can be a cell structure containing matrices of solutions or a solution. The form of the matrix is 6 x _#time steps_.  
  
`evalGeo.m`
  
The function requires as input two solutions.  
The output is a cell structure containing _n=#time steps_ geodesics connecting each point of the inputs at each time step.  

`plotTraj.m`, `plotGeo.m`, `plotEnergy.m`
  
Given the solutions and the geodesics, one can plot and visualized the trajectories, the geodesics and evaluate the energy of the system. 
  
`riemannDistance.m`
  
The input requires two solutions and (optionally) the geodesics.
The output is  a vector of the distance in the Riemaniann Manifold at each time step of the two solutions.  
