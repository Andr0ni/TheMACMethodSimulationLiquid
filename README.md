# TheMACMethodSimulationLiquid
When studying phenomena on the earth's surface using numerical methods for incompressible media, it is necessary to take into account the influence of free surfaces that limit the medium and the earth's gravity. In such cases, the marker and cell method (MAC method) can be used. This method is based on finite difference approximation,
applied to the complete Navier-Stokes equations. The main dependent variables are pressure and two velocity components (for the two-dimensional case).

# Procedure for calculating the method of markers and cells

1. At the initial moment of time, the velocity field is known (as a result of the previous calculation cycle, or from the given initial conditions). It is assumed that this velocity field is conservative. The coordinates of a set of elements are also known - markers they show which area is occupied by liquid and which is empty.

2. The pressure field is calculated in such a way as to ensure that the rate of change of the velocity divergence is also zero everywhere. To do this, the Poisson equation is solved by the iterative method.

3. Two components of acceleration are calculated; their products with increments of time per cycle then give changes in speed, which must be added to the old values.

4. Marker particles move in accordance with the velocity components in their vicinity.

5. Adjustments are made to allow marker particles to pass through cell boundaries. Every time, when this results in fluid entering a previously empty cell or emptying a cell that previously contained fluid, the necessary velocity changes are made.

6. The calculation is completed at the current moment of time and a transition is made to a new iteration in time. It should be noted that the marker particles introduced into this incompressible flow calculation are only used to indicate fluid configuration. They show which cells contain fluid and which cells lie along the free surface. Marker particles also serve as a visualization of the flow, with which you can observe the trajectories and the relative position of fluid elements, while they are not involved in the calculations.

# Examples
![image](https://github.com/Andr0ni/TheMACMethodSimulationLiquid/blob/main/Example.gif) 
