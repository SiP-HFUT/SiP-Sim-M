### Descriptions
Provide simple templates for FEM simulations of single and dual uniform waveguide (WG), for studying **mode characteristics** of various WGs.

### Usage
In both files, the setup script of the "model"  have been coded such that  **waveguide parameters are controlled in "model"**. Note that "model" has been coded such that the simulation region and the  oxide range have also been auto-modified to fit  the simulations. 

### Simulation Region Size
The span of the FED regions are set according to the FDTD rules posted in [Ansys Insight: Key FDTD simulation settings](https://forum.ansys.com/discussion/24201/ansys-insight-key-fdtd-simulation-settings%EF%BC%8C%E5%85%B6%E4%B8%AD%EF%BC%9A),  that 
>PML boundaries are a half wavelength away from the sides of any geometry objects in the simulation. “Wavelength” here refers to the longest wavelength in the source spectrum, taking into account the refractive index of the material between the object and the boundary.

Note that if the above rule also applied to FDE simulation here still remains to be confirmed.




