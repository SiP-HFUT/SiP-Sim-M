### Objective
Provide a simple template for calculating $\kappa$ and $\tau$ for a directional coupler made up by **a straight and a ring waveguide**. This is useful to calculate  $\kappa$ and $\tau$ of some micro-ring resonators.

### Usage
**Modify the parameters in *model*, not in individual objects below *model*** ;  Then run the simulation. Finally observe the S parameters of port1 and port2, and calculate $\kappa$ and $\tau$.

Note that *model* has been scripted such that the simulation environment can be automatically configured to fit the simulated device, 

### Saving the .fsp file with a name including all the ring parameters
Open the .fsp file, modify the device parameters, and then run *sim_file_generator.lsp*. Doing this will save a new .fsp file with a name including all the parameter Info. of the device

###  Manual checks are still needed before running the simulation
Note that although  the simulation environment is set to automatically fit the device,  basic manual checks are still needed before running the simulation, including；
- the wavelength range and the number of the wavelength points 
- selected modes of the ports
- directions of the ports
- source mode
- simulation time (if it is enough)
- index monitor, to check if the layout is meshed correctly, and  if the mesh resolution is enough to recover the original geometry of the layout 
- etc.






