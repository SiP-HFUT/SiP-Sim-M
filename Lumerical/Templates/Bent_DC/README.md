### Objective
Provide a simple template for designing  **bent directional couplers**. 

The device consists of two straight-to-bent-to-straight (SBS) waveguides, with the upper one having a larger radius than the lower one, thus resulting in a gap between them.

### Usage
**Modify the parameters in *model*, not in individual objects below *model*** ;  Then run the simulation. Finally observe the S parameters of port1 and port2, and calculate $\kappa$ and $\tau$.

Note that *model* has been scripted such that the simulation environment can be automatically configured to fit the simulated device, 

Please note that the radii of the transition parts of the two SBSs are needed to be chosen properly to ensure that  the waveguides at the ends of the couplers have a sufficient gap (> 0.6 um for simulation purpose and > 2.5 um for designing real devices) to avoid residual coupling. Generally, setting the  upper SBS's transition radius $\times 1.5$ larger than that of the lower one can achieve this goal. Also note that thoses radius cannot be too small which could lead to losses. 



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



