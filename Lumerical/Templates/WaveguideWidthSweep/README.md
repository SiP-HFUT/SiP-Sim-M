## Objective
The files of the current folder provide a rapid tool for calculating and plotting the effective refractive index (neff) of the modes on silicon-on-insulator waveguides as a function of the waveguide width, using Lumerical. 

This functionality is useful for designing, for example, mode selective coupling devices.
## How to use
 1. create a new **Mode** project and then open and run *wg_2D_draw.lsf* script, which will generate the base of the simulation, including the structure and simulation region. 
 2. open the script *wg_2D_width_sweep.lsf*. In this script,  make some setup according to your application:
	 - the flag variable "*mode_po*" determines the calculated polarization of the modes  (0: only TE modes; 1: only TM modes; 2: both TE and TM modes). 
	 - the list variable "*widths*" determines the swept widths.
	 - the variable "*modes*" decides how many modes to be plotted in the following neff vs width plot. 
	 - the variable "*wavelength*" determines the $\lambda$ of the mode calculations.
 3. run the script *wg_2D_width_sweep.lsf* and observe the results.


**Welcome to report bugs to me at rcheng@hfut.edu.cn.**


