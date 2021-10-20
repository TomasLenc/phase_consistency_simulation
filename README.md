Simulation of how phase consistency affects FFT and ITPC using stimuli from Lenc et al. (2020, CeCorCom). 

To execute the simulation, set your current matlab working directory to `code/` and run the script `code/run_simulation.m`. All the parameters are set directly in this script (and feel free to play around with them).  

The output of the simulation is saved to a folder `figures/` created by the script if it doesn't exist. To open the `.fig` file in matlab, make sure you have `lib/panel.m` on your path and you have run the command below at least once in your matlab session: 
```
panel()
```
Otherwise you'll get an error when trying to open the figure and matlab will not want to close the figure window. If this happens, you can call 
```
close all force
```

Most of the work under the hood is done through the `subject` class which you can find in `lib/@subject`. An example instance of this class is saved with the figure which shows the results of the simulation. This way, one can keep track of parameters that were used to generate each figure. 