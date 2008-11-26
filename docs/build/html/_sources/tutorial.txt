Tutorial
========

This is where the tutorial will live..

Open the directory TOPKAPI_Example/. It contains three sub-directories that
correspond to the three steps of the process presented in the introduction
(create_the_parameter_files, run_the_model and analyse_the_results).  In
addition there is a small file called run_example.py which is a Python command
file, designed to run the tutorial in one action.  The following discussion
takes the reader step-by-step through the process to provide a better idea of
the structure of the processes.

Open a Python shell
Open a command shell for Python. 
If you installed IPython on Windows then click on Start->All programs->IPython.

Run the model
Type the following commands in the Python Shell (the shell commands are
hereafter highlighted in grey).
Import the package TOPKAPI
In [1]:  import TOPKAPI 
In the Python shell, navigate to the TOPKAPI_Example directory by copying the
address in the Windows folder. Typically this can be achieved using the `cd`
command:
In [2]: cd C:\ … \TOPKAPI_Example
Then run the model by typing (or copying without a leading space):
ln [3]: TOPKAPI.model.run(‘run_the_model/TOPKAPI.ini’)
The simulation starts and will last around 3 to 5 minutes depending on your
machine. Some comments are written and the number of time steps (to 169 of 6
hours each) scrolls followed by ‘***** THE END *****’.  Leave the iPython
window open - we will use it shortly.

At the end of the simulation, a folder results/ should have been created in
the run_the_model folder, containing an HDF5 file called 
Example_simulation_results.h5; to open it you will have to associate the file
(if it’s the first *.h5 file you’ve opened) with HDF_Explorer. This file is
composed of three groups: Channel, Overland and Soil. The group Channel
contains the table Qc_out that is the channel flow of each cell (in m3/s)
[to view a file’s content, double click on the little green square under the
yellow label in the left panel of the window - see Figure 7.]   The group
Overland contains the table V_o that is the volume of the overland reservoir
(in m3). The group Soil contains the table V_s that is the volume of the soil
reservoir (in m3). For the three tables Qc_out, V_o and V_s, the 3563 columns
are the cells and the 170 rows are the time steps.
