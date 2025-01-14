This folder contains all scripts used to generate the results and figures for the paper.
The python scripts assume that you have installed Magritte using the last version of the comoving-approx branch.

The subfolder "stability_analysis" contains the julia code for the first few figures. (see readme in that folder for required depencies)
The subfolder "test_approx" contains the script used to test the difference between the approximate and complete implementation of the boundary conditions.

Finally, this folder contains the results of the application section of the paper.
Expected run order: import_and_reduce... -> compute_errors_... -> Plot_relative_...
The jupyter notebooks ending on "_more_refinement.ipynb" use the adaptive angular discretization, while the other notebooks use a uniform angular discretization.
