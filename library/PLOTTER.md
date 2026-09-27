# Plotter Guide

A successful experiment result can be visualized using a generalized plotting script `plotter.py`. It is recommended to generate figures and GIFs using the utilities in the script for any possible type of experiment that SpaceNet and its lunar extension allows.

## Assumptions

Before the plotting is implemented, there are 2 assumptions made by the plotter script

- Experiment output exists within the `output` folder with sufficient non-empty folders i.e. `experiment_info`, `node_indices`, `terrestrial_info`, `satellites_orbits`, `topology_graph`, `connectivity`, `routing` and `optimal_routes`.

- Correct selection of plotter variables as described [here](#definitions).

## Definitions

The experiment output can be referenced simply by inputting the output folder name to `outputfolder_name` variable. Although that is enough for the plotter script to search and determine all the necessary data for basic plotting, there are various ways the script can be tuned to provide a more flexible plotting mechanism. There exists plotter variables to professionally plot the results or use the script for visual debugging!

There are essentially 14 plotter variables excluding the `outputfolder_name` as follows:

- `plot_GS` -> Plots all the groudn stations
- `plot_only_optimal` -> Plots only the optimal path and no more details
- `show_optimal` -> Plots the  optimal path
- `plot_in_3D` -> 2D/3D plotting
- `plot_debug` -> Debugging mode
- `plot_optimal_orbits` -> Plots the orbit corresponding to the optimal paths
- `plot_all_isl` -> Plots all the ISLs of the constellation
- `make_gif` -> Makes the GIF
- `lon0_3d` -> Center longitude for 3D view
- `lat0_3d` -> Center latitude for 3D view
- `ll` -> scaled lower left point of the view
- `ur` -> scaled upper right point of the view
- `time_index` -> Time index for specific time instance debugging/plotting
- `gif_name` -> Name of the GIF
- `ref` -> Reference satellite ID for debugging (only runs when plot_debug=True)

<!-- ## Workflow

Draw a workflow of how the plotter workss!! -->

## Under-development
 1. $$\color{red}\text{Congestion-levels visualization}$$