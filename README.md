# CompuCell3D Simulation of Cell Sorting and Chemotaxis 🧬

This repository contains a CompuCell3D project simulating the interplay between cell-cell adhesion and chemotaxis in driving cell sorting phenomena from various initial spatial configurations.


## 📝 Project Overview

This project uses the Cellular Potts Model (CPM) framework provided by [CompuCell3D](https://compucell3d.org/).

The primary objective is to observe and quantify how chemotactic signals, combined with differential adhesion, influence the self-organization and sorting of cells. The simulation is designed to be flexible, allowing for different initial cell layouts and tunable parameters to explore a range of biophysical scenarios.

### Key Features
- **Differential Adhesion:** Models cell sorting based on surface energy minimization between cell types.
- **Chemotaxis:** Implements directed cell migration along a chemical gradient.
- **Multiple Initial Configurations:** Includes steppables to initialize cells in various patterns (e.g., a central block, left/right halves, scattered).
- **Data Collection:** Automatically logs key metrics like the center of mass for each cell type and inter-cell contact areas over time.
- **Automated Simulation Runs:** Includes a helper script to run batches of simulations with varying parameters.
- **Data Analysis & Visualization:** Provides Python scripts to parse the output data and generate plots for analysis.

P.S. I'm sorry that the code is not that organized (e.g. the folder Code/Code which does NOT help with the clarity ofcourse). Furthermore there are some files through eachother which makes it harder to read. THere were 2 things that lead to this descision. The first one is that the modelling package called CompuCell3D does work in certain ways but does not in other ways. In this way, it worked. Secondary to this is the fact that I wanted to be done with my thesis. Thanks for checking it out!
