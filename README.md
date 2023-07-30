# GEO-Nav: A geometric dataset of voltage-gated sodium channels

This repository shares the output files of Chanalyzer and HOLE, together with the tools that were used to produced the results in [1]. The dataset and the ground truths behind the benchmark can be found at:

https://github.com/rea1991/GEO-Nav

### Before you start
This python implementation is tasked to analyze and compare the Nav channel retrieved by two computational methods (Chanalyzer and HOLE) for the models of the GEO-Nav dataset. The following external modules are required: `networkx` and `matplotlib`.

### Input files
It requires the centerlines (endowed with radius values of the maximal inscribed balls) of the Nav channel of the considered model retrieved by performing Chanalyzer and HOLE.

Example input files can be found in the `example` folder.
To run it, once you have downloaded the repository, open the terminal, reach the folder `GEO-Nav-methods-main` and execute the following command line `python GEO-Nav_evaluation_measures.py`
The results will appear in the folder `example/output`. To analyse output files 1 or 2 from HOLE, set `number` to be 1 or 2. 

### Parameters
- HOLE adopts a process of user-assisted cavity localization. GEO-Nav dataset contains two Nav channel retrieved by HOLE by providing as initial information the two initial points and directions retrieved as the entrances of the channel by Chanalyzer. Parameter `number` can be freely set as 1 or as 2 for choosing among the two Nav channels returned by HOLE.
- Parameter `sC` (analogously, `sH`) sets the number of points adopted for sampling the spheres representing the Nav channel returned by Chanalyzer (analogously, by HOLE).
- In order to run the program on models different than the one considered in the example: comment lines 150 and 164, uncomment lines 151 and 165, and choose the structure to be considered by setting the `model` parameter. An analogous suitable modification has to be performed to run the program on the toy models.  


### Output files
It returns on the terminal the following information about the centerline of the Nav channel of the considered model retrieved by the two considered computational methods (Chanalyzer and HOLE):
- the number of vertices,
- the length,
- the straightness,
- the values of the comparison measures `match` and `d_\rho`.

It returns four OFF files encoding point clouds representing:
- the centerline of the Nav channel of the considered model retrieved by Chanalyzer and colored in accordance with the radius values adopting the `coolwarm` colormap of Matplotlib,
- the Nav channel (as a collection of spheres) of the considered model retrieved by Chanalyzer and colored in accordance with the radius values adopting the `coolwarm` colormap of Matplotlib,
- the centerline of the Nav channel of the considered model retrieved by HOLE and colored in accordance with the radius values adopting the `coolwarm` colormap of Matplotlib,
- the Nav channel (as a collection of spheres) of the considered model retrieved by HOLE and colored in accordance with the radius values adopting the `coolwarm` colormap of Matplotlib.

It return one PNG file representing:
- the graphs of the radius functions of the centerlines of the considered model. For each structure, the centerline retrieved by Chanalyzer is depicted in blue while the one produced by HOLE is in orange. Vertical dashed lines denote the extrema of the interval in which the two centerlines are identified as matched.

Example output files can be found in the `example` folder.
For a better visual evaluation, the `example` folder also contains the OFF file representing the triangulated surface of the considered model.

### References
[1]   A. Raffo, U. Fugacci, "GEO-Nav: A geometric dataset of voltage-gated sodium channels", *Computers & Graphics*, 2023. DOI: 10.1016/j.cag.2023.06.023.
