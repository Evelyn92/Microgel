# Jupyter Notebooks



- main_project_real_data_007.ipynb

    This notebook includes some useful functions such as convex hull calculation, histogram calculation and visualization.

- Fancy_16FM.ipynb: 
  
  Generate fancy visualization of 3d point cloud for 16FM dataset

- Fancy_point_cloud.ipynb: 
  
  Generate fancy visualization of 3d point cloud for PAINT_DiffTemp dataset

- Fancy_real_data_selection.ipynb:
  
  generated 3D visualization for the selected point clouds according to the sample selection information.

- Sample_selection.ipynb:
  
  Present the density png plots from the MatLab-generated density folder.

- 2d_project_3d_visual.ipynb
  
  Generate 2d projection and 3d scatter plot of microgel samples.



# 1 Dataset: PAINT_DiffTemp

The shapes of microgels responds to the stimulate of temperatures.

Check the features in the paper:

[Microgel PAINT – nanoscopic polarity imaging of adaptive microgels without covalent labelling](https://pubs.rsc.org/en/content/articlehtml/2019/sc/c9sc03373d)

## 1.1 Dataset

### Collected real dataset

.\Microgel\4_Microgel_plotter_V1_density\PAINT_DiffTemp\Core-shell\\*C\Microgel_plotter_v2_solvatochromism\

where "\*C" means the dataset for the samples at temperature $*^\circ C$

### Generated dataset

.\Microgel\4_Microgel_plotter_V1_density\results_mtemp\results_mtemp_hilo

Notes:

- "sample_go1.ply" denotes the original real sample 1

- "sample_gg1.ply" refers to the corresponding generated sample 1.

- "sample_info1.ply" includes the information for normalizing "sample_go1.ply" and "sample_gg1.ply"

## 1.2 Matlab files

Generate the 2d density for generated "go" "gg" samples:

.\Microgel\4_Microgel_plotter_V1_density\visualize_hilo_density.m

Rank the correlation values to select the samples:

.\Microgel\4_Microgel_plotter_V1_density\SampleSelection.m



# 2 Dataset: 24

There are four different states of surface of microgels, check the details in the paper

 [Deformation of Microgels at Solid–Liquid Interfaces Visualized in Three-Dimension | Nano Letters](https://pubs.acs.org/doi/abs/10.1021/acs.nanolett.9b03688)



# 3 Dataset: 16FM

Microgel point cloud in 2d coordinate with lifetime values as the third dimension.

### Dataset folders

Collected real dataset:

.\Microgel\16-FM\lifetime

Generated dataset: 

.\Microgel\16-FM\results_lifetime

Notes: 

- "sample_go1.ply" denotes the original real sample 1

- "sample_gg1.ply" refers to the corresponding generated sample 1.

- "sample_info1.ply" includes the information for normalizing "sample_go1.ply" and "sample_gg1.ply"

```python

```

### Matlab files

Extract the samples from the entire slide and create 2d density plot for collected real data:

.\Microgel\16-FM\pc16fm_real.m

![](assets/2023-02-06-04-27-24-den-S1.png)

Create 2d density plot for "go""gg" files

.\Microgel\16-FM\pc16fm_gogg.m

![](assets/2023-02-06-04-23-14-den-S1gg.png)
