## OpenNucleome

OpenNucleome is an open-source software designed for conducting molecular dynamics (MD) simulations of the human nucleus. This software streamlines the process of setting up whole nucleus simulations through just a few lines of Python scripting. OpenNucleome can unveil intricate, high-resolution structural and dynamic chromosome arrangements at a 100 KB resolution. It empowers researchers to track the kinetics of condensate formation and fusion while also exploring the influence of chemical modifications on condensate stability. Furthermore, it facilitates the examination of nuclear envelope deformation's impact on genome organization. The software's modular architecture enhances its adaptability and extensibility. Leveraging the power of OpenMM, a GPU-accelerated MD engine, OpenNucleome ensures efficient simulations.

<p align="center">
<img src="./images/intro_figure.png" width="500px"><img>
</p>

## Citation
Please cite the following paper if you use opennucleome package: 

    "OpenNucleome for high resolution nuclear structural and dynamical modeling", eLife (in press), https://doi.org/10.7554/eLife.93223.1

## Environment

We recommend using openmm 7.5.1 together with openNucleome, which can be installed with the following command: 

```
conda install -c conda-forge openmm=7.5.1
```

Other required packages: numpy, pandas, scipy, scikit-learn, MDAnalysis, mdtraj.

## Installation

The user can either download the package from github, or use pip install:

```
pip install opennucleome
```

## Tutorials 

- Setting up and Performing Nucleus Simulations (tutorials/HFF_100KB)

- Analyzing Simulation Trajectories and Computing Various Genome Wide Properties (tutorials/compute_contact_prob)

- Preparing Initial Configurations for Simulations (tutorials/init_config)

- Optimizing Model Parameters (tutorials/optimization)

## Manual

Detailed documentation of the source code can be found at: https://zhanggroup-mitchemistry.github.io/OpenNucleome/

