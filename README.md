## OpenNucleome

OpenNucleome is an open-source software designed for conducting molecular dynamics (MD) simulations of the human nucleus. This software streamlines the process of setting up whole nucleus simulations through just a few lines of Python scripting. OpenNucleome can unveil intricate, high-resolution structural and dynamic chromosome arrangements at a 100 KB resolution. It empowers researchers to track the kinetics of condensate formation and fusion while also exploring the influence of chemical modifications on condensate stability. Furthermore, it facilitates the examination of nuclear envelope deformation's impact on genome organization. The software's modular architecture enhances its adaptability and extensibility. Leveraging the power of OpenMM, a GPU-accelerated MD engine, OpenNucleome ensures efficient simulations.

<img src="./images/Figure1.png" width="1000px"><img>

## Citation
Please cite the following paper if you use opennucleome package: 

    "OpenNucleome for high resolution nuclear structural and dynamical modeling", doi: XXX 

## Environment

We recommend using openmm 7.5.1 together with openNucleome, which can be installed with the following command: 

```
conda install -c conda-forge openmm=7.5.1
```

Other required packages: numpy, pandas, scikit-learn, MDAnalysis, mdtraj.

## Tutorials 

Detailed instructions for performing simulations of a human nucleus can be found at: tutorials/HFF_100KB. We also provide a tutorial for parameter optimization to develop new nucleus models at: tutorials/optimization.

## Manual

Detailed documentation of the source code can be found at: https://zhanggroup-mitchemistry.github.io/OpenNucleome/

