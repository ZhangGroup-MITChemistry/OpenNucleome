## openNucleome

OpenNucleome is an open-source software designed for conducting molecular dynamics (MD) simulations of the human nucleus. This software streamlines the process of setting up whole nucleus simulations through just a few lines of Python scripting. OpenNucleome can unveil intricate, high-resolution structural and dynamic chromosome arrangements at a 100 KB resolution. It empowers researchers to track the kinetics of condensate formation and fusion while also exploring the influence of chemical modifications on condensate stability. Furthermore, it facilitates the examination of nuclear envelope deformation's impact on genome organization. The software's modular architecture enhances its adaptability and extensibility. Leveraging the power of OpenMM, a GPU-accelerated MD engine, OpenNucleome ensures efficient simulations.

## Manual

The output html manual file is docs/index.html. 

The manual is also shown in: https://zhanggroup-mitchemistry.github.io/OpenNucleome/

Instructions for class methods and functions are also included as comments in the source code. 

## Environment

We recommend using openmm 7.5.1 for using openNucleome, as openNucleome is built based on openmm 7.5.1. 

Install openmm 7.5.1 with the following command: 

```
conda install -c conda-forge openmm=7.5.1
```

Other required packages: numpy, mdanalysis, mdtraj.

## Workflow for creating the whole nucleus model and efficiency across different platforms

<img src="./images/Figure1.png" width="2000px"><img>

## Examples

We provide an example in examples/HFF_100KB, and we also show how to calculate the contact probabilities and compare with experiments to obtain optimized potential parameters in ./analysis_code.
