.. OpenNucleome documentation master file, created by
   sphinx-quickstart on Mon Apr 10 12:26:19 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to openNucleome's documentation!
========================================

OpenNucleome is an open-source software designed for conducting molecular dynamics (MD) simulations of the human nucleus. This software streamlines the process of setting up whole nucleus simulations through just a few lines of Python scripting. OpenNucleome can unveil intricate, high-resolution structural and dynamic chromosome arrangements at a 100 KB resolution. It empowers researchers to track the kinetics of condensate formation and fusion while also exploring the influence of chemical modifications on condensate stability. Furthermore, it facilitates the examination of nuclear envelope deformation's impact on genome organization. The software's modular architecture enhances its adaptability and extensibility. Leveraging the power of OpenMM, a GPU-accelerated MD engine, OpenNucleome ensures efficient simulations.

Environment
===========

We recommend using openmm 7.5.1 for using openNucleome, as openNucleome is built based on openmm 7.5.1.

Install openmm 7.5.1 with the following command: ``conda install -c conda-forge openmm=7.5.1``

Other required packages: numpy, pandas, scikit-learn, MDAnalysis, mdtraj.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   GettingStarted/Installation
   GettingStarted/Introduction

.. toctree::
   :maxdepth: 2
   :caption: OpenNucleome:

   OpenNucleome

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   Tutorials/LangevinSimulation
   Tutorials/BrownianSimulation
   Tutorials/NuclearDeformation
   Tutorials/ConfigsGeneration
   Tutorials/ConfigsRelaxation
   Tutorials/ConfigsSelection
   Tutorials/ComputeContactProb
   Tutorials/ChrChrOptimization
   Tutorials/ChrNLOptimization

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Citations
=========

"OpenNucleome for high resolution nuclear structural and dynamical modeling", eLife (in press), doi: 10.1101/2023.10.16.562451
