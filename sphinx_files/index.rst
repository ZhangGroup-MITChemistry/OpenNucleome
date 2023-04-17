.. OpenNucleome documentation master file, created by
   sphinx-quickstart on Mon Apr 10 12:26:19 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OpenNucleome's documentation!
========================================

OpenNucleome supports OpenMM GPU-accelerated simulations of the three-dimensional architecture of genomes. It is flexible and allows adding the potential for chromosomes-chromosomes and chromosomes-nuclear bodies with custom values. The package dramatically simplifies the simulation setup: users only need a few lines of python code to carry out the genome simulations starting from initial configurations of 46 chromosomes, nucleoli, speckles, and lamina, which form the cell nuclei. The package is integrated with OpenMM, a GPU-accelerated MD simulation engine enabling efficient simulations. We provide tutorials in Jupyter notebooks to demonstrate the various capabilities. We anticipate OpenNucleome to significantly facilitate for simulating the three-dimensinoal architecture of genomes.

Environment
===========

We recommend using openmm 7.5.1 for using OpenNucleome, as OpenNucleome is built based on openmm 7.5.1.

Install openmm 7.5.1 with the following command: ``conda install -c conda-forge openmm=7.5.1``

Other required packages: numpy, mdanalysis, mdtraj.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   manual

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Citations
=========

We will add reference after the paper is formally online.
