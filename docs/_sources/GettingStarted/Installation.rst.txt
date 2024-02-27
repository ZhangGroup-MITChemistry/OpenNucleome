.. _installation:

============
Installation
============

The **openNucleome** library can be installed via `pip <https://pypi.org/>`_.

Install via pip
-----------------

The code below will install **openNucleome** from `PyPI <https://pypi.org/project/opennucleome/>`_.

.. code-block:: bash

    pip install opennucleome

The **openNucleome** library uses `OpenMM <http://openmm.org/>`_ API to run the whole nucleus simulations.
These requirements can be met by installing the OpenMM package from the `conda-forge channel <https://conda-forge.org/>`__:

.. code-block:: bash

    conda install -c conda-forge openmm=7.5.1
    
    
The following are libraries **required** for installing **openNucleome**:

- `Python <https://www.python.org/>`__ (>=3.6)
- `NumPy <https://www.numpy.org/>`__ (>=1.19.5)
- `SciPy <https://www.scipy.org/>`__ (>=1.5.4)
- `pandas <https://pandas.pydata.org/>`__ (>=1.1.5)
- `scikit-learn <https://scikit-learn.org/>`__ (>=0.24.2)
- `MDAnalysis <https://www.mdanalysis.org/>`__ (>=2.0.0)
- `OpenMM <https://openmm.org/>`__ (>=7.5.1)
