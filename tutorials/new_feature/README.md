## New feature

As an example, this folder provides tutorials for including new features into the existing model.

- chromosome.py, nucleolus.py, speckle.py, lamina.py: source code

- whole_nucleus_model.py: users should change some lines for adding new features, and a more comprehensive tutorials will be followed

- newfeature.py: save the potential between new beads and other components. Here, as an example, we add the potential between 1. new beads and new beads; 2. new beads and centromeres

- configs_generation: generate the configuration with new beads

Here, we show how to add new features.

 1. Create the configurations with new beads, a comprehensive tutorial is included in the folder "./configs_generation"

 2. Modify the functions in the script "whole_nucleus_model.py"
  2.1 ajsodg

