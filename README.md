# kidney_glom_studio

This repository provides a model for the HuBMAP kidney functional tissue unit (FTU) and a GUI to run a 2-D simulation and visualize results.

## Dependencies
* C++ compiler with OpenMP support.
* Python and a few modules that aren't provided in the standard Python library. We typically recommend installing the Anaconda Python distribution. However, it is quite large, so if you install another, leaner Python distribution, you will likely need to:
```
pip install matplotlib
pip install scipy
pip install PyQt5
```
## Usage
```
git clone <this repo>
# in kidney_glom_studio/model, compile the model:
make
# and copy the executable over to the studio directory:
cp mymodel ../studio

# in kidney_glom_studio/studio/config, unzip 2 data files
unzip grad_x_pbm_875x750.zip 
unzip grad_y_pbm_875x750.zip 

# Run the studio GUI:
cd ..
python studio.py
# on Run tab, click Run Simulation button
# on Plot tab: select "pbm_gbm_distance" in dropbox widget next to "Substrates", click Play button.
# Uncheck "Substrates" checkbox and click Play button.
```
