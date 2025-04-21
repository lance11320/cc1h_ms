# cc1h_ms

Usage (if need):
Arrange your data first and adjust the path in the scripts (The recommended way is: 

- animal/
  - session/
    - v4-miniscope/
    - behavior_cam/

Run BatchProcessMinian.py (in your minian environment);

Use ManuallyExtraction.m for manually add neurons;

Run FasterLowerMemCostGen.m to extract traces;

Use CheckNeuron.mlapp (Call appdesigner in MATLAB, and run it in appdesigner)

## BatchProcessMinian.py
This is a script version of MiniAn pipeline, which allows us to process multiple recording sessions at one run, with most parameters remaining the same as recommended or default. We found these parameters doing OK with our data.

We add a little part to process the zarr files. Since we are a little more familiar with MATLAB, we have to convert the data into .mat and use MATLAB to process.
## ManuallyExtraction.m
This script allows you to add neuron footprint according to max projection. MiniAn or CNMFe won't do anything for you after all.
## GenerateNeuTraceMat.m
We then use the footprint, which is A.zarr and we convert it into A.mat, to extract raw signal from preprocessed movies.

(We use FasterLowerMemCostGen.m in practice, since it is compatible with 64GB memory or less)

## CheckNeuron.mlapp
A MATLAB app for checking the neuron footprint and trace quality. The check results are marked and saved in quality metrics (0 for bad, 1 for good, 2 for uncertain, we use only 1 class at last).

## Other tips
About MiniAn installation:
1. Open your Anaconda prompt
2. conda create -n minian python=3.8 (note that it has to be 3.8)
3. conda activate minian
4. conda install mamba -c conda-forge
5. mamba install -y -c conda-forge minian
Now you can install it faster. We noticed that it is hard to use conda to install minian due to the slow solving process. Hoping this will make it easier for you to process your data.

## About FP_Process
This is the code for double-fiber photometry data processing (intend for activation at PVN while recording at mPFC neurons)

Use GetDFFdata or FPdata class to read doric V6 data if needed (however, it depends on how Doric arrange their data)

viewFPevent.mlapp is a GUI which might be used to view the doric V6 data and convert it into CSV for being compatible with old pipeline. We did not really use this in research but thought it might help.

## About AnalysisPipeline
Here we put our code of building analysis object, with some of the analysis method inside the class. 

## About Figure5_Code_YingZhou
Here we put the code for Figure5, which focus on registered neurons. These Python scripts include implementations of clustering algorithms and classifier construction. 

This part is contributed by YingZhou from Xiaoxuan Jia's lab.
