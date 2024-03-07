# cc1h_ms

## BatchProcessMinian.py
This is a script version of MiniAn pipeline, which allows us to process multiple recording sessions at one run, with most parameters remaining the same as recommended or default. We found these parameters doing OK with our data.
We add a little part to process the zarr files. Since we are a little more familiar with MATLAB, we have to convert the data into .mat and use MATLAB to process.
## GenerateNeuTraceMat.m
We then use the footprint, which is A.zarr and we convert it into A.mat, to 
