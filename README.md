# RMM4PYTHIA
PYTHIA8 output for RMM (rapidity-mass matrix)

This program generates Pythia8 events and fills rapidity-mass matrix as ROOT tree. The output file is in the "out" directory.
Pythia8 input parameters are in "tev14_SM_A14_NNPDF23LO.py".
 
To compile and execute this program, define $ROOTSYS (ROOT program), $FASTJET (fastjet library) and $PYTHIA (Pythia8 program).
Then compile and run using 


```
./A_RUN_SM 
```

It fills ROOT file in "out". Then validate the RMM using the PyROOT script (validate_data.py). This script removes 9 invariant masses which could be your signal region, and makes the image with average values of RMM.

S.Chekanov (ANL)    
