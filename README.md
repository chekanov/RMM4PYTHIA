# RMM4PYTHIA
PYTHIA8 output for RMM (rapidity-mass matrix)

This program generates Pythia8 events and fills rapidity-mass matrix as ROOT tree. The output file is in the "out" directory.
Pythia8 input parameters are in "tev14_SM_A14_NNPDF23LO.py".
 
To compile and execute this program, define $ROOTSYS (ROOT program), $FASTJET (fastjet library) and $PYTHIADIR (Pythia8 program).
If you do not have these libraries, one can use the singularity image file (/centos7hepsim.img) from https://atlaswww.hep.anl.gov/hepsim/doc/doku.php?id=hepsim:dev_singularity
and then run:

```
singularity exec centos7hepsim.img  bash -l
source /opt/hepsim.sh
```

Then compile main.cxx and run it using:  


```
./A_RUN_SM 
```

It fills a ROOT file in the "out" directory. It has some histogram and a ROOT tree with non-zero values of RMM.
Then validate the RMM using the PyROOT script (validate_data.py). This script removes 9 invariant masses which could be your signal region, and makes the image with average values of RMM.

S.Chekanov (ANL)    
