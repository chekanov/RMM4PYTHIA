#export PYTHIADIR=/afs/cern.ch/user/z/zhangr/wis/bkg_study/pythia8309
#export FASTJET=/afs/cern.ch/user/z/zhangr/wis/bkg_study/fastjet-install
export PYTHIADIR=/cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/
export FASTJET=/cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/
source /afs/cern.ch/user/z/zhangr/wis/bkg_study/root/bin/thisroot.sh

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase # use your path
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

## More memory
#ulimit -S -s unlimited
lsetup "views LCG_102b x86_64-centos7-gcc11-opt"
##lsetup "root 6.10.02-x86_64-slc6-gcc62-opt"
