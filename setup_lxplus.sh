export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase # use your path
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

## More memory
ulimit -S -s unlimited
lsetup "views LCG_102b x86_64-centos7-gcc11-opt"

export PYTHIADIR=/cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/
export FASTJET=/cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/
