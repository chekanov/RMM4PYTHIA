# Pythia8 ATLAS MB Tune A14 with NNPDF23LO tune.
EventsNumber=1000
# sliming the particle record?
ApplyParticleSlim = off
#Collision settings
Random:setSeed = on
Random:seed = 0
Beams:idA = 2212
Beams:idB = 2212
Beams:eCM = 14000.
#physics processes
#HardQCD:all = on
HardQCD:all = off
ParticleDecays:limitTau0 = on
#Makes particles with c*tau>10 mm stable
ParticleDecays:tau0Max = 10
Tune:pp = 14 
Tune:ee = 7
PDF:pSet = LHAPDF6:NNPDF23_lo_as_0130_qed 
PDF:extrapolate = on

SpaceShower:rapidityOrder = on
SigmaProcess:alphaSvalue = 0.140
SpaceShower:pT0Ref = 1.56
SpaceShower:pTmaxFudge = 0.91
SpaceShower:pTdampFudge = 1.05
SpaceShower:alphaSvalue = 0.127
TimeShower:alphaSvalue = 0.127
BeamRemnants:primordialKThard = 1.88
MultipartonInteractions:pT0Ref = 2.09
MultipartonInteractions:alphaSvalue = 0.126
# old method: BeamRemnants:reconnectRange
ColourReconnection:reconnect=on
ColourReconnection:range=1.71
#
HardQCD:all = on
# 
PhaseSpace:pTHatMin = 120      # min pT
PhaseSpace:pTHatMax = 10000   # max pT
PhaseSpace:mHatMin = 250       # min mHat
PhaseSpace:mHatMax = 20000    # max mHat

# fill high-pT tail and add weights to events
#PhaseSpace:bias2Selection = on
#PhaseSpace:bias2SelectionPow = 5.0;

# all top decays
Top:all = on
# Higgs
HiggsSM:all=on
# W/Z
WeakSingleBoson:all = on
WeakDoubleBoson:all = on
WeakBosonAndParton:all = on
# prompt photons
PromptPhoton:all = on
