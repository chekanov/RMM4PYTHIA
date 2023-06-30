// Fill RMM (rapidity-mass matrix) from Pythia8 
// S.Chekanov (ANL)
// 2023 (May) 

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom2.h>
#include <TGraph.h>
#include <TProfile2D.h>
#include <map>
#include <limits>       // std::numeric_limits
#pragma GCC diagnostic ignored "-pedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "LParticle.h"
#include "CParticle.h"
#define MAGENTA "\033[35m"      /* Magenta */
#define RESET   "\033[0m"
#define RED     "\033[31m"      /* Red */

#include <TMath.h>
const double kPI   = TMath::Pi();
const double k2PI  = 2*kPI;

#include "Pythia8/Pythia.h"
using namespace Pythia8;
using namespace fastjet;

// project event
float**  projectevent(const float  CMS, const int maxN, const int maxNumberTypes,
                      const vector<LParticle> missing,
                      const vector<LParticle> jets,
                      const vector<LParticle> bjets,
                      const vector<LParticle> muons,
                      const vector<LParticle> electrons,
                      const vector<LParticle> photons);

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
        std::stringstream ss(s);
        std::string item;
        while(std::getline(ss, item, delim)) {
                elems.push_back(item);
        }
        return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
        std::vector<std::string> elems;
        return split(s, delim, elems);
}


int main(int argc, char* argv[]) {


	// Check that correct number of command-line arguments
	if (argc != 5) {
		cerr << " Unexpected number of command-line arguments. \n You are"
		<< " Arguments: input, output, nr_events, seed (or -1 for CPU seed)\n"
		<< " Program stopped! " << endl;
		return 1;
	}


        const bool   debug=false;

	cout << "HepSim:  Pythia8 Input Configuration =" << argv[1] << endl;
	cout << "HepSim:  ProMC Output =" << argv[2] << endl;
	cout << "HepSim:  Input events =" << argv[3] << endl;
        cout << "HepSim:  Input seed =" << argv[4] << endl;


	string infile("-"), outfile("-"),  sevents("-"), strseed("-");
	infile = argv[1];
	outfile = argv[2];
	sevents = argv[3];
        strseed = argv[4];
	int Ntot = atoi(sevents.c_str());
        int Iseed = atoi( strseed.c_str());
        string sseed="Random:seed = "+strseed;
        if (Iseed<0)  cout << "HepSim:  CPU timestamp will be used for seed" << endl; 
        else          cout << "HepSim:  Pythia uses  =" << strseed << " as a seed " << endl;

	// Generator. Process selection. Tevatron initialization. Histogram.
	Pythia pythia;

	/////////// read config files ////////////////////
	string sets="";
	string sets1="";
	bool   apply_slim=true;

	vector<string> configs;
	string events;
	ifstream myfile;
	myfile.open(infile.c_str(), ios::in);
	if (!myfile) {
		cerr << "Can't open input file:  " << infile.c_str() << endl;
		exit(1);
	} else {
		string line;
		while(getline(myfile,line))
		{
			//the following line trims white space from the beginning of the string
			line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace))));
			if(line[0] == '#') continue;
			if (line.length()<3) continue;
			string tmp=string(line);
			// no empty spaces inside string
			std::string::iterator end_pos = std::remove(tmp.begin(), tmp.end(), ' ');
			tmp.erase(end_pos, tmp.end());
			bool special=false;
			int found1=tmp.find("EventsNumber");
		    if (found1!=(int)std::string::npos) {events=tmp; special=true;}
			int found2=tmp.find("ApplyParticleSlim=on");
			if (found2!=(int)std::string::npos) {apply_slim=true; special=true;}
			int found3=tmp.find("ApplyParticleSlim=off");
			if (found3!=(int)std::string::npos) {apply_slim=false; special=true;}
			if (!special)  {sets1=sets1+tmp+"; "; pythia.readString(line); }
			configs.push_back(line);
		}
		myfile.close();
		vector<string> readnum=split(events,'=');
		//Ntot= atoi(readnum[1].c_str());
		cout << "Reading events. " << events << " Total number is=" << Ntot<< endl;
		for (unsigned int i=0; i<configs.size(); i++) {
			cout << ".. input ="+configs[i] << endl;
			sets=sets+configs[i]+";";
		}
	} // end else


        if (Iseed>0) { 
           pythia.readString("Random:setSeed = on");
           pythia.readString(sseed.c_str());
        }

	pythia.init();

	pythia.settings.listChanged(); // Show changed settings
	double versionNumber = pythia.settings.parm("Pythia:versionNumber");
	pythia.particleData.listChanged(); // Show changed particle data
	std::stringstream s;
	s << versionNumber;
	string version=s.str();


        float CMenergy=(float)pythia.settings.parm("Beams:eCM");

	// fastjet
	const double ptLepton=30;
	const double ptJet=30;
	const double etaJet=2.4;
        const double Rsize = 0.4; // jet size 0.6	
        const double misRateEle=1;  // 1% 
        const double misRateMu=0.1; // 0.1% 
       // https://arxiv.org/pdf/1512.01094.pdf
        //1+0.003*pT  (~1000 GeV is about 4%)
        const double misRateB=1;   // 1% b-tag misidentification rate
        const double Rlep_iso=0.1;
        const double Rjet_iso=0.4;
        const double etaLep=3.0;
        const double angNorm=0.15;

        // https://arxiv.org/pdf/1512.01094.pdf
        //1+0.003*pT  (~1000 GeV is about 4%)

	cout << "min PT lepton=" << ptLepton << endl;
	cout << "min PT jet=" << ptJet << endl;
	cout << "max ETA jet=" << etaJet << endl;
	cout << "jet R=" << Rsize << endl;
	cout << "lepton isolation cone=" << Rlep_iso << endl;
        cout << "jet  dR isolation=" << Rjet_iso << endl;
        cout << "Mis-identification rates in %. Ele:" << misRateEle << " Mu:" << misRateMu << endl;
        cout << "mistag b-jet rate % :" << misRateB << endl;
        cout << "CM energy [GeV] :" << CMenergy << endl;

        const int maxNumber=10; // max number for each object (MET is not counted)
        const int maxTypes=5;   // max numbers of types (met not counted)
        string names[maxTypes+1] = {"MET","j", "b", "#mu", "e", "#gamma"};
        string names_debug[maxTypes+1] = {"e/met","j", "b", "m", "e", "g"};

        cout << "Project using max number of each object=" << maxNumber << endl;
        cout << "Number of particle types==" << maxTypes << endl;
        const int mSize=maxTypes*maxNumber+1;
        // non-zero in triangular matrix
        //int NonZero=((1+mSize)*mSize)/2;



        // jets
        Strategy strategy = fastjet::Best;
        JetDefinition jet_def(fastjet::antikt_algorithm, Rsize, strategy);
        //JetDefinition jet_def(fastjet::kt_algorithm, Rparam, strategy);
        //JetDefinition jet_def(fastjet::cambridge_algorithm, Rparam, strategy);
        //JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, strategy);

        TRandom2 misrate;

	// book a histogram, make sure that the output name is Analysis.root
	cout << "\n -> Output file is =" << outfile.c_str() << endl;
	TFile * RootFile = new TFile(outfile.c_str(), "RECREATE", "Histogram file");
	TH1D * h_debug = new TH1D("debug", "debug", 5, 0, 5.);
        TH1D * h_info = new TH1D("info", "version,CM,seed", 5, 0, 5.);
        TH1D * h_cross = new TH1D("cross", "cross,events,lumi", 10, 0, 10.);
	TH1D * h_weightsFine = new TH1D("weightsFine", "weights fine", 10000, 0, 0.001);
        TH1D * h_weightsCorse = new TH1D("weightsCorse", "weights corse", 100, 0, 5);
	TH1D * h_iso = new TH1D("iso_energy", "isolation fraction", 100, 0, 3.0);
        TH1D * h_met = new TH1D("met_pt", "MET",100,0,1000);
        TH1D * h_bquarks = new TH1D("bquarks", "Nr of b-quarks",20,0,20);
        TH1D * h_jetN = new TH1D("jetN", "Nr of light jets",20,0,20);
        TH1D * h_bjetN = new TH1D("bjetN", "Nr of BJet",20,0,20);
        TH1D * h_photonN = new TH1D("photonN", "Nr of photons",20,0,20);
        TH1D * h_muonN = new TH1D("muonN", "Nr of muons",20,0,20);
        TH1D * h_electronN = new TH1D("electronN", "Nr of electrons",20,0,20);
        TH1D * h_n_lepton = new TH1D("lepton_nr", "Nr of leptons",20,0,20);
        TH1D * h_pt_jet = new TH1D("jet_pt", "pt",100,0,1000);
        TH1D * h_pt_photon = new TH1D("photon_pt", "leading photon pT",100,0,1000);
        TH1D * h_pt_electron = new TH1D("electron_pt", "leading electron pT",100,0,1000);
        TH1D * h_pt_muon = new TH1D("muon_pt", "leading muon pT",100,0,1000);
        TH1D * h_pt_bjet = new TH1D("bjet_pt", "pt",100,0,1000);
        TH1D * h_eta_jet = new TH1D("jet_eta", "eta", 40, -10, 10);
        TH1D * h_eta_bjet = new TH1D("bjet_eta", "eta", 40, -10, 10);
        TH2D * h_proj = new TH2D("projection", "projection", mSize, 0, (double)(mSize), mSize, 0, (double)mSize);
        TH1D * h_dimensions = new TH1D("dimensions", "(1)maxNumber,(2)maxTypes, (3)mSize",5,0,5);
        TH1D * h_process = new TH1D("processID", "processID",5000,0,5000);
        TH1D * cpucores = new TH1D( "cpucores", "CPU cores", 5, 0, 5);
        cpucores->Fill(1,1);


        h_dimensions->Fill(1,(float)maxNumber);
        h_dimensions->Fill(2,(float)maxTypes);
        h_dimensions->Fill(3,(float)mSize);

        TProfile2D * h_prof = new TProfile2D("profile", "profile", mSize, 0, (double)(mSize), mSize, 0, (double)mSize, 0, 1000);


        TTree *  m_tree  = new TTree("inputNN","inputNN");
        m_tree->SetAutoFlush(100000);
        UInt_t m_id;
        std::vector<Double32_t> m_proj;
        std::vector<UInt_t> m_proj_index1;
        std::vector<UInt_t> m_proj_index2;
        UInt_t m_run;
        Long64_t m_event;
        m_tree->Branch("id",  &m_id);
        m_tree->Branch("proj",   &m_proj);
        m_tree->Branch("proj_index1",   &m_proj_index1);
        m_tree->Branch("proj_index2",   &m_proj_index2);
        m_tree->Branch("run",   &m_run);
        m_tree->Branch("event", &m_event);

        double mass_jj, mass_jb, mass_bb;
        m_tree->Branch("mass_jj",   &mass_jj);
        m_tree->Branch("mass_jb",   &mass_jb);
        m_tree->Branch("mass_bb",   &mass_bb);

        std::vector<string> Names1;
        Names1.push_back(names[0]);

        for (int h = 1; h < maxTypes+1; h++) {
                for (int i = 1; i <  maxNumber+1; i++) {
                        ostringstream ss;
                        ss << i;
                        Names1.push_back(names[h]+"_{"+ss.str()+"}");
                }
        }


        std::vector<string> Names2;
        for (unsigned int i=0; i<Names1.size(); i++) {
                cout << "Name=" << i << " " << Names1.at(i) << endl;
                Names2.push_back(Names1.at(i));
        }

        // for plotting
        std::reverse(Names1.begin(), Names1.end());

/*
// range extensions
ranges=[]
R=13156
deltaStart=246
for i in range(60):
    delta=deltaStart+5+i
    R=R+delta
    ranges.append(R)

print ranges
*/
	double mjjBins[] = {99,112,125,138,151,164,177,190, 203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156, 13407, 13659, 13912, 14166, 14421, 14677, 14934, 15192, 15451, 15711, 15972, 16234, 16497, 16761, 17026, 17292, 17559, 17827, 18096, 18366, 18637, 18909, 19182, 19456, 19731, 20007, 20284, 20562, 20841, 21121, 21402, 21684, 21967, 22251, 22536, 22822, 23109, 23397, 23686, 23976, 24267, 24559, 24852, 25146, 25441, 25737, 26034, 26332, 26631, 26931, 27232, 27534, 27837, 28141, 28446, 28752, 29059, 29367, 29676, 29986};  

	const int nBins=sizeof(mjjBins)/sizeof(double);
	double xbins[nBins];
        double xbins_tev[nBins];
        for (int j=0; j<nBins; j++){xbins[j]=mjjBins[j]; xbins_tev[j]=0.001*mjjBins[j]; };

        TH1D * binsM = new TH1D("bins_m", "bins_m", nBins-1, xbins);
	binsM->Sumw2();
        TH1D * binsM_tev = new TH1D("bins_m_tev", "bins_m_tev", nBins-1, xbins_tev);
        binsM_tev->Sumw2();



	for (Int_t j=0; j<nBins-1; j++) {
		float x=xbins[j+1]-xbins[j];
		binsM->Fill(xbins[j]+0.5*x,x);
                float xx=xbins_tev[j+1]-xbins_tev[j];
                binsM_tev->Fill(xbins_tev[j]+0.5*xx,xx);
                // cout << j << " " << xbins[j+1] << " bin size=" << x << endl;
	}
	// set bin errors to 0 (we use it to divide the bin width!)
          for (int i=0 ; i<(binsM->GetNbinsX()); i++) {
                                binsM->SetBinError(  i+1, 0.0);
                                binsM_tev->SetBinError(  i+1, 0.0);
           }



// 9 invariant masses
    TH1D* Mjj=new TH1D( "Mjj", "Jet Jet Mass", nBins-1, xbins);Mjj->Sumw2();
    TH1D* Mjb=new TH1D( "Mjb", "Jet Bjet Mass", nBins-1, xbins);Mjb->Sumw2();
    TH1D* Mbb=new TH1D( "Mbb", "BB jet  Mass", nBins-1, xbins);Mbb->Sumw2();
    TH1D* Mje=new TH1D( "Mje", "jet+e", nBins-1, xbins);Mje->Sumw2();
    TH1D* Mjm=new TH1D( "Mjm", "jet+m", nBins-1, xbins);Mjm->Sumw2();
    TH1D* Mjg=new TH1D( "Mjg", "jet+g", nBins-1, xbins);Mjg->Sumw2();
    TH1D* Mbe=new TH1D( "Mbe", "b+e mass", nBins-1, xbins);Mbe->Sumw2();
    TH1D* Mbm=new TH1D( "Mbm", "b+m mass", nBins-1, xbins);Mbm->Sumw2();
    TH1D* Mbg=new TH1D( "Mbg", "b+g mass", nBins-1, xbins);Mbg->Sumw2();


	// Begin event loop. Generate event. Skip if error. List first one.
	for (int n = 0; n < Ntot; n++) {
		if (!pythia.next()) continue;
		// if (n < 1) {pythia.info.list(); pythia.event.list();}
		// Loop over particles in event. Find last Z0 copy. Fill its pT.

          
                 // process code
                 m_run=(int)pythia.info.code();
                 m_event=n;

/*
 | g g -> g g                                     111 |   1.499e+04 |
 | g g -> q qbar (uds)                            112 |   7.073e+02 |
 | q g -> q g                                     113 |   6.642e+04 |
 | q q(bar)' -> q q(bar)'                         114 |   3.043e+04 |
 | q qbar -> g g                                  115 |   2.145e+02 |
 | q qbar -> q' qbar' (uds)                       116 |   2.299e+02 |
 | g g -> c cbar                                  121 |   2.357e+02 |
 | q qbar -> c cbar                               122 |   7.662e+01 |
 | g g -> b bbar                                  123 |   2.352e+02 |
 | q qbar -> b bbar                               124 |   7.648e+01 |
*/




		// get weights
		double weight=pythia.info.weight();
		h_weightsFine->Fill(weight);
                h_weightsCorse->Fill(weight);
   

		if ( int(Ntot / int(Ntot/10)) == 0) {
			cout << "Number of events  " << n << " passed"  << endl; };

		vector<PseudoJet> avec;
		h_debug->Fill("Generated",1.);
                vector<int> nrid;
                vector<LParticle> candidates;
                vector<LParticle> bquarks;
                double pxsum=0;
                double pysum=0;
                double pzsum=0;


		for (int i =0; i<pythia.event.size(); i++) {
			int status=pythia.event[i].statusHepMC();
			int pdgid=pythia.event[i].id();
			int type=pdgid;
			double ee=pythia.event[i].e();
			double px=pythia.event[i].px();
			double py=pythia.event[i].py();
			double pz=pythia.event[i].pz();
			// double mm=pythia.event[i].m();
			//double xx=pythia.event[i].xProd();
			//double yy=pythia.event[i].yProd();
			//double zz=pythia.event[i].zProd();
			//double tt=pythia.event[i].tProd();
			double pt=sqrt(px*px+py*py);
			double eta=-log(tan(atan2(pt,(double)pz)/2));
			if ( pt < 0.2)                   continue;
			if ( fabs(eta)> 3.0 )            continue;

			// b-quarks
			if (abs(type) ==5 and pt>0.3*ptJet) {
				int charge=1;
				if (pdgid<0) charge=-1;
				LParticle b(px,py,pz,ee,charge);
				b.SetStatus(0);
				b.SetType(type);
				bquarks.push_back(b);
			};


			// only final states, no neutrinos
			if (status !=1) continue;
			if (abs(pdgid)==12 || abs(pdgid)==14 || abs(pdgid)==16 ) continue;



                                // muon (11) or electron (13) or photon (22) or tau (15)
                                if (abs(type) ==11 || abs(type) ==13 || abs(type) ==22 || abs(type) ==15) {
                                        int charge=1;
                                        if (type<0) charge=-1;
                                        LParticle p(px,py,pz,ee,charge);
                                        p.SetStatus( i );
                                        p.SetType(type);
                                        if (pt>ptLepton && TMath::Abs(eta)<etaLep) {
                                                if (debug) cout << i << " type=" << type << " pt=" << pt << "eta=" << eta << endl;
                                                candidates.push_back(p);
                                        }

                                };
                                // int charge=chargemap[pa->pdg_id(j)]; // get charge
                                if ( pt < 0.1)                   continue;
                                if ( fabs(eta)> 3.5 )            continue;
                                avec.push_back(  PseudoJet(px,py,pz,ee)  );
                                nrid.push_back(i);

                                pxsum=pxsum+px;
                                pysum=pysum+py;
                                pzsum=pzsum+pz;

		} // end particle loop


                 if (debug) cout << "Particles filled" << endl;




                        // isolate candidates from hadronic activity making sure that
                        // energy around this candidates is 90% inside the cone of 0.1
                        vector<LParticle> leptons;
                        double IsoEnergy=0.1;
                        // isolate leading lepton
                        for (unsigned int k1 = 0; k1<candidates.size(); k1++) {
                                LParticle xlep=(LParticle)candidates.at(k1);
                                TLorentzVector lep=xlep.GetP();
                                double p_pt=lep.Et();
                                double p_eta=lep.PseudoRapidity();
                                double p_phi=lep.Phi();
                                if (p_phi<0) p_phi=k2PI+p_phi;
                                if (p_pt<ptLepton) continue;
                                double esumP=0;
                                for (unsigned int k2 = 0; k2<avec.size(); k2++) {
                                        PseudoJet part = avec.at(k2);
                                        double pt=part.perp();
                                        double eta=part.pseudorapidity();
                                        double phi=part.phi();
                                        if (phi<0) phi=k2PI+phi;
                                        double deta    = p_eta - eta;
                                        double dphi    = p_phi - phi;
                                        double adphi=TMath::Abs(dphi);
                                        if (adphi>kPI) adphi=k2PI-adphi;
                                        double ddr = TMath::Sqrt(deta*deta+adphi*adphi);
                                        if (ddr<Rlep_iso) esumP=esumP+pt;
                                }
                                double isoFrac=esumP/p_pt;
                                //cout << "Esum=" << esumP << " p_pt = " << p_pt << " iso=" << isoFrac << endl;
                                h_iso->Fill(isoFrac);
                                if (isoFrac<1.0+IsoEnergy) {
                                        leptons.push_back(xlep);
                                };
                        }


                    if (debug) cout << "Leptons/Photons isolated" << endl;


                       // Nr of b-quarks
                        h_bquarks->Fill(bquarks.size());


                        unsigned int nLeptons=leptons.size();
                        if (nLeptons>1) std::sort(leptons.begin(), leptons.end(), greater<LParticle>() ) ;
                        h_n_lepton->Fill(nLeptons);
                        h_debug->Fill("Iso leptons",1.0);

                        // remove isolated leptons and photons  from vectors used for jets


                        vector<PseudoJet> hadrons;

                        for (unsigned int k = 0; k<avec.size(); k++) {
                                PseudoJet part = avec.at(k);
                                int id=nrid.at(k);
                                int isLep=false;
                                for (unsigned int ll=0; ll<leptons.size(); ll++){
                                        LParticle LL=leptons.at(ll);
                                        int id_lep=LL.GetStatus();
                                        if (id_lep == id) isLep=true;
                                }
                                if (!isLep) hadrons.push_back(part);
                        }

                       if (debug) cout << "Start making jets" << endl;

                       // make jets
                        ClusterSequence clust_seq(hadrons, jet_def);
                        vector<PseudoJet> jets_truth = clust_seq.inclusive_jets(ptJet);
                        vector<PseudoJet> sorted_jets = sorted_by_pt(jets_truth);
                        vector<LParticle> jets;  // light-flavoured
                        vector<LParticle> bjets; // jets with b-quarks

                        for (unsigned int k = 0; k<sorted_jets.size(); k++) {
                                double eta=sorted_jets[k].pseudorapidity();
                                if ( fabs(eta)> etaJet )            continue;
                                double phi=sorted_jets[k].phi();
                                double pt = sorted_jets[k].perp();
                                if (pt<ptJet)                  continue;
                                double e = sorted_jets[k].e();

                                // misidentified b-jets using misRateB rate
                                int FakeB=0;
                                if (100.0*misrate.Rndm()< (misRateB+0.003*pt) ) FakeB=1;

                                // find and label b-quark jets
                                int matchB=0;
                                for (unsigned int k3=0; k3<bquarks.size(); k3++) {
                                        LParticle p=(LParticle)bquarks.at(k3);
                                        TLorentzVector L= p.GetP();
                                        double pt_b = L.Perp();
                                        double eta_b = L.PseudoRapidity();
                                        double phi_b = L.Phi();
                                        if (phi_b<0) phi_b=k2PI+phi_b;
                                        double deta=eta_b-eta;
                                        double dphi=phi_b-phi;
                                        double dR=sqrt(deta*deta + dphi*dphi);
                                        if (dR<Rsize && pt_b>0.5*pt) matchB=1;
                                }


                                // light flavor jets
                                TLorentzVector l;
                                l.SetPtEtaPhiE(pt, eta, phi, e);
                                LParticle p;
                                p.SetP(l);
                                p.SetType(matchB); // label b-jet quarks
                                int nmulti=sorted_jets[k].constituents().size();
                                p.SetCharge(nmulti); // multiplicity
                                if (matchB==0 && FakeB==0) {
                                        h_pt_jet->Fill(pt);
                                        h_eta_jet->Fill(eta);
                                        jets.push_back(p); // light-flavored jets
                                } else if (matchB==1 || FakeB==1) {
                                        h_pt_bjet->Fill(pt);
                                        h_eta_bjet->Fill(eta);
                                        if (FakeB==1) p.SetType(-1);
                                        bjets.push_back(p); // b-jets and fake b-jets 
                                }

                     } // end loop over jets 


                        unsigned int nJets=jets.size();
                        if (nJets>1) std::sort(jets.begin(), jets.end(), greater<LParticle>() ) ;

                        unsigned int nJetsB=bjets.size();
                        if (nJetsB>1) std::sort(bjets.begin(), bjets.end(), greater<LParticle>() ) ;

                        h_debug->Fill("jets.",1.0);


                       // overlap removal for 3 first jets (how many tipes)
                        int nMaxJet=jets.size();
                        if (nMaxJet>maxNumber) nMaxJet=maxNumber;

                        vector<LParticle> leptons_iso;
                        for (unsigned int ll=0; ll<nLeptons; ll++){
                                LParticle LL=leptons.at(ll);
                                TLorentzVector LP=LL.GetP();
                                double phi_lep=LP.Phi();
                                double eta_lep=LP.PseudoRapidity();
                                double y_lep=LP.Rapidity();

                                bool found=false; // isolated from 3 leading jets
                                for (int ii=0; ii<nMaxJet; ii++){
                                        LParticle LPP=jets.at(ii);
                                        TLorentzVector LP=LPP.GetP();
                                        double phi_jet=LP.Phi();
                                        double eta_jet=LP.PseudoRapidity();
                                        double y_jet=LP.Rapidity();
                                        double deta=TMath::Abs(y_jet-y_lep);
                                        double dphi=TMath::Abs(phi_jet - phi_lep);
                                        if (dphi>kPI) dphi=k2PI-dphi;
                                        double dR=TMath::Sqrt(deta*deta+dphi*dphi);
                                        if (dR<Rjet_iso) found=true;
                                }
                                if (found) continue;
                                leptons_iso.push_back(LL);
                        }

                        // now classify isolated candidates
                        vector<LParticle> missing;
                        vector<LParticle> electrons;
                        vector<LParticle> muons;
                        vector<LParticle> photons;
                        vector<LParticle> taus;

                       // for different charges
                        vector<LParticle> electronsPlus;
                        vector<LParticle> muonsPlus;
                        vector<LParticle> tausPlus;
                        vector<LParticle> electronsMinus;
                        vector<LParticle> muonsMinus;
                        vector<LParticle> tausMinus;

                        for (unsigned int k1 = 0; k1<leptons_iso.size(); k1++) {
                                LParticle xlep=(LParticle)leptons_iso.at(k1);
                                int type=xlep.GetType();
                                if (abs(type) ==11) electrons.push_back(xlep);
                                if (abs(type) ==13) muons.push_back(xlep);
                                if (abs(type) ==22) photons.push_back(xlep);
                                if (abs(type) ==15) taus.push_back(xlep);
                                // positive charge
                                if (type ==11) electronsPlus.push_back(xlep);
                                if (type ==13) muonsPlus.push_back(xlep);
                                if (type ==15) tausPlus.push_back(xlep);
                                // negative charge
                                if (type ==-11) electronsMinus.push_back(xlep);
                                if (type ==-13) muonsMinus.push_back(xlep);
                                if (type ==-15) tausMinus.push_back(xlep);
                        }


                        // missing ET
                        TLorentzVector l;
                        double met_pt=sqrt(pxsum*pxsum+pysum*pysum);
                        double met_phi = 0;
                        if (pxsum>0)
                                met_phi=atan2(pysum,pxsum);
                        if (pxsum<0)
                                met_phi=atan2(pysum,pxsum) + kPI;
                        // MET below pT jet is not considered
                        if (met_pt<ptJet){
                                met_pt=0;
                                met_phi=0;
                        }

                        l.SetPtEtaPhiM(met_pt, 0, met_phi, 0);
                        h_met->Fill(met_pt);
                        l.SetPz(0);
                        LParticle p;
                        p.SetP(l);
                        p.SetType(1);
                        missing.push_back(p);

                        m_proj_index1.clear();
                        m_proj_index2.clear();
                        m_proj.clear();
                        h_process->Fill((double)m_id);


                        if (muons.size()>0) {
                             LParticle LPP1=muons.at(0);
                             TLorentzVector L=LPP1.GetP();
                             h_pt_muon->Fill(L.Perp());
                          }


                        if (electrons.size()>0) {
                             LParticle LPP1=electrons.at(0);
                             TLorentzVector L=LPP1.GetP();
                             h_pt_electron->Fill(L.Perp());
                          }

                        if (photons.size()>0) {
                             LParticle LPP1=photons.at(0);
                             TLorentzVector L=LPP1.GetP();
                             h_pt_photon->Fill(L.Perp());
                          }


                        h_jetN->Fill((float)jets.size());
                        h_bjetN->Fill((float)bjets.size());
                        h_photonN->Fill((float)photons.size());
                        h_muonN->Fill((float)muons.size());
                        h_electronN->Fill((float)electrons.size());

                        // sort in increasing pT
                        if (jets.size()>1)  std::sort(jets.begin(), jets.end(), greater<LParticle>() ) ;
                        if (bjets.size()>1)  std::sort(bjets.begin(), bjets.end(), greater<LParticle>() ) ;
                        if (muons.size()>1) std::sort(muons.begin(), muons.end(), greater<LParticle>() ) ;
                        if (electrons.size()>1) std::sort(electrons.begin(), electrons.end(), greater<LParticle>() ) ;
                        if (taus.size()>1) std::sort(taus.begin(), taus.end(), greater<LParticle>() ) ;
                        if (photons.size()>1) std::sort(photons.begin(), photons.end(), greater<LParticle>() ) ;

                        // plus
                        if (muonsPlus.size()>1) std::sort(muonsPlus.begin(), muonsPlus.end(), greater<LParticle>() ) ;
                        if (electronsPlus.size()>1) std::sort(electronsPlus.begin(), electronsPlus.end(), greater<LParticle>() ) ;
                        if (tausPlus.size()>1) std::sort(tausPlus.begin(), tausPlus.end(), greater<LParticle>() ) ;

                        // minus
                        if (muonsMinus.size()>1) std::sort(muonsMinus.begin(), muonsMinus.end(), greater<LParticle>() ) ;
                        if (electronsMinus.size()>1) std::sort(electronsMinus.begin(), electronsMinus.end(), greater<LParticle>() ) ;
                        if (tausMinus.size()>1) std::sort(tausMinus.begin(), tausMinus.end(), greater<LParticle>() ) ;


                        // we accept non-empty events only!
                        if (met_pt==0 && jets.size()==0 && muons.size()==0 && electrons.size() == 0 && photons.size()==0 && bjets.size()==0) continue;


                       // matrix has size:
                        if (debug) {

                                cout << "# Nr jets=" << jets.size()<< " muons=" <<  muons.size() << " ele=" << electrons.size() << " pho=" << photons.size() << endl;

                                cout << " Nr of Photons=" << photons.size() << endl;
                                for (unsigned int i1=0; i1<photons.size(); i1++){
                                        LParticle LPP1=photons.at(i1);
                                        TLorentzVector LP1=LPP1.GetP();
                                        cout << i1<< " pt=" << LP1.Et() << " eta=" << LP1.Eta() << " energy=" << LP1.E() << endl;
                                }

                                cout << "\n Nr of Jets=" << jets.size() << endl;
                                for (unsigned int i1=0; i1<jets.size(); i1++){
                                        LParticle LPP1=jets.at(i1);
                                        TLorentzVector LP1=LPP1.GetP();
                                        cout << i1 << " pt=" << LP1.Et() << " eta=" << LP1.Eta() << " energy=" << LP1.E() << endl;
                                }




                        };



     // fill here some masses for debugging
     // jj
     for (unsigned int i=0; i< jets.size(); i++){
        if ( jets.size()>1){
            LParticle p1= jets.at(0);
            LParticle p2= jets.at(1);
            TLorentzVector PP=p1.GetP()+p2.GetP();
            mass_jj=PP.M();
            Mjj->Fill(mass_jj, weight);
        }
     }

      // j+b
     for (unsigned int i=0; i< jets.size(); i++){
        if ( bjets.size()>0){
            LParticle p1= jets.at(0);
            LParticle p2= bjets.at(0);
            TLorentzVector PP=p1.GetP()+p2.GetP();
            mass_jb=PP.M();
            Mjb->Fill(mass_jb, weight);
        }
     }

     // b+b
     for (unsigned int i=0; i< bjets.size(); i++){
        if ( bjets.size()>1){
            LParticle p1= bjets.at(0);
            LParticle p2= bjets.at(1);
            TLorentzVector PP=p1.GetP()+p2.GetP();
            mass_bb=PP.M(); 
            Mbb->Fill(mass_bb, weight);
        }
     }

     // je
     for (unsigned int i=0; i< jets.size(); i++){
        if ( electrons.size()>0){
            LParticle p1= jets.at(0);
            LParticle p2= electrons.at(0);
            TLorentzVector PP=p1.GetP()+p2.GetP();
            Mje->Fill(PP.M(), weight);
        }
     }

     // jm 
     for (unsigned int i=0; i< jets.size(); i++){
        if ( muons.size()>0){
            LParticle p1= jets.at(0);
            LParticle p2= muons.at(0);
            TLorentzVector PP=p1.GetP()+p2.GetP();
            Mjm->Fill(PP.M(), weight);
        }
     }

     // jg 
     for (unsigned int i=0; i< jets.size(); i++){
        if ( photons.size()>0){
            LParticle p1= jets.at(0);
            LParticle p2= photons.at(0);
            TLorentzVector PP=p1.GetP()+p2.GetP();
            Mjg->Fill(PP.M(), weight);
        }
     }

     // be 
     for (unsigned int i=0; i< bjets.size(); i++){
        if ( electrons.size()>0){
            LParticle p1= bjets.at(0);
            LParticle p2= electrons.at(0);
            TLorentzVector PP=p1.GetP()+p2.GetP();
            Mbe->Fill(PP.M(), weight);
        }
     }

     // bm 
     for (unsigned int i=0; i< bjets.size(); i++){
        if ( muons.size()>0){
            LParticle p1= bjets.at(0);
            LParticle p2= muons.at(0);
            TLorentzVector PP=p1.GetP()+p2.GetP();
            Mbm->Fill(PP.M(), weight);
        }
     }

     // bg 
     for (unsigned int i=0; i< bjets.size(); i++){
        if ( photons.size()>0){
            LParticle p1= bjets.at(0);
            LParticle p2= photons.at(0);
            TLorentzVector PP=p1.GetP()+p2.GetP();
            Mbg->Fill(PP.M(), weight);
        }
     }



                      // Build RMM 
                      float** projArray =  projectevent(CMenergy, maxNumber, maxTypes, missing, jets, bjets, muons, electrons, photons);



                        // triangular matrix is a special kind of square matrix.
                        // upper triangular matrix or right triangular matrix.
                        // non-zero in triangular matrix
                        if (debug) {
                                cout << "  " << names_debug[0] << "          ";
                                for (int h = 1; h < maxTypes+1; h++)
                                        for (int i = 0; i <  maxNumber; i++)
                                                cout << names_debug[h] << i << "       ";
                                cout << endl;
                                for (int h = 0; h < mSize; h++) {
                                        cout << h << " ";
                                        for (int w = 0; w < mSize; w++)
                                        {
                                                if (h != w) printf("%.2e ", float(projArray[h][w]));
                                                else  {
                                                        //cout << std::setprecision(1) << std::scientific <<  float(projArray[h][w]);
                                                        //printf("%.2e ", float(projArray[h][w]));
                                                        cout << RED; printf("%.2e ", float(projArray[h][w])); cout << RESET;
                                                        // else printf("\033[1;31%.2e033[0m", float(projArray[h][w]));
                                                }

                                        }
                                        printf("\n");
                                }
                        }; // end debug



   int non_empty=0;
   for (int w1 = 0; w1 < mSize; w1++) {
          for (int h1 = 0; h1 < mSize; h1++) {
                  float dd=projArray[w1][h1];
                  if (h1<w1)   dd=dd*angNorm; // decrease angles by 0.15 
                  int i1=h1;
                  int i2=mSize-w1-1;
                  h_proj->Fill((Names2.at(i1)).c_str(),  (Names1.at(i2)).c_str(),dd);
                  h_prof->Fill((Names2.at(i1)).c_str(),  (Names1.at(i2)).c_str(),dd);
                  if (dd>0) {
                   m_proj.push_back(dd);
                   m_proj_index1.push_back( w1 );
                   m_proj_index2.push_back( h1 );
                   non_empty++;
                  }

               }}


            m_tree->Fill();

            for (int w1 = 0; w1 < mSize; w1++)  delete[] projArray[w1];
            delete[] projArray;


// ------------------------------------------------------------------
         } // endl loop over events


	// To check which changes have actually taken effect
	pythia.settings.listChanged();
	// pythia.particleData.listChanged();
	pythia.particleData.list(25);
	// ParticleDataTable::listAll()
	// ParticleDataTable::list(25);


	pythia.stat();


	// Output histograms
	double sigmapb = pythia.info.sigmaGen() * 1.0E9;
	double sigmapb_err = pythia.info.sigmaErr() * 1.0E9;

	cout << "== Run statistics: " << endl;
        cout << "== Seed number   =" <<  Iseed << endl;
	cout << "== Cross section    =" <<  sigmapb << " +- " << sigmapb_err << " pb" << endl;
	cout << "== Generated Events =" <<  Ntot << endl;
	double lumi=(Ntot/sigmapb);
	cout << "== Luminosity       =" <<  lumi  << " pb-1" << endl;
	cout << "\n\n";

	h_cross->Fill("Nr jobs",1.0); // calibration check
	h_cross->Fill("Cross section [pb]",sigmapb);
	h_cross->Fill("Events",Ntot);
	h_cross->Fill("Lumi [pb]",lumi);
        h_cross->Fill("weightSum",pythia.info.weightSum());
        h_cross->Fill("mergingWeight",pythia.info.mergingWeight());
	h_cross->Fill("ptlepton",ptLepton);
	h_cross->Fill("ptjet",ptJet);

   
        h_info->Fill("PythiaVErsion",versionNumber); // calibration check
        h_info->Fill("CM",pythia.info.eCM()); // calibration check
        h_info->Fill("Random Seed",Iseed);



	// also make a graph
        /*
	TGraph *info = new TGraph();
	info->SetPoint(1, 1.0, 1.0); // calibration
	info->SetPoint(2, 2.0,versionNumber);
	info->SetPoint(3, 3.0,sigmapb);
	info->SetPoint(4, 4.0,Ntot);
	info->SetPoint(5, 5.0,lumi);
	info->SetPoint(6, 6.0,ptLepton);
	info->SetPoint(7, 7.0,ptJet);
	info->SetPoint(8, 8.0,etaJet);
	info->SetPoint(9, 9.0,R);
	info->SetPoint(10, 10.0,Rlepton);
	info->Write("info"); // save the graph
        */


	RootFile->Write();
	RootFile->Print();
	RootFile->Close();


	return 0;
}
