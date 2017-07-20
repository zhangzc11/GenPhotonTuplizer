
#ifndef GENPHOTONTUPLIZER_H
#define GENPHOTONTUPLIZER_H


/***********************c++*********************/
#include <memory>
#include <string>
#include <vector>
#include <tuple>

/***********************root********************/
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TTree.h"

/***********************cmssw*******************/
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"



#include "Geometry/Records/interface/CaloGeometryRecord.h"



using namespace std;
using namespace edm;

#define OBJECTARRAYSIZE 1000
#define GENPARTICLEARRAYSIZE 500

class GenPhotonTuplizer : public edm::EDAnalyzer {
public:
 //analyzer constructor and destructor
  	explicit GenPhotonTuplizer(const edm::ParameterSet&);
  	 ~GenPhotonTuplizer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   
private:
      	virtual void beginJob() ;
      	virtual void analyze(const edm::Event&, const edm::EventSetup&);
      	virtual void endJob() ;

      	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      	virtual void endRun(edm::Run const&, edm::EventSetup const&);
   
 //specific functions 
   	void loadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup); //call at the beginning of each event to get input handles from the python config
	virtual void resetBranches(); // clear all variables
   	virtual void setBranches(); // set branch of ntuple
	virtual void fillGenParticles();
	const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle);
	const reco::Candidate* findOriginalMotherWithSameID(const reco::Candidate *particle);

 //output TTree and file
      	TTree *GenEvents;

 //variables to be saved in the pi0 ntuple
      	uint    runNum;
      	uint    lumiNum;
      	uint    eventNum;

	int 	nGenParticles;
	vector<int> * genParticleMotherId;
	vector<int> * genParticleMotherIndex;
	vector<int> * genParticleId;
	vector<bool> * genParticleIsPromptDecayed;
	vector<bool> * genParticleIsPromptFinalState;
	vector<int> * genParticleStatus;
	vector<float> * genParticleE;
	vector<float> * genParticlePt;
	vector<float> * genParticleEta;
	vector<float> * genParticlePhi;

	int nGenJets;
	vector<float> * genJetE;
	vector<float> * genJetPt;
	vector<float> * genJetEta;
	vector<float> * genJetPhi;
	

	
 //input tags
 edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
 edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
 //input collections
 edm::Handle<reco::GenParticleCollection> genParticles;
 edm::Handle<reco::GenJetCollection> genJets;

 //cuts and options read from cfg file	
 
 
};

#endif
