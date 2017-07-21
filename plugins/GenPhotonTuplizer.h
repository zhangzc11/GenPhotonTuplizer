
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

#include "DataFormats/Math/interface/deltaR.h"

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
	virtual void fillGenJets();
	virtual void fillGenHgg();
	const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle);
	const reco::Candidate* findOriginalMotherWithSameID(const reco::Candidate *particle);

 //output TTree and histograms
      	TTree *GenEvents;

	TH1D *sumWeights;
	TH1D *sumWeights_Hgg;

 //variables to be saved in the pi0 ntuple
      	uint    runNum;
      	uint    lumiNum;
      	uint    eventNum;

	float genWeight;

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
	int nGenJets_cut;
	vector<float> * genJetE;
	vector<float> * genJetPt;
	vector<float> * genJetEta;
	vector<float> * genJetPhi;
	vector<bool> * genJetIsFromHiggsPhoton;


	float genHiggs_pho1_E;		
	float genHiggs_pho1_Pt;		
	float genHiggs_pho1_Eta;		
	float genHiggs_pho1_Phi;		
	float genHiggs_pho1_Iso;		

	float genHiggs_pho2_E;		
	float genHiggs_pho2_Pt;		
	float genHiggs_pho2_Eta;		
	float genHiggs_pho2_Phi;		
	float genHiggs_pho2_Iso;		

	float genHiggs_M;
	float genHiggs_Pt;
	float genHiggs_Eta;
	float genHiggs_Phi;

	
 //input tags
 	edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
 	edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
 	edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
	
 //input collections
 	edm::Handle<reco::GenParticleCollection> genParticles;
 	edm::Handle<reco::GenJetCollection> genJets;
	edm::Handle<GenEventInfoProduct> genInfo;

 //cuts and options
 	const float isoConeSize = 0.3;
	 
 
};

#endif
