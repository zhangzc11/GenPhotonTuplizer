// -*- C++ -*-
/*
 *   Description: reconstruction of pi0/eta and save the ntuple
*/
// Author: Zhicai Zhang
// Created: Fri Mar 24 11:01:09 CET 2017


#include "GenPhotonTuplizer.h"

using namespace std;

#define DEBUG


GenPhotonTuplizer::GenPhotonTuplizer(const edm::ParameterSet& iConfig)
{
	//get parameters from iConfig
 	genParticlesToken_ 	= consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
 	genJetsToken_ 		= consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets")); 

	edm::Service<TFileService> fs;
	GenEvents = fs->make<TTree>("GenEvents", "gen events");	

}


GenPhotonTuplizer::~GenPhotonTuplizer()
{

}


//------ Method called for each event ------//
void GenPhotonTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


        resetBranches();
        loadEvent(iEvent, iSetup);

//fill event info
        runNum = iEvent.id().run();
        lumiNum = iEvent.luminosityBlock();
        eventNum = iEvent.id().event();

//genJets
	for(const reco::GenJet &j : *genJets){
        genJetE->push_back(j.energy());
        genJetPt->push_back(j.pt());
        genJetEta->push_back(j.eta());
        genJetPhi->push_back(j.phi());
        nGenJets++;
    	}
//genParticles
	fillGenParticles();

//fill output ntuple
	GenEvents->Fill();
}


//load all collection from cfg file
void GenPhotonTuplizer::loadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	iEvent.getByToken(genParticlesToken_,genParticles);
	iEvent.getByToken(genJetsToken_,genJets);

}



// ------------ method called once each job just before starting event loop  ------------
void GenPhotonTuplizer::beginJob()
{


	setBranches();


}

//------ Method called once each job just after ending the event loop ------//
void GenPhotonTuplizer::endJob()
{

}

// ------------ method called when starting to processes a run  ------------
void GenPhotonTuplizer::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{

}

// ------------ method called when ending the processing of a run  ------------
void GenPhotonTuplizer::endRun(edm::Run const&, edm::EventSetup const&)
{

}


// set branch address
void GenPhotonTuplizer::setBranches()
{



	genJetE = new std::vector<float>;
	genJetPt = new std::vector<float>;
	genJetEta = new std::vector<float>;
	genJetPhi = new std::vector<float>;

	genParticleMotherId = new std::vector<int>;
	genParticleMotherIndex = new std::vector<int>;
	genParticleId = new std::vector<int>;
	genParticleStatus = new std::vector<int>;
	genParticleIsPromptFinalState = new std::vector<bool>;
	genParticleIsPromptDecayed = new std::vector<bool>;
	genParticleE = new std::vector<float>;
	genParticlePt = new std::vector<float>;
	genParticleEta = new std::vector<float>;
	genParticlePhi = new std::vector<float>;

	genJetE->clear();
	genJetPt->clear();
	genJetEta->clear();
	genJetPhi->clear();
	
	genParticleMotherId->clear();
	genParticleMotherIndex->clear();
	genParticleId->clear();
	genParticleStatus->clear();
	genParticleIsPromptFinalState->clear();
	genParticleIsPromptDecayed->clear();
	genParticleE->clear();
	genParticlePt->clear();
	genParticleEta->clear();
	genParticlePhi->clear();
	
        //diphoton ntuple
  	GenEvents->Branch("runNum", &runNum, "runNum/i");
  	GenEvents->Branch("lumiNum", &lumiNum, "lumiNum/i");
  	GenEvents->Branch("eventNum", &eventNum, "eventNum/i");
  	GenEvents->Branch("nGenParticles", &nGenParticles, "nGenParticles/I");
  	GenEvents->Branch("nGenJets", &nGenJets, "nGenJets/I");
	GenEvents->Branch("genJetE", "vector<float>", &genJetE);
	GenEvents->Branch("genJetPt", "vector<float>", &genJetPt);
	GenEvents->Branch("genJetEta", "vector<float>", &genJetEta);
	GenEvents->Branch("genJetPhi", "vector<float>", &genJetPhi);
  	
	GenEvents->Branch("nGenParticles", &nGenParticles, "nGenParticles/I");
	GenEvents->Branch("genParticleMotherId", "vector<int>", &genParticleMotherId);
	GenEvents->Branch("genParticleMotherIndex", "vector<int>", &genParticleMotherIndex);
	GenEvents->Branch("genParticleId", "vector<int>", &genParticleId);
	GenEvents->Branch("genParticleStatus", "vector<int>", &genParticleStatus);
	GenEvents->Branch("genParticleIsPromptFinalState", "vector<bool>", &genParticleIsPromptFinalState);
	GenEvents->Branch("genParticleIsPromptDecayed", "vector<bool>", &genParticleIsPromptDecayed);
	GenEvents->Branch("genParticleE", "vector<float>", &genParticleE);
	GenEvents->Branch("genParticlePt", "vector<float>", &genParticlePt);
	GenEvents->Branch("genParticleEta", "vector<float>", &genParticleEta);
	GenEvents->Branch("genParticlePhi", "vector<float>", &genParticlePhi);
  	

}

// clear the content of all variables in the output ntuple
void GenPhotonTuplizer::resetBranches()
{
	runNum = -1;
	lumiNum = -1;
	eventNum = -1;
	

	nGenJets = 0;
	genJetE->clear();
	genJetPt->clear();
	genJetEta->clear();
	genJetPhi->clear();
	

	nGenParticles = 0;
	genParticleMotherId->clear();
	genParticleMotherIndex->clear();
	genParticleId->clear();
	genParticleStatus->clear();
	genParticleIsPromptFinalState->clear();
	genParticleIsPromptDecayed->clear();
	genParticleE->clear();
	genParticlePt->clear();
	genParticleEta->clear();
	genParticlePhi->clear();
}


void GenPhotonTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      	edm::ParameterSetDescription desc;
      	desc.setUnknown();
      	descriptions.addDefault(desc);
}


void GenPhotonTuplizer::fillGenParticles()
{
	std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding
	std::vector<bool> v_isPromptFinalState;
  	std::vector<bool> v_isPromptDecayed;
	
	for(size_t i=0; i<genParticles->size();i++)
	{
    		if(
       			(abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6
        		&& ( (*genParticles)[i].status() < 30)) 
       			|| (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
       			|| (abs((*genParticles)[i].pdgId()) == 21
           		&& (*genParticles)[i].status() < 30)
       			|| (abs((*genParticles)[i].pdgId()) >= 22 && abs((*genParticles)[i].pdgId()) <= 25
           		&& ( (*genParticles)[i].status() < 30))
       			|| (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
       			|| (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039)
       		){
      		prunedV.push_back(&(*genParticles)[i]);
        	//v_isPromptFinalState.push_back((*genParticles)[i].isPromptFinalState());
        	//v_isPromptDecayed.push_back((*genParticles)[i].isPromptDecayed());
    		}
  	}

	nGenParticles = prunedV.size();
	
	for(unsigned int i = 0; i < prunedV.size(); i++)
	{
		genParticleId->push_back(prunedV[i]->pdgId());
    		genParticleStatus->push_back(prunedV[i]->status());
    	//	genParticleIsPromptDecayed->push_back(v_isPromptDecayed[i]);
    	//	genParticleIsPromptFinalState->push_back(v_isPromptFinalState[i]);
    		genParticleE->push_back(prunedV[i]->energy());
    		genParticlePt->push_back(prunedV[i]->pt());
		genParticleEta->push_back(prunedV[i]->eta());
		genParticlePhi->push_back(prunedV[i]->phi());

		int temp_MotherId = 0;
		int temp_MotherIndex = -1;

		if(prunedV[i]->numberOfMothers() > 0)
		{
			//find the ID of the first mother that has a different ID than the particle itself
			const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV[i]);
			if (firstMotherWithDifferentID) 
			{
				temp_MotherId = firstMotherWithDifferentID->pdgId();
			}
			//find the mother and keep going up the mother chain if the ID's are the same
			const reco::Candidate* originalMotherWithSameID = findOriginalMotherWithSameID(prunedV[i]);
			for(unsigned int j = 0; j < prunedV.size(); j++)
			{
				if(prunedV[j] == originalMotherWithSameID)
				{
					temp_MotherIndex = j;
					break;
				}
			}
		}
		else 
		{
      			temp_MotherIndex = -1;
    		}
		genParticleMotherId->push_back(temp_MotherId);
		genParticleMotherIndex->push_back(temp_MotherIndex);
	}

}

const reco::Candidate* GenPhotonTuplizer::findFirstMotherWithDifferentID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is this the first parent with a different ID? If yes, return, otherwise
  // go deeper into recursion
   if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->pdgId() == particle->mother(0)->pdgId()
        && particle->mother(0)->status() != 11  // prevent infinite loop for sherpa documentation gluons
        ) {
      return findFirstMotherWithDifferentID(particle->mother(0));
    } else {
      return particle->mother(0);
    }
  }

  return 0;
}

const reco::Candidate* GenPhotonTuplizer::findOriginalMotherWithSameID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }
// Is there another parent with the same ID? If yes, go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->mother(0)->numberOfMothers() == 0 ||
        particle->mother(0)->status() == 11 ||  // prevent infinite loop for sherpa documentation gluons
        (particle->mother(0)->numberOfMothers() > 0 && particle->mother(0)->mother(0)->pdgId() != particle->mother(0)->pdgId())
        ) {
      return particle->mother(0);
    } else {
      return findOriginalMotherWithSameID(particle->mother(0));
    }
  }

  return 0;
}


//define this as a plug-in
DEFINE_FWK_MODULE(GenPhotonTuplizer);

