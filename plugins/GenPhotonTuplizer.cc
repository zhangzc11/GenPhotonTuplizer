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
 	genParticlesToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));

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

	nGenParticles = genParticles->size();

//fill output ntuple
	GenEvents->Fill();
}


//load all collection from cfg file
void GenPhotonTuplizer::loadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	iEvent.getByToken(genParticlesToken_,genParticles);
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
	
        //diphoton ntuple
  	GenEvents->Branch("runNum", &runNum, "runNum/i");
  	GenEvents->Branch("lumiNum", &lumiNum, "lumiNum/i");
  	GenEvents->Branch("eventNum", &eventNum, "eventNum/i");
  	GenEvents->Branch("nGenParticles", &nGenParticles, "nGenParticles/I");
  	
}

// clear the content of all variables in the output ntuple
void GenPhotonTuplizer::resetBranches()
{
	runNum = -1;
	lumiNum = -1;
	eventNum = -1;
	nGenParticles = 0;
}


void GenPhotonTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.setUnknown();
      descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenPhotonTuplizer);

