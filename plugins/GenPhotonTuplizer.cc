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
	genInfoToken_		= consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"));
	lheInfoToken_		= consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheInfo"));

	edm::Service<TFileService> fs;
	GenEvents = fs->make<TTree>("GenEvents", "gen events");	
	sumWeights = fs->make<TH1D>("sumWeights",";;sumWeights;",1,-0.5,0.5);
	sumWeights_Hgg = fs->make<TH1D>("sumWeights_Hgg",";;sumWeights_Hgg;",1,-0.5,0.5);//sum of weights of events with Hgg
	
	sumScaleWeights = fs->make<TH1D>("sumScaleWeights",";;sumScaleWeights;",9,-0.5,8.5);
	sumPdfWeights = fs->make<TH1D>("sumPdfWeights",";;sumPdfWeights;",100,-0.5,99.5);

	sumScaleWeights_Hgg = fs->make<TH1D>("sumScaleWeights_Hgg",";;sumScaleWeights_Hgg;",9,-0.5,8.5);
	sumPdfWeights_Hgg = fs->make<TH1D>("sumPdfWeights_Hgg",";;sumPdfWeights_Hgg;",100,-0.5,99.5);

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
	genWeight = genInfo->weight();
	sumWeights->Fill(0.,genWeight);

	double nomlheweight = lheInfo->weights()[0].wgt;
	
	if (lheInfo->weights().size()>=9) 
	{  
	  	for (unsigned int iwgt=0; iwgt<9; ++iwgt) 
		{
	    		double wgtval = lheInfo->weights()[iwgt].wgt*genWeight/nomlheweight;
	    		scaleWeights->push_back(wgtval);
			sumScaleWeights->Fill(double(iwgt),(*scaleWeights)[iwgt]);
	  	}
	}
	
	if (lheInfo->weights().size()>=108) 
	{
		
	  	for (unsigned int iwgt=9; iwgt<109; ++iwgt) 
		{
			double wgtval = lheInfo->weights()[iwgt].wgt/nomlheweight;
	    		pdfWeights->push_back(wgtval);
			sumPdfWeights->Fill(double(iwgt-9),(*pdfWeights)[iwgt-9]);
		}
	}

//genJets
	fillGenJets();
	
//genParticles
	fillGenParticles();

//fill Hgg events - GEN level
	fillGenHgg();	

//fill output ntuple
	GenEvents->Fill();
}


//load all collection from cfg file
void GenPhotonTuplizer::loadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	iEvent.getByToken(genParticlesToken_,genParticles);
	iEvent.getByToken(genJetsToken_,genJets);
	iEvent.getByToken(genInfoToken_,genInfo);
	iEvent.getByToken(lheInfoToken_, lheInfo);

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
	genJetIsFromHiggsPhoton = new std::vector<bool>;
	genJetNumberOfDaughters = new std::vector<int>;

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
	
	scaleWeights = new std::vector<float>; scaleWeights->clear();
	pdfWeights = new std::vector<float>; pdfWeights->clear();

	genJetE->clear();
	genJetPt->clear();
	genJetEta->clear();
	genJetPhi->clear();
	genJetIsFromHiggsPhoton->clear();
	genJetNumberOfDaughters->clear();
	
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
	GenEvents->Branch("genWeight", &genWeight, "genWeight/F");
	GenEvents->Branch("scaleWeights", "vector<float>",&scaleWeights);
	GenEvents->Branch("pdfWeights", "vector<float>",&pdfWeights);
  	GenEvents->Branch("nGenParticles", &nGenParticles, "nGenParticles/I");
  	GenEvents->Branch("nGenJets", &nGenJets, "nGenJets/I");
  	GenEvents->Branch("nGenJets_cut", &nGenJets_cut, "nGenJets_cut/I");
  	GenEvents->Branch("nGenJets_cut_noND", &nGenJets_cut_noND, "nGenJets_cut_noND/I");
	GenEvents->Branch("genJetE", "vector<float>", &genJetE);
	GenEvents->Branch("genJetPt", "vector<float>", &genJetPt);
	GenEvents->Branch("genJetEta", "vector<float>", &genJetEta);
	GenEvents->Branch("genJetPhi", "vector<float>", &genJetPhi);
	GenEvents->Branch("genJetIsFromHiggsPhoton", "vector<bool>", &genJetIsFromHiggsPhoton);
	GenEvents->Branch("genJetNumberOfDaughters", "vector<int>", &genJetNumberOfDaughters);
  	
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
  
		
	GenEvents->Branch("genHiggs_pho1_E", &genHiggs_pho1_E, "genHiggs_pho1_E/F");
	GenEvents->Branch("genHiggs_pho1_Pt", &genHiggs_pho1_Pt, "genHiggs_pho1_Pt/F");
	GenEvents->Branch("genHiggs_pho1_Eta", &genHiggs_pho1_Eta, "genHiggs_pho1_Eta/F");
	GenEvents->Branch("genHiggs_pho1_Phi", &genHiggs_pho1_Phi, "genHiggs_pho1_Phi/F");
	GenEvents->Branch("genHiggs_pho1_Iso", &genHiggs_pho1_Iso, "genHiggs_pho1_Iso/F");

	GenEvents->Branch("genHiggs_pho2_E", &genHiggs_pho2_E, "genHiggs_pho2_E/F");
	GenEvents->Branch("genHiggs_pho2_Pt", &genHiggs_pho2_Pt, "genHiggs_pho2_Pt/F");
	GenEvents->Branch("genHiggs_pho2_Eta", &genHiggs_pho2_Eta, "genHiggs_pho2_Eta/F");
	GenEvents->Branch("genHiggs_pho2_Phi", &genHiggs_pho2_Phi, "genHiggs_pho2_Phi/F");
	GenEvents->Branch("genHiggs_pho2_Iso", &genHiggs_pho2_Iso, "genHiggs_pho2_Iso/F");

	GenEvents->Branch("genHiggs_M", &genHiggs_M, "genHiggs_M/F");
	GenEvents->Branch("genHiggs_Pt", &genHiggs_Pt, "genHiggs_Pt/F");
	GenEvents->Branch("genHiggs_Eta", &genHiggs_Eta, "genHiggs_Eta/F");
	GenEvents->Branch("genHiggs_Phi", &genHiggs_Phi, "genHiggs_Phi/F");

}

// clear the content of all variables in the output ntuple
void GenPhotonTuplizer::resetBranches()
{
	runNum = -1;
	lumiNum = -1;
	eventNum = -1;

	genWeight = 1;	
	scaleWeights->clear();
	pdfWeights->clear();

	nGenJets = 0;
	nGenJets_cut = 0;
	nGenJets_cut_noND = 0;
	genJetE->clear();
	genJetPt->clear();
	genJetEta->clear();
	genJetPhi->clear();
	genJetIsFromHiggsPhoton->clear();
	genJetNumberOfDaughters->clear();
	

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

	genHiggs_pho1_Pt = 0.0;
	genHiggs_pho1_Eta = 0.0;
	genHiggs_pho1_Phi = 0.0;
	genHiggs_pho1_E = 0.0;
	genHiggs_pho1_Iso = 0.0;
	
	genHiggs_pho2_Pt = 0.0;
	genHiggs_pho2_Eta = 0.0;
	genHiggs_pho2_Phi = 0.0;
	genHiggs_pho2_E = 0.0;
	genHiggs_pho2_Iso = 0.0;

	genHiggs_M = 0.0;
	genHiggs_Pt = 0.0;
	genHiggs_Eta = 0.0;
	genHiggs_Phi = 0.0;
	
}


void GenPhotonTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      	edm::ParameterSetDescription desc;
      	desc.setUnknown();
      	descriptions.addDefault(desc);
}

void GenPhotonTuplizer::fillGenHgg()
{
	std::map<float, int> index_pt_map;
	
	for(int i=0;i<nGenParticles;i++)
	{
		if(genParticleId->at(i) == 22 && genParticleMotherId->at(i) == 25 && genParticleStatus->at(i)==1) 
		{
			index_pt_map.insert(std::pair<float, int>(-1.0*genParticlePt->at(i), i));
		}
	}

	if(index_pt_map.size() < 2) return;
	
	int ind_pho1 = 0;
	int ind_pho2 = 0;
	int ind_temp = 0;
	for (auto tmp :index_pt_map)
	{
		if(ind_temp == 0) ind_pho1 = tmp.second;
		else if(ind_temp == 1) ind_pho2 = tmp.second;
		else break;
		ind_temp ++;
	}

	TLorentzVector pho1(0,0,0,0);
	TLorentzVector pho2(0,0,0,0);
	TLorentzVector higgs(0,0,0,0);
	
	pho1.SetPtEtaPhiE(genParticlePt->at(ind_pho1), genParticleEta->at(ind_pho1),  genParticlePhi->at(ind_pho1),  genParticleE->at(ind_pho1));
	pho2.SetPtEtaPhiE(genParticlePt->at(ind_pho2), genParticleEta->at(ind_pho2),  genParticlePhi->at(ind_pho2),  genParticleE->at(ind_pho2));
	
	higgs = pho1 + pho2;
	
	genHiggs_pho1_Pt = pho1.Pt();	
	genHiggs_pho1_Eta = pho1.Eta();	
	genHiggs_pho1_Phi = pho1.Phi();	
	genHiggs_pho1_E = pho1.E();	
	
	genHiggs_pho2_Pt = pho2.Pt();	
	genHiggs_pho2_Eta = pho2.Eta();	
	genHiggs_pho2_Phi = pho2.Phi();	
	genHiggs_pho2_E = pho2.E();	

	genHiggs_M = higgs.M();
	genHiggs_Pt = higgs.Pt();
	genHiggs_Eta = higgs.Eta();
	genHiggs_Phi = higgs.Phi();
	
	if(genHiggs_M>10.0)
	{
		sumWeights_Hgg->Fill(0.,genWeight);	
		for (unsigned int iwgt=0; iwgt<scaleWeights->size(); ++iwgt) 
		{
      			sumScaleWeights_Hgg->Fill(double(iwgt),(*scaleWeights)[iwgt]);
		}
	
		for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) 
		{
      			sumPdfWeights_Hgg->Fill(double(iwgt),(*pdfWeights)[iwgt]);
		}
	}
//gen level isolation
	float iso_sum_pho1 = 0.0;
	float iso_sum_pho2 = 0.0;
	
	for(size_t i=0; i<genParticles->size();i++)
        {
		if((*genParticles)[i].status()==1)
		{
			if(reco::deltaR((*genParticles)[i].eta(), (*genParticles)[i].phi(), pho1.Eta(), pho1.Phi()) < isoConeSize && fabs(pho1.Pt() - (*genParticles)[i].pt()) > 0.001 ) iso_sum_pho1 += (*genParticles)[i].et();			
			if(reco::deltaR((*genParticles)[i].eta(), (*genParticles)[i].phi(), pho2.Eta(), pho2.Phi()) < isoConeSize && fabs(pho2.Pt() - (*genParticles)[i].pt()) > 0.001 ) iso_sum_pho2 += (*genParticles)[i].et();			
		}
	}
	
	genHiggs_pho1_Iso = iso_sum_pho1;
	genHiggs_pho2_Iso = iso_sum_pho2;

}

void GenPhotonTuplizer::fillGenJets()
{

	for(const reco::GenJet &j : *genJets)
	{
		genJetE->push_back(j.energy());
		genJetPt->push_back(j.pt());
		genJetEta->push_back(j.eta());
		genJetPhi->push_back(j.phi());
		nGenJets++;
		bool fromGenPhoton = false;
		int numberOfDaughters = j.numberOfDaughters();
		genJetNumberOfDaughters->push_back(numberOfDaughters);

		std::vector< const reco::GenParticle * > jetGenParticles = j.getGenConstituents();
		for(size_t i=0; i<jetGenParticles.size();i++)
		{
			const reco::Candidate * prunedV = 0;
			prunedV =  jetGenParticles[i];
			if(prunedV->pdgId() == 22)
			{

				if(prunedV->numberOfMothers() > 0)
				{
					const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV);
					if (firstMotherWithDifferentID && firstMotherWithDifferentID->pdgId() == 25)
					{
						fromGenPhoton = true;	
						break;
					}
				}

			}
		}
	
		genJetIsFromHiggsPhoton->push_back(fromGenPhoton);

		if((!fromGenPhoton) && fabs(j.eta())<2.5 && j.pt()>30 && numberOfDaughters>5)
		{	
			nGenJets_cut ++;
		}
		if((!fromGenPhoton) && fabs(j.eta())<2.5 && j.pt()>30)
		{	
			nGenJets_cut_noND ++;
		}

    	}

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

