/*
Simple macro showing how to access branches from the delphes output root file,
loop over events.
root -l examples/BranchAnalysis.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

void BranchAnalysis(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  // Book histograms
  TList *l = new TList();
  TH1 *PartID = new TH1F("PID", "PID", 100, 1000, 4000);

  // Chi Histograms
  TH1 *chiStatus= new TH1F("chi_Status", "chi_Status", 100, 0, 10);
  TH1 *chiPT = new TH1F("chi_PT", "chi_PT", 100, 0, 1000);
  TH1 *chiE = new TH1F("chi_E", "chi_E", 100, 0, 2000);
  TH1 *chiEta = new TH1F("chi_Eta", "chi_Eta", 100, -7, 7);
  TH1 *chiPhi = new TH1F("chi_phi", "chi_phi", 100, -3.5, 3.5);

  // Scalar Histograms
  TH1 *scalarStatus= new TH1F("scalar_Status", "scalar_Status", 100, 0, 10);
  TH1 *scalarPT = new TH1F("scalar_PT", "scalar_PT", 100, 0, 1000);
  TH1 *scalarE = new TH1F("scalar_E", "scalar_E", 100, 0, 2000);
  TH1 *scalarEta = new TH1F("scalar_Eta", "scalar_Eta", 100, -7, 7);
  TH1 *scalarPhi = new TH1F("scalar_phi", "scalar_phi", 100, -3.5, 3.5);

  // Psi Histograms
  TH1 *psiStatus= new TH1F("psi_Status", "psi_Status", 100, 0, 10);
  TH1 *psiPT = new TH1F("psi_PT", "psi_PT", 100, 0, 1000);
  TH1 *psiE = new TH1F("psi_E", "psi_E", 100, 0, 2000);
  TH1 *psiEta = new TH1F("psi_Eta", "psi_Eta", 100, -7, 7);
  TH1 *psiPhi = new TH1F("psi_phi", "psi_phi", 100, -3.5, 3.5);

  // Jet Histograms
  TH1 *JetPT = new TH1F("Jet_PT", "Jet_PT", 100, 0, 1000);
  TH1 *JetEta = new TH1F("Jet_Eta", "Jet_Eta", 100, -7, 7);
  TH1 *JetPhi = new TH1F("Jet_phi", "Jet_phi", 100, -3.5, 3.5);

  // Electron Histograms
  TH1 *elecPT = new TH1F("elec_PT", "elec_PT", 100, 0, 1000);
  TH1 *elecEta = new TH1F("elec_Eta", "elec_Eta", 100, -7, 7);
  TH1 *elecPhi = new TH1F("elec_phi", "elec_phi", 100, -3.5, 3.5);

  // Muon Histograms
  TH1 *MuonStatus= new TH1F("Muon_Status", "Muon_Status", 100, 0, 10);
  TH1 *MuonE = new TH1F("Muon_E", "Muon_E", 100, 0, 2000);
  TH1 *MuonPT = new TH1F("Muon_PT", "Muon_PT", 100, 0, 1000);
  TH1 *MuonEta = new TH1F("Muon_Eta", "Muon_Eta", 100, -7, 7);
  TH1 *MuonPhi = new TH1F("Muon_phi", "Muon_phi", 100, -3.5, 3.5);

  //NN1 Histograms
  TH1 *NN1Status= new TH1F("NN1_Status", "NN1_Status", 100, 0, 10);
  TH1 *NN1E = new TH1F("NN1_E", "NN1_E", 100, 0, 2000);
  TH1 *NN1PT = new TH1F("NN1_PT", "NN1_PT", 100, 0, 1000);
  TH1 *NN1Eta = new TH1F("NN1_Eta", "NN1_Eta", 100, -7, 7);
  TH1 *NN1Phi = new TH1F("NN1_phi", "NN1_phi", 100, -3.5, 3.5);

  //NN2 Histograms
  TH1 *NN2Status= new TH1F("NN2_Status", "NN2_Status", 100, 0, 10);
  TH1 *NN2E = new TH1F("NN2_E", "NN2_E", 100, 0, 2000);
  TH1 *NN2PT = new TH1F("NN2_PT", "NN2_PT", 100, 0, 1000);
  TH1 *NN2Eta = new TH1F("NN2_Eta", "NN2_Eta", 100, -7, 7);
  TH1 *NN2Phi = new TH1F("NN2_phi", "NN2_phi", 100, -3.5, 3.5);

  //NN3 Histograms
  TH1 *NN3Status= new TH1F("NN3_Status", "NN3_Status", 100, 0, 10);
  TH1 *NN3E = new TH1F("NN3_E", "NN3_E", 100, 0, 2000);
  TH1 *NN3PT = new TH1F("NN3_PT", "NN3_PT", 100, 0, 1000);
  TH1 *NN3Eta = new TH1F("NN3_Eta", "NN3_Eta", 100, -7, 7);
  TH1 *NN3Phi = new TH1F("NN3_phi", "NN3_phi", 100, -3.5, 3.5);

  // Add all histograms to list
  l->Add(PartID);

  l->Add(chiStatus);
  l->Add(chiE);
  l->Add(chiPT);
  l->Add(chiEta);
  l->Add(chiPhi);

  l->Add(scalarStatus);
  l->Add(scalarE);
  l->Add(scalarPT);
  l->Add(scalarEta);
  l->Add(scalarPhi);

  l->Add(psiStatus);
  l->Add(psiE);
  l->Add(psiPT);
  l->Add(psiEta);
  l->Add(psiPhi);

  l->Add(JetPT);
  l->Add(JetEta);
  l->Add(JetPhi);

  l->Add(elecPT);
  l->Add(elecEta);
  l->Add(elecPhi);


  l->Add(MuonPT);
  l->Add(MuonEta);
  l->Add(MuonPhi);

  l->Add(NN1PT);
  l->Add(NN1Eta);
  l->Add(NN1Phi);

  l->Add(NN2PT);
  l->Add(NN2Eta);
  l->Add(NN2Phi);

  l->Add(NN3PT);
  l->Add(NN3Eta);
  l->Add(NN3Phi);
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);


    //HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    //LHEFEvent *event = (LHEFEvent*) branchEvent -> At(0);
    //Float_t weight = event->Weight;

      if (branchParticle->GetEntries() > 0)
      {
          cout << "Number of Particles: "<< branchParticle->GetEntries() << endl;
          for( unsigned int a = 0; a < branchParticle->GetEntries(); a = a + 1 )
          {
            GenParticle *particle =  (GenParticle*) branchParticle->At(a);
            PartID->Fill(particle->PID);
            cout << "PID: "<<particle->PID << endl;
            if (3==(particle->Status)) // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate
            {
            if (200001 == abs(particle->PID))
            {
              chiStatus->Fill(particle->Status);
              chiE->Fill(particle->E);
              chiPT->Fill(particle->PT);
              chiEta->Fill(particle->Eta);
              chiPhi->Fill(particle->Phi);
            }
            if (5000001 == abs(particle->PID))
            {
              scalarStatus->Fill(particle->Status);
              scalarE->Fill(particle->E);
              scalarPT->Fill(particle->PT);
              scalarEta->Fill(particle->Eta);
              scalarPhi->Fill(particle->Phi);
            }
            if (200000 == abs(particle->PID))
            {
              psiStatus->Fill(particle->Status);
              psiE->Fill(particle->E);
              psiPT->Fill(particle->PT);
              psiEta->Fill(particle->Eta);
              psiPhi->Fill(particle->Phi);
            }
            if (13 == abs(particle->PID))
            {
              MuonStatus->Fill(particle->Status);
              MuonE->Fill(particle->E);
              MuonPT->Fill(particle->PT);
              MuonEta->Fill(particle->Eta);
              MuonPhi->Fill(particle->Phi);
            }
            if (200002 == abs(particle->PID))
            {
              NN1Status->Fill(particle->Status);
              NN1E->Fill(particle->E);
              NN1PT->Fill(particle->PT);
              NN1Eta->Fill(particle->Eta);
              NN1Phi->Fill(particle->Phi);
            }
            if (200003 == abs(particle->PID))
            {
              NN2Status->Fill(particle->Status);
              NN2E->Fill(particle->E);
              NN2PT->Fill(particle->PT);
              NN2Eta->Fill(particle->Eta);
              NN2Phi->Fill(particle->Phi);
            }
            if (200004 == abs(particle->PID))
            {
              NN3Status->Fill(particle->Status);
              NN3E->Fill(particle->E);
              NN3PT->Fill(particle->PT);
              NN3Eta->Fill(particle->Eta);
              NN3Phi->Fill(particle->Phi);
            }

          }
          }
        }

        // If event contains at least 1 jet
        if(branchJet->GetEntries() > 0)
        {
          for( unsigned int b = 0; b < branchJet->GetEntries(); b = b + 1 )
          {
          // Take first jet
            Jet *jet = (Jet*) branchJet->At(b);
              ///if (jet != nullptr)

              // Plot jet transverse momentum
              JetPT->Fill(jet->PT);
              JetEta->Fill(jet->Eta);
              JetPhi->Fill(jet->Phi);
          
          }
        }

        // If event contains at least 1 electrons
        if(branchElectron->GetEntries() > 0)
        {
          for( unsigned int c = 0; c < branchElectron->GetEntries(); c = c + 1 )
          {
          // Take first two electrons
          Electron *elec = (Electron*) branchElectron->At(c);
          elecPT->Fill(elec->PT);
          elecEta->Fill(elec->Eta);
          elecPhi->Fill(elec->Phi);
          }
        }
        /*
        if(branchMuon->GetEntries() > 0)
        {
          for( unsigned int b = 0; b < sizeof(branchMuon); b = b + 1 )
          {
          // Take first Muon
            Muon *muon = (Muon*) branchMuon->At(b);
              if (muon != nullptr)
              {
              // Plot Muon transverse momentum
              MuonPT->Fill(muon->PT);
              MuonEta->Fill(muon->Eta);
              MuonPhi->Fill(muon->Phi);
          }
          }
        }
        */





  }


  TFile *f = new TFile("ntuple_delphes.root","RECREATE");
  l->Write("ntuple_delphes", TObject::kSingleKey);
  f->ls();

}
