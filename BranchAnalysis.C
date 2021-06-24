/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.
root -l examples/Example1.C'("delphes_output.root")'
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

  // Book histograms
  TList *l = new TList();
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0);
  TH1 *histMass = new TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0);
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

  // Add all histograms to list
  l->Add(histJetPT);
  l->Add(histMass);
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
          for( unsigned int a = 0; a < sizeof(branchParticle); a = a + 1 )
          {
            GenParticle *particle =  (GenParticle*) branchParticle->At(a);
            PartID->Fill(particle->PID);
            cout << "PID: "<<particle->PID << endl;
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
          }
        }

  }


  TFile *f = new TFile("ntuple_delphes.root","RECREATE");
  l->Write("ntuple_delphes", TObject::kSingleKey);
  f->ls();

}
