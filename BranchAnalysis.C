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
  TH1 *chiStatus= new TH1F("chi_Status", "chi_Status", 100, 0, 100);
  TH1 *chiPT = new TH1F("chi_PT", "chi_PT", 100, 0, 1000000);
  TH1 *chiE = new TH1F("chi_E", "chi_E", 100, 0, 1000000);
  TH1 *chiEta = new TH1F("chi_Eta", "chi_Eta", 100, 0, 7);
  TH1 *chiPhi = new TH1F("chi_phi", "chi_phi", 100, 0, 7);

  // Add all histograms to list
  l->Add(histJetPT);
  l->Add(histMass);
  l->Add(PartID);

  l->Add(chiStatus);
  l->Add(chiE);
  l->Add(chiPT);
  l->Add(chiEta);
  l->Add(chiPhi);

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
          GenParticle *particle =  (GenParticle*) branchParticle->At(0);
          PartID->Fill(particle->PID);
          cout << "PID: "<<particle->PID << endl;
          if (1000024 == abs(particle->PID))
          {
            chiStatus->Fill(particle->Status);
            chiE->Fill(particle->E);
            chiPT->Fill(particle->PT);
            chiEta->Fill(particle->Eta);
            chiPhi->Fill(particle->Phi);
          }
        }


    // If event contains at least 1 jet
    if(branchJet->GetEntries() > 0)
    {
      // Take first jet
      Jet *jet = (Jet*) branchJet->At(0);

      // Plot jet transverse momentum
      histJetPT->Fill(jet->PT);

      // Print jet transverse momentum
      //cout << "Jet pt: "<<jet->PT << endl;
    }

    Electron *elec1, *elec2;

    // If event contains at least 2 electrons
    if(branchElectron->GetEntries() > 1)
    {
      // Take first two electrons
      elec1 = (Electron *) branchElectron->At(0);
      elec2 = (Electron *) branchElectron->At(1);

      // Plot their invariant mass
      histMass->Fill(((elec1->P4()) + (elec2->P4())).M());
    }
  }


  TFile *f = new TFile("histlist.root","RECREATE");
  l->Write("histlist", TObject::kSingleKey);
  f->ls();

}
