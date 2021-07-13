/*
Simple macro showing how to access branches from the delphes output root file,
loop over events.
root -l examples/BranchAnalysis.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include <iostream>
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#endif

//------------------------------------------------------------------------------
double Delta_R_Function(TLorentzVector c1,TLorentzVector c2,TLorentzVector m1,TLorentzVector m2 )
{
  float DR_pr1 = c1.DeltaR(m1);
  float DR_pr2 = c1.DeltaR(m2);
  float DR_pr3 = c2.DeltaR(m1);
  float DR_pr4 = c2.DeltaR(m2);
  if (DR_pr1 + DR_pr4 <  DR_pr2 + DR_pr3)
  {
    if (DR_pr1 < DR_pr4){
      return DR_pr1;
    }
    else
    {
      return DR_pr4;
    }
  }
  else
  {
    if (DR_pr2 < DR_pr3){
      return DR_pr2;
    }
    else
    {
      return DR_pr3;
    }
  }
}

void AnalysisCPP(const char *inputFile)
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




  // Analysis Tree Variables
  vector<int> chi_status_int;
  vector<float> chi_pt_float;
  vector<float> chi_E_float;
  vector<float> chi_Eta_float;
  vector<float> chi_Phi_float;

  vector<int> scalar_status_int;
  vector<float> scalar_pt_float;
  vector<float> scalar_E_float;
  vector<float> scalar_Eta_float;
  vector<float> scalar_Phi_float;

  vector<int> psi_status_int;
  vector<float> psi_pt_float;
  vector<float> psi_E_float;
  vector<float> psi_Eta_float;
  vector<float> psi_Phi_float;

  vector<float> Jet_pt_float;
  vector<float> Jet_Eta_float;
  vector<float> Jet_Phi_float;

  vector<int> elec_N_int;
  vector<int> elec_status_int;
  vector<float> elec_pt_float;
  vector<float> elec_Eta_float;
  vector<float> elec_Phi_float;

  vector<int> Muon_status_int;
  vector<int> Muon_N_int;
  vector<float> Muon_pt_float;
  vector<float> Muon_Eta_float;
  vector<float> Muon_Phi_float;
  vector<float> Muon_E_float;

  vector<int> tau_status_int;
  vector<int> tau_N_int;
  vector<float> tau_pt_float;
  vector<float> tau_Eta_float;
  vector<float> tau_Phi_float;

  vector<float> NN1_pt_float;
  vector<float> NN1_Eta_float;
  vector<float> NN1_Phi_float;

  vector<float> NN2_pt_float;
  vector<float> NN2_Eta_float;
  vector<float> NN2_Phi_float;

  vector<float> NN3_pt_float;
  vector<float> NN3_Eta_float;
  vector<float> NN3_Phi_float;

  // Delta R Analysis Variable
  vector<float> Muon_DeltaR;



  // Book histograms
  TFile *f = new TFile("tcl.root","RECREATE");
  auto T = new TTree("analysisTree","test");
  auto T_Delta = new TTree("Delta_Analysis","test");

  // Analysis Tree Branches
  T->Branch("chi_status",&chi_status_int);
  T->Branch("chi_pt",&chi_pt_float);
  T->Branch("chi_E",&chi_E_float);
  T->Branch("chi_Eta",&chi_Eta_float);
  T->Branch("chi_Phi",&chi_Phi_float);

  T->Branch("scalar_status",&scalar_status_int);
  T->Branch("scalar_pt",&scalar_pt_float);
  T->Branch("scalar_E",&scalar_E_float);
  T->Branch("scalar_Eta",&scalar_Eta_float);
  T->Branch("scalar_Phi",&scalar_Phi_float);

  T->Branch("psi_status",&psi_status_int);
  T->Branch("psi_pt",&psi_pt_float);
  T->Branch("psi_E",&psi_E_float);
  T->Branch("psi_Eta",&psi_Eta_float);
  T->Branch("psi_Phi",&psi_Phi_float);

  T->Branch("Jet_pt",&Jet_pt_float);
  T->Branch("Jet_Eta",&Jet_Eta_float);
  T->Branch("Jet_Phi",&Jet_Phi_float);

  T->Branch("elec_Status",&elec_status_int);
  T->Branch("elec_N",&elec_N_int);
  T->Branch("elec_pt",&elec_pt_float);
  T->Branch("elec_Eta",&elec_Eta_float);
  T->Branch("elec_Phi",&elec_Phi_float);

  T->Branch("Muon_Status",&Muon_status_int);
  T->Branch("Muon_N",&Muon_N_int);
  T->Branch("Muon_pt",&Muon_pt_float);
  T->Branch("Muon_Eta",&Muon_Eta_float);
  T->Branch("Muon_Phi",&Muon_Phi_float);

  T->Branch("tau_Status",&tau_status_int);
  T->Branch("tau_N",&tau_N_int);
  T->Branch("tau_pt",&tau_pt_float);
  T->Branch("tau_Eta",&tau_Eta_float);
  T->Branch("tau_Phi",&tau_Phi_float);

  T->Branch("NN1_pt",&NN1_pt_float);
  T->Branch("NN1_Eta",&NN1_Eta_float);
  T->Branch("NN1_Phi",&NN1_Phi_float);

  T->Branch("NN2_pt",&NN2_pt_float);
  T->Branch("NN2_Eta",&NN2_Eta_float);
  T->Branch("NN2_Phi",&NN2_Phi_float);

  T->Branch("NN3_pt",&NN3_pt_float);
  T->Branch("NN3_Eta",&NN3_Eta_float);
  T->Branch("NN3_Phi",&NN3_Phi_float);


  // Delta R Analysis Branches

  T_Delta->Branch("Muon_Delta_R", &Muon_DeltaR);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
      if (branchParticle->GetEntries() > 0)
      {
          //cout << "Number of Particles: "<< branchParticle->GetEntries() << endl;
          for( unsigned int a = 0; a < branchParticle->GetEntries(); a = a + 1 )
          {
            GenParticle *particle =  (GenParticle*) branchParticle->At(a);
            //PartID->Fill(particle->PID);
            //cout << "PID: "<<particle->PID << endl;
            //if (3==(particle->Status)) // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate
            //{
            if (200001 == abs(particle->PID))
            {
              chi_status_int.push_back(particle->Status);
              if (3 == (particle->Status))
              {
              chi_pt_float.push_back(particle->PT);
              chi_E_float.push_back(particle->E);
              chi_Eta_float.push_back(particle->Eta);
              chi_Phi_float.push_back(particle->Phi);
            }
            }
            if (5000001 == abs(particle->PID))
            {
              scalar_status_int.push_back(particle->Status);
              if (3 == (particle->Status))
              {
              scalar_pt_float.push_back(particle->PT);
              scalar_E_float.push_back(particle->E);
              scalar_Eta_float.push_back(particle->Eta);
              scalar_Phi_float.push_back(particle->Phi);
            }
            }
            if (200000 == abs(particle->PID))
            {
              psi_status_int.push_back(particle->Status);
              if (3 == (particle->Status))
              {
              psi_pt_float.push_back(particle->PT);
              psi_E_float.push_back(particle->E);
              psi_Eta_float.push_back(particle->Eta);
              psi_Phi_float.push_back(particle->Phi);
              }
            }
            if (13 == abs(particle->PID))
            {
              Muon_status_int.push_back(particle->Status);
              Muon_pt_float.push_back(particle->PT);
              Muon_Eta_float.push_back(particle->Eta);
              Muon_Phi_float.push_back(particle->Phi);
              Muon_E_float.push_back(particle->E);
            }
            if (15 == abs(particle->PID))
            {
              tau_status_int.push_back(particle->Status);
              tau_pt_float.push_back(particle->PT);
              tau_Eta_float.push_back(particle->Eta);
              tau_Phi_float.push_back(particle->Phi);
            }
            if (11 == abs(particle->PID))
            {
              elec_status_int.push_back(particle->Status);
              elec_pt_float.push_back(particle->PT);
              elec_Eta_float.push_back(particle->Eta);
              elec_Phi_float.push_back(particle->Phi);
            }
            if (200002 == abs(particle->PID))
            {
              NN1_pt_float.push_back(particle->PT);
              NN1_Eta_float.push_back(particle->Eta);
              NN1_Phi_float.push_back(particle->Phi);
            }
            if (200003 == abs(particle->PID))
            {
              NN2_pt_float.push_back(particle->PT);
              NN2_Eta_float.push_back(particle->Eta);
              NN2_Phi_float.push_back(particle->Phi);
            }
            if (200004 == abs(particle->PID))
            {
              NN3_pt_float.push_back(particle->PT);
              NN3_Eta_float.push_back(particle->Eta);
              NN3_Phi_float.push_back(particle->Phi);
            }
          }
          //}
          ////////
        }
        if(branchJet->GetEntries() > 0)
        {
          //cout << "Number of Jets: "<< branchJet->GetEntries() << endl;
          for( unsigned int b = 0; b < branchJet->GetEntries(); b = b + 1 )
          {
          // Take first jet
            Jet *jet = (Jet*) branchJet->At(b);
              // Plot jet transverse momentum
              Jet_pt_float.push_back(jet->PT);
              Jet_Eta_float.push_back(jet->Eta);
              Jet_Phi_float.push_back(jet->Phi);
          }
        }
        // Calculate the number of particles per event
        if (Muon_pt_float.size() > 0){
          Muon_N_int.push_back(Muon_pt_float.size());
        }
        if (elec_pt_float.size() > 0){
          elec_N_int.push_back(elec_pt_float.size());
        }
        if (tau_pt_float.size() > 0){
          tau_N_int.push_back(tau_pt_float.size());
        }

        // Delta R Analysis
        if ( Muon_pt_float.size() == 2 and chi_pt_float.size() == 2 )
        {
          //cout << "Creating TLorentzVector"<< endl;
          TLorentzVector c1;
          TLorentzVector c2;
          TLorentzVector m1;
          TLorentzVector m2;

          c1.SetPtEtaPhiE(chi_pt_float[0],chi_Eta_float[0],chi_Phi_float[0],chi_E_float[0]);
          c2.SetPtEtaPhiE(chi_pt_float[1],chi_Eta_float[1],chi_Phi_float[1],chi_E_float[1]);
          m1.SetPtEtaPhiE(Muon_pt_float[0],Muon_Eta_float[0],Muon_Phi_float[0],Muon_E_float[0]);
          m2.SetPtEtaPhiE(Muon_pt_float[1],Muon_Eta_float[1],Muon_Phi_float[1],Muon_E_float[1]);

          Muon_DeltaR.push_back(Delta_R_Function(c1,c2,m1,m2));
          }



        T->Fill();
        T_Delta->Fill();

        chi_pt_float.clear();
        chi_status_int.clear();
        chi_pt_float.clear();
        chi_E_float.clear();
        chi_Eta_float.clear();
        chi_Phi_float.clear();

        scalar_status_int.clear();
        scalar_pt_float.clear();
        scalar_E_float.clear();
        scalar_Eta_float.clear();
        scalar_Phi_float.clear();

        psi_status_int.clear();
        psi_pt_float.clear();
        psi_E_float.clear();
        psi_Eta_float.clear();
        psi_Phi_float.clear();

        Jet_pt_float.clear();
        Jet_Eta_float.clear();
        Jet_Phi_float.clear();

        elec_status_int.clear();
        elec_N_int.clear();
        elec_pt_float.clear();
        elec_Eta_float.clear();
        elec_Phi_float.clear();


        Muon_status_int.clear();
        Muon_N_int.clear();
        Muon_pt_float.clear();
        Muon_Eta_float.clear();
        Muon_Phi_float.clear();
        Muon_E_float.clear();

        tau_status_int.clear();
        tau_N_int.clear();
        tau_pt_float.clear();
        tau_Eta_float.clear();
        tau_Phi_float.clear();

        NN1_pt_float.clear();
        NN1_Eta_float.clear();
        NN1_Phi_float.clear();

        NN2_pt_float.clear();
        NN2_Eta_float.clear();
        NN2_Phi_float.clear();

        NN3_pt_float.clear();
        NN3_Eta_float.clear();
        NN3_Phi_float.clear();


        Muon_DeltaR.clear();
  }
  T->Print();
  T_Delta->Print();
  f->Write();
}