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
float PairsDeltaR(vector<GenParticle*> pt1 , vector<GenParticle*> pt2){
  // Calculates the delta R values between pairs of particles
      TLorentzVector c1;
      TLorentzVector c2;
      TLorentzVector m1;
      TLorentzVector m2;
      //cout << pt1[0]->Mass << endl;
      c1.SetPtEtaPhiE(pt1[0]->PT,pt1[0]->Eta,pt1[0]->Phi,pt1[0]->Mass);
      c2.SetPtEtaPhiE(pt1[1]->PT,pt1[1]->Eta,pt1[1]->Phi,pt1[1]->Mass);
      m1.SetPtEtaPhiE(pt2[0]->PT,pt2[0]->Eta,pt2[0]->Phi,pt2[0]->Mass);
      m2.SetPtEtaPhiE(pt2[1]->PT,pt2[1]->Eta,pt2[1]->Phi,pt2[1]->Mass);
      float DR_pr1 = c1.DeltaR(m1);
      float DR_pr2 = c1.DeltaR(m2);
      float DR_pr3 = c2.DeltaR(m1);
      float DR_pr4 = c2.DeltaR(m2);
      if (DR_pr1 + DR_pr4 <  DR_pr2 + DR_pr3)
      {
        if (DR_pr1 < DR_pr4){return DR_pr1;}
        else  {return DR_pr4;}
      }
      else
      {
        if (DR_pr2 < DR_pr3){return DR_pr2;}
        else   {  return DR_pr3;  }
      }
}

float PairsMuonDeltaR(vector<Muon*> pt1 , vector<GenParticle*> pt2){
  // Calculates the delta R values between pairs of particles
      TLorentzVector c1;
      TLorentzVector c2;
      TLorentzVector m1;
      TLorentzVector m2;

      c1.SetPtEtaPhiE(pt1[0]->PT,pt1[0]->Eta,pt1[0]->Phi,0.10566);
      c2.SetPtEtaPhiE(pt1[1]->PT,pt1[1]->Eta,pt1[1]->Phi,0.10566);
      m1.SetPtEtaPhiE(pt2[0]->PT,pt2[0]->Eta,pt2[0]->Phi,pt2[0]->Mass);
      m2.SetPtEtaPhiE(pt2[1]->PT,pt2[1]->Eta,pt2[1]->Phi,pt2[1]->Mass);
      float DR_pr1 = c1.DeltaR(m1);
      float DR_pr2 = c1.DeltaR(m2);
      float DR_pr3 = c2.DeltaR(m1);
      float DR_pr4 = c2.DeltaR(m2);
      if (DR_pr1 + DR_pr4 <  DR_pr2 + DR_pr3)
      {
        if (DR_pr1 < DR_pr4){return DR_pr1;}
        else  {return DR_pr4;}
      }
      else
      {
        if (DR_pr2 < DR_pr3){return DR_pr2;}
        else   {  return DR_pr3;  }
      }
}

float PairsElectronDeltaR(vector<Electron*> pt1 , vector<GenParticle*> pt2){
  // Calculates the delta R values between pairs of particles
      TLorentzVector c1;
      TLorentzVector c2;
      TLorentzVector m1;
      TLorentzVector m2;
      c1.SetPtEtaPhiE(pt1[0]->PT,pt1[0]->Eta,pt1[0]->Phi,0.000511);
      c2.SetPtEtaPhiE(pt1[1]->PT,pt1[1]->Eta,pt1[1]->Phi,0.000511);
      m1.SetPtEtaPhiE(pt2[0]->PT,pt2[0]->Eta,pt2[0]->Phi,pt2[0]->Mass);
      m2.SetPtEtaPhiE(pt2[1]->PT,pt2[1]->Eta,pt2[1]->Phi,pt2[1]->Mass);
      float DR_pr1 = c1.DeltaR(m1);
      float DR_pr2 = c1.DeltaR(m2);
      float DR_pr3 = c2.DeltaR(m1);
      float DR_pr4 = c2.DeltaR(m2);
      if (DR_pr1 + DR_pr4 <  DR_pr2 + DR_pr3)
      {
        if (DR_pr1 < DR_pr4){return DR_pr1;}
        else  {return DR_pr4;}
      }
      else
      {
        if (DR_pr2 < DR_pr3){return DR_pr2;}
        else   {  return DR_pr3;  }
      }
}
vector<float> TrackDeltaR(vector<GenParticle*> pt1 , vector<Track*> trks){
  // Function that calculates the Delta R value between the charged particles (chi/scalar) and the detector tracks
  vector<float> TrackDeltaR_Values;
  TrackDeltaR_Values.clear();
  for( unsigned int a = 0; a < pt1.size(); a = a + 1 ){
    for( unsigned int b = 0; b < trks.size(); b = b + 1 ){
      int ptq_charge = pt1[a]->PID/abs(pt1[a]->PID);
      //cout << "Tracks charge: "<< trks[b]->Charge << " particle charge: " << ptq_charge << " , PID: "<< pt1[a]->PID << endl;
      if (ptq_charge == trks[b]->Charge){
        TLorentzVector particle1;
        TLorentzVector track;
        particle1.SetPtEtaPhiM(pt1[a]->PT,pt1[a]->Eta,pt1[a]->Phi,pt1[a]->Mass);
        track.SetPtEtaPhiM(trks[b]->PT,trks[b]->Eta,trks[b]->Phi,trks[b]->Mass);
        if (particle1.DeltaR(track) < 1){
          TrackDeltaR_Values.push_back(particle1.DeltaR(track));
          //cout << trks[b]->PID << endl;
        }
      }
    }
  }
  return TrackDeltaR_Values;
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
  TClonesArray *branchTrack = treeReader->UseBranch("Track");

  // Analysis Tree Variables
  vector<GenParticle*> chis;
  vector<int> chi_status_int;
  vector<float> chi_pt_float;
  vector<float> chi_E_float;
  vector<float> chi_Eta_float;
  vector<float> chi_Phi_float;

  vector<GenParticle*> scalars;
  vector<int> scalar_status_int;
  vector<float> scalar_pt_float;
  vector<float> scalar_E_float;
  vector<float> scalar_Eta_float;
  vector<float> scalar_Phi_float;

  vector<GenParticle*> psis;
  vector<int> psi_status_int;
  vector<float> psi_pt_float;
  vector<float> psi_E_float;
  vector<float> psi_Eta_float;
  vector<float> psi_Phi_float;

  vector<GenParticle*> elecs_gen;
  vector<int> elec_Gen_N_int;
  vector<int> elec_Gen_status_int;
  vector<float> elec_Gen_E_float;
  vector<float> elec_Gen_pt_float;
  vector<float> elec_Gen_Eta_float;
  vector<float> elec_Gen_Phi_float;

  vector<GenParticle*> muons_gen;
  vector<int> Muon_Gen_status_int;
  vector<int> Muon_Gen_N_int;
  vector<float> Muon_Gen_pt_float;
  vector<float> Muon_Gen_Eta_float;
  vector<float> Muon_Gen_Phi_float;
  vector<float> Muon_Gen_E_float;

  vector<GenParticle*> taus_gen;
  vector<int> tau_Gen_status_int;
  vector<int> tau_Gen_N_int;
  vector<float> tau_Gen_E_float;
  vector<float> tau_Gen_pt_float;
  vector<float> tau_Gen_Eta_float;
  vector<float> tau_Gen_Phi_float;

  vector<float> NN1_pt_float;
  vector<float> NN1_Eta_float;
  vector<float> NN1_Phi_float;

  vector<float> NN2_pt_float;
  vector<float> NN2_Eta_float;
  vector<float> NN2_Phi_float;

  vector<float> NN3_pt_float;
  vector<float> NN3_Eta_float;
  vector<float> NN3_Phi_float;

  // Detector Variables
  vector<Electron*> elecs;
  vector<int> elec_N_int;
  vector<float> elec_pt_float;
  vector<float> elec_Eta_float;
  vector<float> elec_Phi_float;

  vector<Muon*> muons;
  vector<int> muon_N_int;
  vector<float> Muon_pt_float;
  vector<float> Muon_Eta_float;
  vector<float> Muon_Phi_float;

  vector<float> Jet_pt_float;
  vector<float> Jet_Eta_float;
  vector<float> Jet_Phi_float;

  vector<float> Track_pt_float;
  vector<float> Track_Eta_float;
  vector<float> Track_Phi_float;

  // Delta R Analysis Variable
  vector<float> chi_Muon_DeltaR;
  vector<float> chi_elec_DeltaR;
  vector<float> chi_tau_DeltaR;

  vector<float> scalar_Muon_DeltaR;
  vector<float> scalar_elec_DeltaR;
  vector<float> scalar_tau_DeltaR;

  vector<float> Psi_Muon_DeltaR;
  vector<float> Psi_elec_DeltaR;
  vector<float> Psi_tau_DeltaR;

  vector<float> Muon_Muon_Gen_Delta_R;
  vector<float> elec_elec_Gen_Delta_R;

  vector<Track*> tracks;

  vector<float> chi_Track_DeltaR;
  vector<float> scalar_Track_DeltaR;


  // Book histograms
  TFile *f = new TFile("ntuple_delphes.root","RECREATE");
  auto T = new TTree("GenParticles","test");
  auto T_Det = new TTree("Particles","test");
  auto T_Delta = new TTree("DeltaR_Analysis","test");

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

  T->Branch("elec_Gen_Status",&elec_Gen_status_int);
  T->Branch("elec_Gen_N",&elec_Gen_N_int);
  T->Branch("elec_Gen_pt",&elec_Gen_pt_float);
  T->Branch("elec_Gen_Eta",&elec_Gen_Eta_float);
  T->Branch("elec_Gen_Phi",&elec_Gen_Phi_float);

  T->Branch("Muon_Gen_Status",&Muon_Gen_status_int);
  T->Branch("Muon_Gen_N",&Muon_Gen_N_int);
  T->Branch("Muon_Gen_pt",&Muon_Gen_pt_float);
  T->Branch("Muon_Gen_Eta",&Muon_Gen_Eta_float);
  T->Branch("Muon_Gen_Phi",&Muon_Gen_Phi_float);

  T->Branch("tau_Gen_Status",&tau_Gen_status_int);
  T->Branch("tau_Gen_N",&tau_Gen_N_int);
  T->Branch("tau_Gen_pt",&tau_Gen_pt_float);
  T->Branch("tau_Gen_Eta",&tau_Gen_Eta_float);
  T->Branch("tau_Gen_Phi",&tau_Gen_Phi_float);

  T->Branch("NN1_pt",&NN1_pt_float);
  T->Branch("NN1_Eta",&NN1_Eta_float);
  T->Branch("NN1_Phi",&NN1_Phi_float);

  T->Branch("NN2_pt",&NN2_pt_float);
  T->Branch("NN2_Eta",&NN2_Eta_float);
  T->Branch("NN2_Phi",&NN2_Phi_float);

  T->Branch("NN3_pt",&NN3_pt_float);
  T->Branch("NN3_Eta",&NN3_Eta_float);
  T->Branch("NN3_Phi",&NN3_Phi_float);



  // Detector Branches
  T_Det->Branch("elec_pt",&elec_pt_float);
  T_Det->Branch("elec_N",&elec_N_int);
  T_Det->Branch("elec_Eta",&elec_Eta_float);
  T_Det->Branch("elec_Phi",&elec_Phi_float);

  T_Det->Branch("Muon_pt",&Muon_pt_float);
  T_Det->Branch("Muon_N",&muon_N_int);
  T_Det->Branch("Muon_Eta",&Muon_Eta_float);
  T_Det->Branch("Muon_Phi",&Muon_Phi_float);


  T_Det->Branch("Jet_pt",&Jet_pt_float);
  T_Det->Branch("Jet_Eta",&Jet_Eta_float);
  T_Det->Branch("Jet_Phi",&Jet_Phi_float);

  T_Det->Branch("Track_pt",&Track_pt_float);
  T_Det->Branch("Track_Eta",&Track_Eta_float);
  T_Det->Branch("Track_Phi",&Track_Phi_float);



  // Delta R Analysis Branches
  T_Delta->Branch("Chi_Muon_Delta_R", &chi_Muon_DeltaR);
  T_Delta->Branch("Chi_elec_Delta_R", &chi_elec_DeltaR);
  T_Delta->Branch("Chi_tau_Delta_R", &chi_tau_DeltaR);

  T_Delta->Branch("scalar_Muon_Delta_R", &scalar_Muon_DeltaR);
  T_Delta->Branch("scalar_elec_Delta_R", &scalar_elec_DeltaR);
  T_Delta->Branch("scalar_tau_Delta_R", &scalar_tau_DeltaR);

  T_Delta->Branch("Psi_Muon_Delta_R", &Psi_Muon_DeltaR);
  T_Delta->Branch("Psi_elec_Delta_R", &Psi_elec_DeltaR);
  T_Delta->Branch("Psi_tau_Delta_R", &Psi_tau_DeltaR);

  T_Delta->Branch("Muon_Muon_Gen_Delta_R", &Muon_Muon_Gen_Delta_R);
  T_Delta->Branch("elec_elec_Gen_Delta_R", &elec_elec_Gen_Delta_R);

  T_Delta->Branch("Chi_Track_Delta_R", &chi_Track_DeltaR);
  T_Delta->Branch("Scalar_Track_Delta_R", &scalar_Track_DeltaR);


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
            // Chi
            if (200001 == abs(particle->PID))
            {
              chi_status_int.push_back(particle->Status);
              if (2 == (particle->Status)) // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate
              {
              chis.push_back(particle);
              chi_pt_float.push_back(particle->PT);
              chi_E_float.push_back(particle->E);
              chi_Eta_float.push_back(particle->Eta);
              chi_Phi_float.push_back(particle->Phi);
              }
            }
            // Scalar
            if (5000001 == abs(particle->PID))
            {
              scalar_status_int.push_back(particle->Status);
              if (2 == (particle->Status)) // For Three Body Process Change to 2, for Two body 1
              {
              scalars.push_back(particle);
              scalar_pt_float.push_back(particle->PT);
              scalar_E_float.push_back(particle->E);
              scalar_Eta_float.push_back(particle->Eta);
              scalar_Phi_float.push_back(particle->Phi);
              }
            }
            // Psi
            if (200000 == abs(particle->PID))
            {
              psi_status_int.push_back(particle->Status);
              if (1 == (particle->Status))
              {
              psis.push_back(particle);
              psi_pt_float.push_back(particle->PT);
              psi_E_float.push_back(particle->E);
              psi_Eta_float.push_back(particle->Eta);
              psi_Phi_float.push_back(particle->Phi);
              }
            }
            // Muon
            if (13 == abs(particle->PID))
            {

              muons_gen.push_back(particle);
              Muon_Gen_status_int.push_back(particle->Status);
              Muon_Gen_pt_float.push_back(particle->PT);
              Muon_Gen_Eta_float.push_back(particle->Eta);
              Muon_Gen_Phi_float.push_back(particle->Phi);
              Muon_Gen_E_float.push_back(particle->E);
            }
            // Tau
            if (15 == abs(particle->PID))
            {
              taus_gen.push_back(particle);
              tau_Gen_status_int.push_back(particle->Status);
              tau_Gen_pt_float.push_back(particle->PT);
              tau_Gen_Eta_float.push_back(particle->Eta);
              tau_Gen_Phi_float.push_back(particle->Phi);
              tau_Gen_E_float.push_back(particle->E);
            }
            // Electron
            if (11 == abs(particle->PID))
            {
              elecs_gen.push_back(particle);
              elec_Gen_status_int.push_back(particle->Status);
              elec_Gen_pt_float.push_back(particle->PT);
              elec_Gen_Eta_float.push_back(particle->Eta);
              elec_Gen_Phi_float.push_back(particle->Phi);
              elec_Gen_E_float.push_back(particle->E);
            }
            // Heavy Neutrino Electron
            if (200002 == abs(particle->PID))
            {
              NN1_pt_float.push_back(particle->PT);
              NN1_Eta_float.push_back(particle->Eta);
              NN1_Phi_float.push_back(particle->Phi);
            }
            // Heavy Neutrino Muon
            if (200003 == abs(particle->PID))
            {
              NN2_pt_float.push_back(particle->PT);
              NN2_Eta_float.push_back(particle->Eta);
              NN2_Phi_float.push_back(particle->Phi);
            }
            // Heavy Neutrino Tau
            if (200004 == abs(particle->PID))
            {
              NN3_pt_float.push_back(particle->PT);
              NN3_Eta_float.push_back(particle->Eta);
              NN3_Phi_float.push_back(particle->Phi);
            }
          }
        }

        // Non Generator Particles
        if(branchJet->GetEntries() > 0)  {
          for( unsigned int b = 0; b < branchJet->GetEntries(); b = b + 1 ){
            //cout << "Jet" << endl;
            Jet *jet = (Jet*) branchJet->At(b);
            Jet_pt_float.push_back(jet->PT);
            Jet_Eta_float.push_back(jet->Eta);
            Jet_Phi_float.push_back(jet->Phi);
          }
        }
        if(branchMuon->GetEntries() > 0)  {
          for( unsigned int q = 0; q < branchMuon->GetEntries(); q = q + 1 ){
            //cout << "Muon" << endl;
            Muon *muon = (Muon*) branchMuon->At(q);
            muons.push_back(muon);
            Muon_pt_float.push_back(muon->PT);
            Muon_Eta_float.push_back(muon->Eta);
            Muon_Phi_float.push_back(muon->Phi);
          }
        }
        if(branchElectron->GetEntries() > 0)  {
          for( unsigned int b = 0; b < branchElectron->GetEntries(); b = b + 1 ){
            Electron *elec = (Electron*) branchElectron->At(b);
            elecs.push_back(elec);
            elec_pt_float.push_back(elec->PT);
            elec_Eta_float.push_back(elec->Eta);
            elec_Phi_float.push_back(elec->Phi);
          }
        }
        if(branchTrack->GetEntries()> 0 ){
          for( unsigned int b = 0; b < branchTrack->GetEntries(); b = b + 1 ){
            Track *trk = (Track*) branchTrack->At(b);
            tracks.push_back(trk);
            Track_pt_float.push_back(trk->PT);
            Track_Eta_float.push_back(trk->Eta);
            Track_Phi_float.push_back(trk->Phi);
          }
        }

        // Calculate the number of particles per event

        Muon_Gen_N_int.push_back(Muon_Gen_pt_float.size());
        elec_Gen_N_int.push_back(elec_Gen_pt_float.size());
        tau_Gen_N_int.push_back(tau_Gen_pt_float.size());
        muon_N_int.push_back(Muon_pt_float.size());
        elec_N_int.push_back(elec_pt_float.size());

        // Delta R Analysis
        // Chis
        if (chis.size() == 2 and muons_gen.size() == 2 ){
          chi_Muon_DeltaR.push_back(PairsDeltaR(chis,muons_gen));
        }
        if (chis.size() == 2 and elecs_gen.size() == 2 ){
          chi_elec_DeltaR.push_back(PairsDeltaR(chis,elecs_gen));
        }
        if (chis.size() == 2 and taus_gen.size() == 2 ){
          chi_tau_DeltaR.push_back(PairsDeltaR(chis,taus_gen));
        }
        // Scalars
        if (scalars.size() == 2 and muons_gen.size() == 2 ){
          scalar_Muon_DeltaR.push_back(PairsDeltaR(scalars,muons_gen));
        }
        if (scalars.size() == 2 and elecs_gen.size() == 2 ){
          scalar_elec_DeltaR.push_back(PairsDeltaR(scalars,elecs_gen));
        }
        if (scalars.size() == 2 and taus_gen.size() == 2 ){
          scalar_tau_DeltaR.push_back(PairsDeltaR(scalars,taus_gen));
        }
        // Psis
        if (psis.size() == 2 and muons_gen.size() == 2 ){
          Psi_Muon_DeltaR.push_back(PairsDeltaR(psis,muons_gen));
        }
        if (psis.size() == 2 and elecs_gen.size() == 2 ){
          Psi_elec_DeltaR.push_back(PairsDeltaR(psis,elecs_gen));
        }
        if (psis.size() == 2 and taus_gen.size() == 2 ){
          Psi_tau_DeltaR.push_back(PairsDeltaR(psis,taus_gen));
        }


        // Lepton Pair Matching

        if (muons.size() == 2 and muons_gen.size() == 2 ){
          Muon_Muon_Gen_Delta_R.push_back(PairsMuonDeltaR(muons,muons_gen));
        }
        if (elecs.size() == 2 and elecs_gen.size() == 2 ){
          elec_elec_Gen_Delta_R.push_back(PairsElectronDeltaR(elecs,elecs_gen));
        }

        // Calculates the deltaR value between chis, scalars and the tracks
        //cout << "Number of Tracks: " << tracks.size() << endl;
        chi_Track_DeltaR = TrackDeltaR(chis,tracks);
        scalar_Track_DeltaR = TrackDeltaR(scalars,tracks);

        // Fill the histograms
        T->Fill();
        T_Det->Fill();
        T_Delta->Fill();

        // GenParticle Clears
        chis.clear();
        chi_pt_float.clear();
        chi_status_int.clear();
        chi_pt_float.clear();
        chi_E_float.clear();
        chi_Eta_float.clear();
        chi_Phi_float.clear();

        scalars.clear();
        scalar_status_int.clear();
        scalar_pt_float.clear();
        scalar_E_float.clear();
        scalar_Eta_float.clear();
        scalar_Phi_float.clear();

        psis.clear();
        psi_status_int.clear();
        psi_pt_float.clear();
        psi_E_float.clear();
        psi_Eta_float.clear();
        psi_Phi_float.clear();

        elecs_gen.clear();
        elec_Gen_status_int.clear();
        elec_Gen_N_int.clear();
        elec_Gen_pt_float.clear();
        elec_Gen_Eta_float.clear();
        elec_Gen_Phi_float.clear();
        elec_Gen_E_float.clear();

        muons_gen.clear();
        Muon_Gen_status_int.clear();
        Muon_Gen_N_int.clear();
        Muon_Gen_pt_float.clear();
        Muon_Gen_Eta_float.clear();
        Muon_Gen_Phi_float.clear();
        Muon_Gen_E_float.clear();

        taus_gen.clear();
        tau_Gen_status_int.clear();
        tau_Gen_N_int.clear();
        tau_Gen_pt_float.clear();
        tau_Gen_Eta_float.clear();
        tau_Gen_Phi_float.clear();
        tau_Gen_E_float.clear();

        NN1_pt_float.clear();
        NN1_Eta_float.clear();
        NN1_Phi_float.clear();

        NN2_pt_float.clear();
        NN2_Eta_float.clear();
        NN2_Phi_float.clear();

        NN3_pt_float.clear();
        NN3_Eta_float.clear();
        NN3_Phi_float.clear();

        //Detector Clears
        muons.clear();
        muon_N_int.clear();
        Muon_pt_float.clear();
        Muon_Eta_float.clear();
        Muon_Phi_float.clear();

        elecs.clear();
        elec_N_int.clear();
        elec_pt_float.clear();
        elec_Eta_float.clear();
        elec_Phi_float.clear();

        Jet_pt_float.clear();
        Jet_Eta_float.clear();
        Jet_Phi_float.clear();

        Track_pt_float.clear();
        Track_Eta_float.clear();
        Track_Phi_float.clear();

        // Delta R Clears
        chi_Muon_DeltaR.clear();
        chi_elec_DeltaR.clear();
        chi_tau_DeltaR.clear();

        scalar_Muon_DeltaR.clear();
        scalar_elec_DeltaR.clear();
        scalar_tau_DeltaR.clear();

        Psi_Muon_DeltaR.clear();
        Psi_elec_DeltaR.clear();
        Psi_tau_DeltaR.clear();

        Muon_Muon_Gen_Delta_R.clear();
        elec_elec_Gen_Delta_R.clear();



        tracks.clear();
        chi_Track_DeltaR.clear();
        scalar_Track_DeltaR.clear();
  }
  // Print and Write Files
  T->Print();
  T_Delta->Print();
  T_Det->Print();
  f->Write();
}
