// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {


  /// ATLAS Wee Wemu Wmumu analysis at Z TeV
  class WWbb : public Analysis {
  private:
    double lepton_etamax,lepton_ptmin;
    double jet_etamax,jet_ptmin;
    double lepton_jet_isolation_dR, lepton_iso_dR, lepton_iso_frac;
    FourMomentum MET_4v;
    Particle lepton_m, lepton_p, nu_m, nu_p;
    double m_ll, m_trans_llMET, m_Wm, m_Wp, MET;
    Jets alljets, lightjets, bjets_central, bjets_forward;
    double m_tm, m_tp;
    // map<jet,GenParticle>    bflavours
    Jet bjet_p, bjet_m;
    VetoedFinalState jetinput;

    double massJJ_min_WBF, deltayJJ_min_WBF;
    double m_trans_llMET_min_WBF, m_ll_min_WBF;
    double ptlep1_min_WBF,ptlep2_min_WBF,MET_min_WBF;

    Histo1DPtr cuts_WBF, njets_before_WBF, njets_after_WBF;

  public:

    /// Default constructor
    WWbb() : Analysis("WWbb"),
    lepton_etamax(2.4), lepton_ptmin(25.*GeV),
    jet_etamax(4.5), jet_ptmin(25.*GeV),
    lepton_jet_isolation_dR(0.4),
    lepton_iso_dR(0.4), lepton_iso_frac(0.1)
    {}


    void init() {
      
      FinalState fs;
      
      ////////////////////////////////////////////////////////
      // NEUTRINOS
      ////////////////////////////////////////////////////////

      LeadingParticlesFinalState neutrinos(FinalState(-50, 50, 0*GeV));
      neutrinos.addParticleIdPair(12);
      neutrinos.addParticleIdPair(14);
      addProjection(neutrinos, "neutrinos");

      ////////////////////////////////////////////////////////
      // BARE LEPTONS
      ////////////////////////////////////////////////////////
      
      LeadingParticlesFinalState muon_bare(FinalState(-2.6, 2.6, 0*GeV));
      muon_bare.addParticleIdPair(PID::MUON);
      addProjection(muon_bare, "muons");

      LeadingParticlesFinalState electron_bare(FinalState(-2.6, 2.6, 0*GeV));
      electron_bare.addParticleIdPair(PID::ELECTRON);
      addProjection(electron_bare, "electrons");

      ////////////////////////////////////////////////////////
      // PHOTONS
      ////////////////////////////////////////////////////////

      IdentifiedFinalState Photon(fs);
      Photon.acceptIdPair(PID::PHOTON);
      
      ////////////////////////////////////////////////////////
      // DRESSED LEPTONS
      //    3.arg: 0.1      = dR(lep,phot)
      //    4.arg: true     = do clustering
      //    7.arg: false    = ignore photons from hadron or tau
      //
      //////////////////////////////////////////////////////////
      
      Cut etaRanges_leptons = Cuts::abseta < lepton_etamax && Cuts::pT > lepton_ptmin;
      
      DressedLeptons muon_dressed(Photon, muon_bare, 0.1, etaRanges_leptons);
      addProjection(muon_dressed, "muon_dressed");
      
      DressedLeptons electron_dressed(Photon, electron_bare, 0.1, etaRanges_leptons);
      addProjection(electron_dressed, "electron_dressed");


      ////////////////////////////////////////////////////////
      // JETS
      ////////////////////////////////////////////////////////

      // jetinput.addVetoOnThisFinalState(muon_bare);
      jetinput.addVetoOnThisFinalState(neutrinos);

      FastJets jetpro(jetinput, FastJets::ANTIKT, 0.4);
      addProjection(jetpro, "jet");
      
      ////////////////////////////////////////////////////////
      // MET
      ////////////////////////////////////////////////////////

      MissingMomentum met(fs);
      addProjection(met, "MET");


      initialize_Histos();
      
    }

    ////////////////////////////////////////////////////////
    // initialize_Histos
    ////////////////////////////////////////////////////////
    void initialize_Histos(){
      
      // put global stuff here
      // ...
      
      // initialize analysis dependent stuff
      initialize_Histos_WW();
      initialize_Histos_WBF();
      initialize_Histos_HH();
      initialize_Histos_BL();
      
    };

    /// Do the analysis
    void analyze(const Event& event) {

      ////////////////////////////////////////////////////////
      // MET
      ////////////////////////////////////////////////////////
      const MissingMomentum& MET_tmp = applyProjection<MissingMomentum>(event, "MET");
      MET_4v = -MET_tmp.visibleMomentum();
      MET = MET_4v.pT();

      ////////////////////////////////////////////////////////
      // leptons
      ////////////////////////////////////////////////////////
      const  vector<DressedLepton>& muon_dressed = applyProjection<DressedLeptons>(event, "muon_dressed").dressedLeptons();
      vector<DressedLepton> muon_dummy;
      muon_dummy.insert(muon_dummy.end(), muon_dressed.begin(), muon_dressed.end());
      
      const  vector<DressedLepton>& electron_dressed = applyProjection<DressedLeptons>(event, "electron_dressed").dressedLeptons();
      vector<DressedLepton> electron_dummy;
      electron_dummy.insert(electron_dummy.end(), electron_dressed.begin(), electron_dressed.end());


      ////////////////////////////////////////////////////////
      // isolated leptons
      ////////////////////////////////////////////////////////
      vector<DressedLepton> muon_isolated, electron_isolated;
      foreach (DressedLepton& l1, muon_dummy) {
        // double ptcone=0;
        // foreach (Particle& jet, jetpro) { // Is this the way to loop on the particles in the jet input? using the jetinput because the neutrino are already excluded. 
          // double deltaRlpart = deltaR(l1.momentum(), jet.momentum(), RAPIDITY);
          // if (deltaRlpart<lepton_iso_dR) ptcone=ptcone+jet.momentum().pT();
        // }
        // if (ptcone<(1.+lepton_iso_frac)*l1.pT()) 
          muon_isolated.push_back(l1);//keep e with highest pT
      }
      foreach (DressedLepton& l1, electron_dummy) {
        electron_isolated.push_back(l1);//keep e with highest pT
      }
      
      /////////////////////////////////////////////////////////////////////////
      // JETS
      /////////////////////////////////////////////////////////////////////////
      Jets alljets, vetojets;
      foreach (const Jet& j, applyProjection<FastJets>(event, "jet").jetsByPt(jet_ptmin)) {
        if (j.absrap() > jet_etamax ) continue;
        alljets.push_back(j);
        bool deltaRcontrol = true;
        foreach (DressedLepton& fl,muon_isolated) {
          double deltaRjets = deltaR(fl.constituentLepton().momentum(), j.momentum(), RAPIDITY);
          if (deltaRjets <= lepton_jet_isolation_dR) deltaRcontrol = false; //false if at least one electron is in the overlap region
        }
        foreach (DressedLepton& fl,electron_isolated) {
          double deltaRjets = deltaR(fl.constituentLepton().momentum(), j.momentum(), RAPIDITY);
          if (deltaRjets <= lepton_jet_isolation_dR) deltaRcontrol = false; //false if at least one electron is in the overlap region
        }
        if (deltaRcontrol) vetojets.push_back(j);
      }
      
      
      // DUMMY STUFF
      if(muon_isolated.size()==0 || electron_isolated.size()==0 ) return;
      lepton_m = muon_isolated[0];
      lepton_p = electron_isolated[0];
      
      ////////////////////////////////////////////////////////
      // RUN ANALYSES
      ////////////////////////////////////////////////////////

      analyze_WW(event);
      
      analyze_WBF(event);
      
      analyze_HH(event);
      
      analyze_BL(event);
      
      
    }


    ////////////////////////////////////////////////////////
    // WW
    ////////////////////////////////////////////////////////
    void initialize_Histos_WW(){
      
    }
    
    void analyze_WW(const Event& event){
      
    }

    ////////////////////////////////////////////////////////
    // WBF
    ////////////////////////////////////////////////////////
    void initialize_Histos_WBF(){
      
      cuts_WBF     = bookHisto1D("cuts_WBF",5,-0.5,4.5);
      njets_before_WBF = bookHisto1D("njets_before_WBF",10,-0.5,9.5);
      njets_after_WBF  = bookHisto1D("njets_after_WBF",10,-0.5,9.5);

    }
    
    void analyze_WBF(const Event& event){

      const double weight = event.weight();
      njets_before_WBF->fill(alljets.size(),weight);

      // cut on 2 jets in opposite hemispheres with minimal mass and rap distance
      if (alljets.size()<2) return;
      double y0(alljets[0].momentum().eta());
      double y1(alljets[1].momentum().eta());
      double massJJ((alljets[0].momentum()+alljets[1].momentum()).mass());
      if (y0*y1>0. || fabs(y0-y1)<deltayJJ_min_WBF || massJJ<massJJ_min_WBF) return;
      cuts_WBF->fill(2,weight);

      // cuts on 2 lepton MET system
      if (m_trans_llMET<m_trans_llMET_min_WBF || m_ll<m_ll_min_WBF) return;
      double ptm(lepton_m.momentum().pT()),ptp(lepton_p.momentum().pT());
      double ptlep1(ptm>ptp?ptm:ptp), ptlep2(ptm>ptp?ptp:ptm);
      if (ptlep1<ptlep1_min_WBF || ptlep2<ptlep2_min_WBF || MET<MET_min_WBF) return;
      cuts_WBF->fill(3,weight);
      
      njets_after_WBF->fill(alljets.size(),weight);

      // veto if tag jets are central (i.e. tagged) bjets
      // if (bjets_central.find(alljets[0])!=bjet_central.end() ||
      // bjets_central.find(alljets[1])!=bjet_central.end()) return
      // cuts_WBF->fill(4,weight);
      
    }

    ////////////////////////////////////////////////////////
    // HH
    ////////////////////////////////////////////////////////
    void initialize_Histos_HH(){
      
    }
    
    void analyze_HH(const Event& event){
      
    }

    ////////////////////////////////////////////////////////
    // BL
    ////////////////////////////////////////////////////////
    void initialize_Histos_BL(){
      
    }
    
    void analyze_BL(const Event& event){
      
    }

    /// Finalize
    void finalize() {
    }


  private:
    
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(WWbb);

}
