// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// ATLAS Wee Wemu Wmumu analysis at Z TeV
  class WWbb : public Analysis {
  public:

    /// Default constructor
    WWbb()
    : Analysis("WWbb")
    {    }

    const MissingMomentum& MET_4v;
    DressedLepton lepton_m, lepton_p;

    // GenParticle  lepton_m, lepton_m, nu_m, nu_p            
    // double         m_ll, m_trans_llMET, m_Wm, m_Wp, MET
    // JetVector    lightjets, bjets_central, bjets_forward        
    // double         m_tm, m_tp
    // map<jet,GenParticle>    bflavours
    // Jet         bjet_p, bjet_m

    // Hist1D         cuts_WW, cuts_VBF, cuts_HH, cuts_Mass, cuts_BL (Mass_BL)

    void init() {
      FinalState fs;

      Cut etaRanges_EL = (Cuts::abseta < 1.37 || Cuts::absetaIn(1.52, 2.47)) && Cuts::pT > 20*GeV;
      Cut etaRanges_MU = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;

      MissingMomentum met(fs);
      addProjection(met, "MET");

      IdentifiedFinalState Photon(fs);
      Photon.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState bare_EL(fs);
      bare_EL.acceptIdPair(PID::ELECTRON);

      IdentifiedFinalState bare_MU(fs);
      bare_MU.acceptIdPair(PID::MUON);

      IdentifiedFinalState neutrinoFS(fs);
      neutrinoFS.acceptNeutrinos();
      addProjection(neutrinoFS, "Neutrinos");

      ////////////////////////////////////////////////////////
      // DRESSED LEPTONS
      //    3.arg: 0.1      = dR(lep,phot)
      //    4.arg: true     = do clustering
      //    7.arg: false    = ignore photons from hadron or tau
      //
      //////////////////////////////////////////////////////////
      DressedLeptons electronFS(Photon, bare_EL, 0.1, etaRanges_EL);
      addProjection(electronFS, "ELECTRON_FS");

      DressedLeptons muonFS(Photon, bare_MU, 0.1, etaRanges_MU);
      addProjection(muonFS, "MUON_FS");

      VetoedFinalState jetinput;
      jetinput.addVetoOnThisFinalState(bare_MU);
      jetinput.addVetoOnThisFinalState(neutrinoFS);

      FastJets jetpro(jetinput, FastJets::ANTIKT, 0.4);
      addProjection(jetpro, "jet");


      initialize_Histos();
      
    }

    void initialize_Histos(){
      
      // put global stuff here
      // ...
      
      // initialize analysis dependent stuff
      initialize_Histos_WW();
      initialize_Histos_WBF();
      initialize_Histos_HH();
      initialize_Histos_BL();
      
    };

    void initialize_Histos_WW(){
      
    }
    
    void analyze_WW(){
      
    }

    void initialize_Histos_WBF(){
      
    }
    
    void analyze_WBF(){
      
    }

    void initialize_Histos_HH(){
      
    }
    
    void analyze_HH(){
      
    }

    void initialize_Histos_BL(){
      
    }
    
    void analyze_BL(){
      
    }

    /// Do the analysis
    void analyze(const Event& e) {

      analyze_WW();
      
      analyze_WBF();
      
      analyze_HH();
      
      analyze_BL();
      
    }


    /// Finalize
    void finalize() {
    }


  private:
    
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(WWbb);

}
