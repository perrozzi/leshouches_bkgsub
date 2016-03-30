// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// ATLAS Wee Wemu Wmumu analysis at Z TeV
  class WWbb : public Analysis {
  private:
    double lepton_etamax,lepton_ptmin;
    double jet_etamax,jet_ptmin;
    double lepton_jet_isolation_dR, lepton_iso_dR, lepton_iso_frac;
    double bhad_ptmin;
    FourMomentum MET_4v;
    Particle lepton_m, lepton_p;//, nu_m, nu_p;
    Particles nu_m, nu_p;
    double m_ll, m_trans_llMET, m_Wm, m_Wp, MET;
    Jets alljets, lightjets, bjets_central, bjets_forward;
    double m_tm, m_tp;
    double m_bl_p, m_bl_m;
    // map<jet,GenParticle>    bflavours
    Jets bjet_p, bjet_m;


    double massJJ_min_WBF, deltayJJ_min_WBF;
    double m_trans_llMET_min_WBF, m_ll_min_WBF;
    double ptlep1_min_WBF,ptlep2_min_WBF,MET_min_WBF;


    double METrel_min_hh, mass_ll_min_hh, bjets_central_min_hh, m_trans_llMET_min_hh, m_trans_llMET_max_hh, mBB_min_hh, mBB_max_hh;
    double m_bl_min_hh,m_bl_max_hh;

    double m_ll_min_ww, MET_min_ww, alljets_max_ww;

    Histo1DPtr h_nbhadrons;
    Histo1DPtr h_pt_Wm, h_pt_Wp;
    Histo1DPtr h_rapidity_Wm, h_rapidity_Wp;
    Histo1DPtr h_mass_Wm, h_mass_Wp;
    Histo1DPtr h_massZoom_Wm, h_massZoom_Wp;
    Histo1DPtr h_mass_tm, h_mass_tp;
    Histo1DPtr h_massZoom_tm, h_massZoom_tp;
    Histo1DPtr h_mass_bl_m, h_mass_bl_p;
    Histo1DPtr h_massZoom_bl_m, h_massZoom_bl_p;

    Histo1DPtr presel_njets_after, presel_nbjets_central_after, presel_nbjets_forward_after;
    Histo1DPtr presel_jets_pt_after, presel_bjets_pt_central_after, presel_bjets_pt_forward_after;
    Histo1DPtr presel_jets_eta_after, presel_bjets_eta_central_after, presel_bjets_eta_forward_after;
    Histo1DPtr presel_cuts;
    Histo1DPtr WBF_cuts, WBF_njets_before, WBF_njets_after;
    Histo1DPtr hh_cuts, hh_njets_before, hh_njets_after;
    Histo2DPtr hh_njets_nbjets_before, hh_njets_nbjets_after;
    Histo1DPtr ww_cuts, ww_njets_before;

  public:

    /// Default constructor
    WWbb() : Analysis("WWbb"),
    lepton_etamax(2.4), lepton_ptmin(25.*GeV),
    jet_etamax(4.5), jet_ptmin(25.*GeV),
    lepton_jet_isolation_dR(0.4),
    lepton_iso_dR(0.4), lepton_iso_frac(0.1),
    bhad_ptmin(5.*GeV), 
    METrel_min_hh(25*GeV), mass_ll_min_hh(10*GeV), bjets_central_min_hh(2),
    m_trans_llMET_min_hh(100*GeV), m_trans_llMET_max_hh(150*GeV), mBB_min_hh(100*GeV), mBB_max_hh(150*GeV),
    m_bl_min_hh(100*GeV), m_bl_max_hh(180*GeV),
    m_ll_min_ww(10*GeV), MET_min_ww(15*GeV), alljets_max_ww(0)	 

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
      VetoedFinalState jetinput;
      // jetinput.addVetoOnThisFinalState(muon_bare);
      jetinput.addVetoOnThisFinalState(neutrinos);

      FastJets jetpro(jetinput, FastJets::ANTIKT, 0.4);
      addProjection(jetpro, "jet");

      addProjection(jetinput,"vfs");
      
      ////////////////////////////////////////////////////////
      // MET
      ////////////////////////////////////////////////////////

      MissingMomentum met(fs);
      addProjection(met, "MET");

      ////////////////////////////////////////////////////////
      // Unstable Particles
      ////////////////////////////////////////////////////////
      UnstableFinalState ufs(Cuts::pT > bhad_ptmin);
      addProjection(ufs, "UFS");

      initialize_Histos();
      
    }

    ////////////////////////////////////////////////////////
    // initialize_Histos
    ////////////////////////////////////////////////////////
    void initialize_Histos(){
      presel_cuts     = bookHisto1D("presel_cuts",20,-0.5,19.5);

      h_nbhadrons = bookHisto1D("h_nbhadrons",10,0,10);
      h_pt_Wm = bookHisto1D("h_pt_Wm",100,0*GeV,100*GeV);
      h_pt_Wp = bookHisto1D("h_pt_Wp",100,0*GeV,100*GeV);
      h_rapidity_Wm = bookHisto1D("h_rapidity_Wm",50,-5,5);
      h_rapidity_Wp = bookHisto1D("h_rapidity_Wp",50,-5,5);
      h_mass_Wm = bookHisto1D("h_mass_Wm",24,0*GeV,120*GeV);
      h_mass_Wp = bookHisto1D("h_mass_Wp",24,0*GeV,120*GeV);
      h_massZoom_Wm = bookHisto1D("h_massZoom_Wm",40,60*GeV,100*GeV);
      h_massZoom_Wp = bookHisto1D("h_massZoom_Wp",40,60*GeV,100*GeV);

      h_mass_tm = bookHisto1D("h_mass_tm",40,0*GeV,200*GeV);
      h_mass_tp = bookHisto1D("h_mass_tp",40,0*GeV,200*GeV);
      h_massZoom_tm = bookHisto1D("h_massZoom_tm",40,150*GeV,190*GeV);
      h_massZoom_tp = bookHisto1D("h_massZoom_tp",40,150*GeV,190*GeV);

      h_mass_bl_m = bookHisto1D("h_mass_bl_m",40,0*GeV,200*GeV);
      h_mass_bl_p = bookHisto1D("h_mass_bl_p",40,0*GeV,200*GeV);
      h_massZoom_bl_m = bookHisto1D("h_massZoom_bl_m",30,60*GeV,190*GeV);
      h_massZoom_bl_p = bookHisto1D("h_massZoom_bl_p",30,60*GeV,190*GeV);

      presel_njets_after = bookHisto1D("presel_njets_after",10,-0.5,9.5);
      presel_nbjets_central_after = bookHisto1D("presel_nbjets_central_after",10,-0.5,9.5);
      presel_nbjets_forward_after = bookHisto1D("presel_nbjets_forward_after",10,-0.5,9.5);
      presel_jets_pt_after = bookHisto1D("presel_jets_pt_after",100,0,100);
      presel_bjets_pt_central_after = bookHisto1D("presel_bjets_pt_central_after",100,0*GeV,100*GeV);
      presel_bjets_pt_forward_after = bookHisto1D("presel_bjets_pt_forward_after",100,0*GeV,100*GeV);
      presel_jets_eta_after = bookHisto1D("presel_jets_eta_after",50,-5,5);
      presel_bjets_eta_central_after = bookHisto1D("presel_bjets_eta_central_after",50,-5,5);
      presel_bjets_eta_forward_after = bookHisto1D("presel_bjets_eta_forward_after",50,-5,5);


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
      
      const double weight = event.weight();
      presel_cuts->fill(0,weight);
      WBF_cuts->fill(0,weight);
      hh_cuts->fill(0,weight);
      ww_cuts->fill(0,weight);

      alljets.clear();
      lightjets.clear();
      bjets_central.clear();
      bjets_forward.clear();     
      nu_m.clear();
      nu_p.clear();
      bjet_m.clear();
      bjet_p.clear();

      ////////////////////////////////////////////////////////
      // Visible particles for the isolation
      //////////////////////////////////////////////////////// 
      Particles vfs_particles=applyProjection<VetoedFinalState>(event,"vfs").particles();

      ////////////////////////////////////////////////////////
      // unstable particles for the b-hadron identification
      ////////////////////////////////////////////////////////
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(event, "UFS");

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

      //cout << " ------------------------------------" << endl;
      //cout << " n_e = " << electron_dummy.size() << "    n_mu=  " << muon_dummy.size() << endl;
      if ((electron_dummy.size()==1 && muon_dummy.size()==1) /*&& (electron_dummy[0].charge()*muon_dummy[0].charge()==-1)*/) presel_cuts->fill(1,weight);
      if ((electron_dummy.size()==1 && muon_dummy.size()==1) && (electron_dummy[0].charge()*muon_dummy[0].charge()==-1)) presel_cuts->fill(2,weight);



      ////////////////////////////////////////////////////////
      // isolated leptons
      ////////////////////////////////////////////////////////
      vector<DressedLepton> muon_isolated, electron_isolated;
      foreach (DressedLepton& l1, muon_dummy) {
        double ptcone=0;
        foreach (Particle& p, vfs_particles) {
          double deltaRlpart = deltaR(l1.momentum(), p.momentum(), RAPIDITY);
          if (deltaRlpart<lepton_iso_dR) ptcone=ptcone+p.momentum().pT();
        }
        if (ptcone<(1.+lepton_iso_frac)*l1.pT()) 
        muon_isolated.push_back(l1);//keep mu with highest pT
      }
      foreach (DressedLepton& l1, electron_dummy) {
        double ptcone=0;         
        foreach (Particle& p, vfs_particles) { 
          double deltaRlpart = deltaR(l1.momentum(), p.momentum(), RAPIDITY);
          if (deltaRlpart<lepton_iso_dR) ptcone=ptcone+p.momentum().pT();
        }
        if (ptcone<(1.+lepton_iso_frac)*l1.pT())   
        electron_isolated.push_back(l1);//keep e with highest pT
      }
      
      /////////////////////////////////////////////////////////////////////////
      // JETS
      /////////////////////////////////////////////////////////////////////////
      foreach (const Jet& j, applyProjection<FastJets>(event, "jet").jetsByPt(jet_ptmin)) {
        if (j.absrap() > jet_etamax ) continue;
        bool deltaRcontrol = true;
        foreach (DressedLepton& fl,muon_isolated) {
          double deltaRjets = deltaR(fl.constituentLepton().momentum(), j.momentum(), RAPIDITY);
          if (deltaRjets <= lepton_jet_isolation_dR) deltaRcontrol = false; //false if at least one electron is in the overlap region
        }
        foreach (DressedLepton& fl,electron_isolated) {
          double deltaRjets = deltaR(fl.constituentLepton().momentum(), j.momentum(), RAPIDITY);
          if (deltaRjets <= lepton_jet_isolation_dR) deltaRcontrol = false; //false if at least one electron is in the overlap region
        }
        if (deltaRcontrol) alljets.push_back(j);
      }

      /////////////////////////////////////////////////////////////////////////
      // B-Hadrons
      /////////////////////////////////////////////////////////////////////////
      Particles BHadrons;
      // Loop over the unstable particles
      foreach (const Particle& p, ufs.particles()) {
        const PdgId pid = p.pid();

        // Look for particles with a bottom quark
        if (PID::hasBottom(pid)) {

          bool good_B = false;
          const GenParticle* pgen = p.genParticle();
          const GenVertex* vgen = pgen -> end_vertex();

          // Loop over the decay products of each unstable particle.
          // Look for a couple of B hadrons.
          for (GenVertex::particles_out_const_iterator it = vgen->particles_out_const_begin(); it !=  vgen->particles_out_const_end(); ++it) {
            // If the particle produced has a bottom quark do not count it and go to the next loop cycle.
            if (!( PID::hasBottom( (*it)->pdg_id() ) ) ) {
              good_B = true;
              continue;
            } else {
              good_B = false;
              break;
            }
          }
          if (good_B ) BHadrons.push_back( p );
        }
        else continue;
      }
      h_nbhadrons->fill(BHadrons.size(),weight);
      //cout << BHadrons.size()<<endl;

      /////////////////////////////////////////////////////////////////////////
      // Jet B-Label
      /////////////////////////////////////////////////////////////////////////     
      vector<int> bjets_central_index, bjets_forward_index;
      foreach (const Jet& j, alljets) {
        double dRmin=1000;
        int bindex=0;
        int bfound=-1;
        foreach (const Particle& b, BHadrons){
          bool notusedbhad=true;
          if (find(bjets_central_index.begin(),bjets_central_index.end(),bindex)!=bjets_central_index.end())notusedbhad=false;
          if (find(bjets_forward_index.begin(),bjets_forward_index.end(),bindex)!=bjets_forward_index.end())notusedbhad=false;
          if (notusedbhad){
            double dR=deltaR(b.momentum(), j.momentum(), RAPIDITY);;
            if (dR<dRmin){
              dRmin=dR;
              bfound=bindex;
            }
          }
          bindex++;
        }
        if (dRmin<0.4 && j.absrap()<2.4){bjets_central.push_back(j); bjets_central_index.push_back(bfound);}
        if (dRmin<0.4 && j.absrap()>=2.4){bjets_forward.push_back(j); bjets_forward_index.push_back(bfound);}
        if (dRmin>=0.4)lightjets.push_back(j);
      }

      /////////////////////////////////////////////////////////////////////////
      // Event Preselection
      /////////////////////////////////////////////////////////////////////////

      //cout <<"  n_e =  "<< electron_isolated.size() << "    n_mu=  " << muon_isolated.size() << endl;

      if (electron_isolated.size()!=1 || muon_isolated.size()!=1) return;
      presel_cuts->fill(3,weight);
      if (electron_isolated[0].charge()*muon_isolated[0].charge()!=-1) return;
      presel_cuts->fill(4,weight);

      if (electron_isolated[0].charge()>0){
        presel_cuts->fill(5,weight);
        lepton_m = muon_isolated[0];
        lepton_p = electron_isolated[0];
      }
      else{
        presel_cuts->fill(6,weight);
        lepton_p = muon_isolated[0];
        lepton_m = electron_isolated[0];
      }

      m_ll = (lepton_m.momentum() + lepton_p.momentum()).mass();
      FourMomentum METll_4v = lepton_m.momentum() + lepton_p.momentum() + MET_4v;
      // CHECK THAT LONGITUDINAL MOMENTUM IS REMOVED
      m_trans_llMET = sqrt(METll_4v.Et2() - METll_4v.pT2());


      Particles Neutrinos= applyProjection< LeadingParticlesFinalState>(event,"neutrinos").particlesByPt();
      foreach (const Particle& inu, Neutrinos){
        if (inu.pid()+lepton_p.pid()==1) {
          nu_p.push_back(inu); //W+;
          break;
        }
      }
      foreach (const Particle& inu, Neutrinos){
        if (inu.pid()+lepton_m.pid()==-1) {
          nu_m.push_back(inu); //W-;
          break;
        }
      }

      for (unsigned int ibjet=0; ibjet<bjets_central.size();ibjet++){
        if (BHadrons[bjets_central_index[ibjet]].charge()>0) {
          bjet_m.push_back(bjets_central[ibjet]);
          break;
        }
      }

      for (unsigned int ibjet=0; ibjet<bjets_central.size();ibjet++){
        if (BHadrons[bjets_central_index[ibjet]].charge()<0) {
          bjet_p.push_back(bjets_central[ibjet]);
          break;
        }
      }

      if(bjets_central.size()>1 && bjet_m.size()+bjet_p.size()==1) {
        if(bjet_m.size()==0 && bjets_central[0].pT()==bjet_p[0].pT()) bjet_m.push_back(bjets_central[1]);
        else if (bjet_m.size()==0) bjet_m.push_back(bjets_central[0]);
        if(bjet_p.size()==0 && bjets_central[0].pT()==bjet_m[0].pT()) bjet_p.push_back(bjets_central[1]);
        else if (bjet_p.size()==0) bjet_p.push_back(bjets_central[0]);	         
      }
      

      m_bl_m=0;
      m_bl_p=0;

      if(nu_m.size()>0) {
        presel_cuts->fill(7,weight);
        m_Wm = (lepton_m.momentum() + nu_m[0].momentum()).mass();
        if(bjet_m.size()>0 && bjets_central.size()>1){
          presel_cuts->fill(9,weight);
          m_tm = (lepton_m.momentum() + bjet_m[0].momentum() + nu_m[0].momentum()).mass();
          m_bl_m = (lepton_m.momentum() + bjet_m[0].momentum()).mass();

          h_mass_bl_m->fill(m_bl_m,weight);
          h_massZoom_bl_m->fill(m_bl_m,weight);
          h_mass_tm->fill(m_tm,weight);
          h_massZoom_tm->fill(m_tm,weight);
        }
        h_pt_Wm->fill((lepton_m.momentum() + nu_m[0].momentum()).pt(),weight);
        h_rapidity_Wm->fill((lepton_m.momentum() + nu_m[0].momentum()).rapidity(),weight);
        h_mass_Wm->fill(m_Wm,weight);
        h_massZoom_Wm->fill(m_Wm,weight);
      }

      if(nu_p.size()>0) {
        presel_cuts->fill(8,weight);
        m_Wp = (lepton_p.momentum() + nu_p[0].momentum()).mass();
        if(bjet_p.size()>0 && bjets_central.size()>1){
          presel_cuts->fill(10,weight);
          m_tp = (lepton_p.momentum() + bjet_p[0].momentum() + nu_p[0].momentum()).mass();
          m_bl_p = (lepton_p.momentum() + bjet_p[0].momentum()).mass();

          h_mass_bl_p->fill(m_bl_p,weight);
          h_massZoom_bl_p->fill(m_bl_p,weight);
          h_mass_tp->fill(m_tp,weight);
          h_massZoom_tp->fill(m_tp,weight);
        }
        h_pt_Wp->fill((lepton_p.momentum() + nu_p[0].momentum()).pt(),weight);
        h_rapidity_Wp->fill((lepton_p.momentum() + nu_p[0].momentum()).rapidity(),weight);
        h_mass_Wp->fill(m_Wp,weight);
        h_massZoom_Wp->fill(m_Wp,weight);
      }

      presel_njets_after->fill(alljets.size(),weight);
      presel_nbjets_central_after->fill(bjets_central.size(),weight);
      presel_nbjets_forward_after->fill(bjets_forward.size(),weight);
      
      if(alljets.size()>0){
        presel_jets_pt_after->fill(alljets[0].pt(),weight);
        presel_jets_eta_after->fill(alljets[0].eta(),weight);
      }
      if(bjets_central.size()>0){
        presel_bjets_pt_central_after->fill(bjets_central[0].pt(),weight);
        presel_bjets_eta_central_after->fill(bjets_central[0].eta(),weight);
      }
      if(bjets_forward.size()>0){
        presel_bjets_pt_forward_after->fill(bjets_forward[0].pt(),weight);
        presel_bjets_eta_forward_after->fill(bjets_forward[0].eta(),weight);
      }

      ////////////////////////////////////////////////////////
      // RUN ANALYSES
      ////////////////////////////////////////////////////////
      WBF_cuts->fill(1,weight);
      hh_cuts->fill(1,weight);
      ww_cuts->fill(1,weight);

      analyze_WW(event);
      
      analyze_WBF(event);
      
      analyze_HH(event);
      
      analyze_BL(event);
      
      
    }


    ////////////////////////////////////////////////////////
    // WW
    ////////////////////////////////////////////////////////
    void initialize_Histos_WW(){
      ww_cuts     = bookHisto1D("ww_cuts", 20,-0.5,19.5);
      ww_njets_before = bookHisto1D("ww_njets_before",10,-0.5,9.5);     
    }
    
    void analyze_WW(const Event& event){
      const double weight = event.weight();
      ww_njets_before->fill(alljets.size(),weight);     
      ww_cuts->fill(2,weight);

      if(m_ll < m_ll_min_ww) return; // m_ll_min_ww(10*GeV)
      ww_cuts->fill(3,weight);
      
      if(MET < MET_min_ww) return; // MET_min_ww(15*GeV)
      ww_cuts->fill(4,weight);

      if(alljets.size() > alljets_max_ww) return; // alljets_max_ww(0)
      ww_cuts->fill(5,weight);

      // no MPT cut
      // no dPhi(MET,MPT) cut           
    }

    ////////////////////////////////////////////////////////
    // WBF
    ////////////////////////////////////////////////////////
    void initialize_Histos_WBF(){
      
      WBF_cuts     = bookHisto1D("WBF_cuts",5,-0.5,4.5);
      WBF_njets_before = bookHisto1D("WBF_njets_before",10,-0.5,9.5);
      WBF_njets_after  = bookHisto1D("WBF_njets_after",10,-0.5,9.5);

    }
    
    void analyze_WBF(const Event& event){

      const double weight = event.weight();
      WBF_njets_before->fill(alljets.size(),weight);

      // cut on 2 jets in opposite hemispheres with minimal mass and rap distance
      if (alljets.size()<2) return;
      double y0(alljets[0].momentum().eta());
      double y1(alljets[1].momentum().eta());
      double massJJ((alljets[0].momentum()+alljets[1].momentum()).mass());
      if (y0*y1>0. || fabs(y0-y1)<deltayJJ_min_WBF || massJJ<massJJ_min_WBF) return;
      WBF_cuts->fill(2,weight);

      // cuts on 2 lepton MET system
      if (m_trans_llMET<m_trans_llMET_min_WBF || m_ll<m_ll_min_WBF) return;
      double ptm(lepton_m.momentum().pT()),ptp(lepton_p.momentum().pT());
      double ptlep1(ptm>ptp?ptm:ptp), ptlep2(ptm>ptp?ptp:ptm);
      if (ptlep1<ptlep1_min_WBF || ptlep2<ptlep2_min_WBF || MET<MET_min_WBF) return;
      WBF_cuts->fill(3,weight);
      
      WBF_njets_after->fill(alljets.size(),weight);

      // veto if tag jets are central (i.e. tagged) bjets
      // if (bjets_central.find(alljets[0])!=bjet_central.end() ||
      // bjets_central.find(alljets[1])!=bjet_central.end()) return
      // WBF_cuts->fill(4,weight);
      
    }

    ////////////////////////////////////////////////////////
    // HH
    ////////////////////////////////////////////////////////
    void initialize_Histos_HH(){
      hh_cuts     = bookHisto1D("hh_cuts", 20,-0.5,19.5);
      hh_njets_before = bookHisto1D("hh_njets_before",10,-0.5,9.5);
      hh_njets_after = bookHisto1D("hh_njets_after",10,-0.5,9.5);
      hh_njets_nbjets_before     = bookHisto2D("hh_njets_nbjets_before",10,-0.5,9.5, 10,-0.5,9.5);
      hh_njets_nbjets_after     = bookHisto2D("hh_njets_nbjets_after",10,-0.5,9.5, 10,-0.5,9.5);      
    }
    
    void analyze_HH(const Event& event){
      const double weight = event.weight();
      hh_njets_before->fill(alljets.size(),weight);
      hh_njets_nbjets_before->fill(alljets.size(), bjets_central.size(), weight);

      // Objects used: MET, lepton_m, lepton_p, alljets, bjets_central, m_ll, m_trans_llMET

      // Missing ETrel
      //            *"mismom"   for "delta_phi" >= (0.5*pi)
      //            *"mismom.pT()*sin(delta_phi)"   for "delta_phi" < (0.5*pi)
      double METrel, delta_phi;
      METrel = 0, delta_phi = 0;
      vector<double> vL_MET_angle, vJet_MET_angle;
      vL_MET_angle.push_back(fabs(deltaPhi(lepton_m.momentum(), MET_4v)));
      vL_MET_angle.push_back(fabs(deltaPhi(lepton_p.momentum(), MET_4v)));
      foreach (double& lM, vL_MET_angle) if (lM > M_PI) lM = 2*M_PI - lM;
      std::sort(vL_MET_angle.begin(), vL_MET_angle.end());
      if (alljets.size() == 0) delta_phi = vL_MET_angle[0];
      if (alljets.size() > 0) {
        foreach (Jet& vj, alljets) {
          double jet_MET_angle = fabs(deltaPhi(vj.momentum(), MET_4v));
          if (jet_MET_angle > M_PI) jet_MET_angle = 2*M_PI - jet_MET_angle;
          vJet_MET_angle.push_back(jet_MET_angle);
        }
        std::sort(vJet_MET_angle.begin(), vJet_MET_angle.end());
        if (vL_MET_angle[0] <= vJet_MET_angle[0]) delta_phi = vL_MET_angle[0];
        if (vL_MET_angle[0] > vJet_MET_angle[0]) delta_phi = vJet_MET_angle[0];
      }  
      if (delta_phi >= (0.5*M_PI)) delta_phi = 0.5*M_PI;
      METrel = MET_4v.pT()*sin(delta_phi);  

      // Selection
      if(/*METrel*/MET <= METrel_min_hh) return; // METrel_min_hh = 25*GeV
      hh_cuts->fill(2,weight);

      if(m_ll <= mass_ll_min_hh) return; // mass_ll_min_hh = 10*GeV
      hh_cuts->fill(3,weight);

      if(bjets_central.size() < bjets_central_min_hh) return; // bjets_central_min_hh = 2
      hh_cuts->fill(4,weight);
      hh_njets_after->fill(alljets.size(),weight);
      hh_njets_nbjets_after->fill(alljets.size(), bjets_central.size(), weight);

      // 2 leading b-jets invariant mass
      double mBB = (bjets_central[0].momentum() + bjets_central[1].momentum()).mass();
      
      if(m_trans_llMET < m_trans_llMET_min_hh || m_trans_llMET > m_trans_llMET_max_hh) return; // m_trans_llMET_min_hh = 100*GeV; m_trans_llMET_max_hh = 150*GeV
      hh_cuts->fill(5,weight);

      if(mBB < mBB_min_hh || mBB > mBB_max_hh) return; // mBB_min_hh = 100*GeV; mBB_max_hh = 150*GeV
      hh_cuts->fill(6,weight);

      if (m_bl_m < m_bl_max_hh && m_bl_m > m_bl_min_hh) return; // m_bl_min_hh = 100*GeV ; m_bl_max_hh = 180*GeV
      if (m_bl_p < m_bl_max_hh && m_bl_p > m_bl_min_hh) return; // m_bl_min_hh = 100*GeV ; m_bl_max_hh = 180*GeV
      hh_cuts->fill(7,weight);


      
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
      
      // cout << "crossSection()/picobarn/sumOfWeights()= " << crossSection() << " " << picobarn << " " << sumOfWeights() << endl;
      scale(h_nbhadrons, crossSection()/picobarn/sumOfWeights());
      scale(h_pt_Wm, crossSection()/picobarn/sumOfWeights());
      scale(h_pt_Wp, crossSection()/picobarn/sumOfWeights());
      scale(h_rapidity_Wm, crossSection()/picobarn/sumOfWeights());
      scale(h_rapidity_Wp, crossSection()/picobarn/sumOfWeights());
      scale(h_mass_Wm, crossSection()/picobarn/sumOfWeights());
      scale(h_mass_Wp, crossSection()/picobarn/sumOfWeights());
      scale(h_massZoom_Wm, crossSection()/picobarn/sumOfWeights());
      scale(h_massZoom_Wp, crossSection()/picobarn/sumOfWeights());
      scale(h_mass_tm, crossSection()/picobarn/sumOfWeights());
      scale(h_mass_tp, crossSection()/picobarn/sumOfWeights());
      scale(h_massZoom_tm, crossSection()/picobarn/sumOfWeights());
      scale(h_massZoom_tp, crossSection()/picobarn/sumOfWeights());
      scale(h_mass_bl_m, crossSection()/picobarn/sumOfWeights());
      scale(h_mass_bl_p, crossSection()/picobarn/sumOfWeights());
      scale(h_massZoom_bl_m, crossSection()/picobarn/sumOfWeights());
      scale(h_massZoom_bl_p, crossSection()/picobarn/sumOfWeights());
      scale(presel_njets_after, crossSection()/picobarn/sumOfWeights());
      scale(presel_nbjets_central_after, crossSection()/picobarn/sumOfWeights());
      scale(presel_nbjets_forward_after, crossSection()/picobarn/sumOfWeights());
      scale(presel_jets_pt_after, crossSection()/picobarn/sumOfWeights());
      scale(presel_bjets_pt_central_after, crossSection()/picobarn/sumOfWeights());
      scale(presel_bjets_pt_forward_after, crossSection()/picobarn/sumOfWeights());
      scale(presel_jets_eta_after, crossSection()/picobarn/sumOfWeights());
      scale(presel_bjets_eta_central_after, crossSection()/picobarn/sumOfWeights());
      scale(presel_bjets_eta_forward_after, crossSection()/picobarn/sumOfWeights());
      scale(WBF_njets_before, crossSection()/picobarn/sumOfWeights());
      scale(WBF_njets_after, crossSection()/picobarn/sumOfWeights());
      scale(hh_njets_before, crossSection()/picobarn/sumOfWeights());
      scale(hh_njets_after, crossSection()/picobarn/sumOfWeights());
      // scale(hh_njets_nbjets_before, crossSection()/picobarn/sumOfWeights());
      // scale(hh_njets_nbjets_after, crossSection()/picobarn/sumOfWeights());
      
      scale(presel_cuts, crossSection()/picobarn/sumOfWeights());
      scale(WBF_cuts, crossSection()/picobarn/sumOfWeights());
      scale(hh_cuts, crossSection()/picobarn/sumOfWeights());
      
    }


  private:
    
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(WWbb);

}
