
//IMPORTANT: Add .cc extension: WWGen.cc


#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "TFile.h"
#include "TH1.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <iostream>
#include <fstream>

#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

#include "Rivet/Projections/DressedLeptons.hh"

using namespace std;

namespace Rivet {


  

  class WWGen_cut_dl : public Analysis {
  public:

    /// Constructor
    WWGen_cut_dl()
      : Analysis("WWGen_cut_dl")
    { ce=0;
      cmu=0;  
      numev=0;
      numlep=0;
      numneu=0;
      numwm=0;
      numwp=0;
    }

  public:
       
    /// Book histograms and initialise projections before the run
    void init() {


      
      
      /// @todo Initialise and register projections here
      addProjection(VisibleFinalState(-50,50),"vfs");

      LeadingParticlesFinalState leadingLeptons(FinalState(-2.50, 2.50, 0*GeV));
      leadingLeptons.addParticleIdPair(PID::ELECTRON);
      leadingLeptons.addParticleIdPair(PID::MUON);
      //      leadingLeptons.acceptIdPair(15);
      addProjection(leadingLeptons, "leptons");

      LeadingParticlesFinalState leadingNeutrinos(FinalState(-5, 5, 0*GeV));
      leadingNeutrinos.addParticleIdPair(12);
      leadingNeutrinos.addParticleIdPair(14);
      //      leadingNeutrinos.acceptIdPair(16);
      addProjection(leadingNeutrinos, "neutrinos");

      addProjection(FastJets(FinalState(-6,6), FastJets::ANTIKT, 0.5), "jets");

      IdentifiedFinalState photon(FinalState(-5, 5, 0*GeV));
      photon.acceptId(22);
      addProjection(photon, "Photon");

      DressedLeptons dl ( photon,  leadingLeptons, 0.1,true,-5,5,0*GeV,false);
      addProjection(dl, "DressedLeptons");




      //Higgs
      _h_mh  = bookHisto1D("mh", 130, 50, 180);
      _h_mh_0j  = bookHisto1D("mh_0j", 130, 50, 180);
      _h_mh_1j  = bookHisto1D("mh_1j", 130, 50, 180);
      _h_mh_2j  = bookHisto1D("mh_2j", 130, 50, 180);
      _h_mth  = bookHisto1D("mth", 25, 0, 250);
      _h_mth_0j  = bookHisto1D("mth_0j", 25, 0, 250);
      _h_mth_1j  = bookHisto1D("mth_1j", 25, 0, 250);
      _h_mth_2j  = bookHisto1D("mth_2j", 25, 0, 250);  
      _h_pth  = bookHisto1D("pth", 25, 0, 250); //pt higgs
      _h_etah  = bookHisto1D("etah", 25, 0, 4); //eta higgs 
      
      //leptons
      _h_mll  = bookHisto1D("mll", 25, 0, 250);
      _h_mll_0j  = bookHisto1D("mll_0j", 25, 0, 250);
      _h_mll_1j  = bookHisto1D("mll_1j", 25, 0, 250);
      _h_mll_2j  = bookHisto1D("mll_2j", 25, 0, 250);      
      _h_pt_1l  = bookHisto1D("pt_1l", 25, 0, 250); //pt leading lepton
      _h_pt_2l  = bookHisto1D("pt_2l", 25, 0, 250); //pt sub-leading lepton
      _h_pt_ll  = bookHisto1D("pt_ll", 25, 0, 250); //pt leading + sub-leading leptons
      _h_eta_1l  = bookHisto1D("eta_1l", 20, -10, 10); //eta leading lepton
      _h_eta_2l  = bookHisto1D("eta_2l", 20, -10, 10); //eta sub-leading lepton
      _h_deta_l  = bookHisto1D("deta_l", 25, 0, 10); //Delta-eta
      _h_dphi_l  = bookHisto1D("dphi_l", 25, 0, 4); //Delta-phi leptons     

      //W+- mass
      _h_mwm = bookHisto1D("mwm", 160, 0, 160);
      _h_mwp = bookHisto1D("mwp", 160, 0, 160);

      _h_mwm_0j = bookHisto1D("mwm_0j", 160, 0, 160);
      _h_mwm_1j = bookHisto1D("mwm_1j", 160, 0, 160);
      _h_mwm_2j = bookHisto1D("mwm_2j", 160, 0, 160);

      _h_mwp_0j = bookHisto1D("mwp_0j", 160, 0, 160);
      _h_mwp_1j = bookHisto1D("mwp_1j", 160, 0, 160);
      _h_mwp_2j = bookHisto1D("mwp_2j", 160, 0, 160);

      //jets
      _h_mjj  = bookHisto1D("mjj", 25, 0, 250);//massa dei due jets     
      _h_pt_1j  = bookHisto1D("pt_1j", 25, 0, 250); //pt leading jet
      _h_pt_2j  = bookHisto1D("pt_2j", 25, 0, 250); //pt sub-leading jet
      _h_eta_1j  = bookHisto1D("eta_1j", 25, 0, 5); //eta leading lepton
      _h_eta_2j  = bookHisto1D("eta_2j", 25, 0, 5); //eta sub-leading lepton      
      _h_deta_j  = bookHisto1D("deta_j", 25, 0, 10); //Delta-eta jets
      _h_njet = bookHisto1D("njet", 4, -0.5, 3.5);
      _h_njet_onPeak = bookHisto1D("njet_onPeak", 4, -0.5, 3.5);
      _h_njet_offPeak = bookHisto1D("njet_offPeak", 4, -0.5, 3.5);
                 
      //Weight
      _h_weight= bookHisto1D("weight", 300, -300, 300);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double weight = event.weight();

      std::vector<Rivet::DressedLepton> dressed= applyProjection<  DressedLeptons >(event,"DressedLeptons").dressedLeptons();           
      //std::vector<Rivet::Particles> dressed= applyProjection<  DressedLeptons >(event,"DressedLeptons").dressedLeptons();  
      std::sort(dressed.begin(), dressed.end(), cmpMomByPt);
      
      Particles cand_leptons;
      foreach (const Particle & particle, dressed){//ByPt
	//    foreach ( const Particle & particle, applyProjection< LeadingParticlesFinalState>(event,"leptons").particlesByPt() ){//ByPt
	if (cand_leptons.size() == 0 && particle.pt() > 25*GeV && fabs(particle.eta()) < 2.4){
          cand_leptons.push_back(particle);}	 
        else if (cand_leptons.size() == 1 && particle.pt() > 10*GeV && fabs(particle.eta()) < 2.4){  
          cand_leptons.push_back(particle);	 
	}
      }
      
      if (cand_leptons.size() <2){ 
	//     	std::cout << "not enough leptons: STRANGE!!" << std::endl;
	vetoEvent;
      }  
	
      Particles vfs_particles =	applyProjection<VisibleFinalState>(event, "vfs").particles();
      FourMomentum pTmiss;
      foreach ( const Particle & p, vfs_particles ) {
	pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();
      
      Particles cand_neutrinos;
      foreach ( const Particle & neutrino, applyProjection< LeadingParticlesFinalState>(event,"neutrinos").particlesByPt() ){//ByPt
	if (cand_neutrinos.size() < 2)
          cand_neutrinos.push_back(neutrino);
      }

      if (cand_neutrinos.size()<2  ){ 
	//	std::cout << "not enough neutrinos: STRANGE!!" << std::endl;
       	vetoEvent;          
      }  
        

      Jets cand_jets;
      foreach (const Jet& jet,
	       applyProjection<FastJets>(event, "jets").jetsByPt(30.0*GeV) ) {
	if ( fabs( jet.eta() ) < 4.9 ) {
	  bool isolated = true;
	  foreach (const Particle& neutrino, cand_neutrinos){
	    if (deltaR(neutrino, jet) < 0.5){
	      isolated = false;
	      break;
	    }  
	  }
	  if (isolated) {
	    foreach (const Particle& lepton, cand_leptons){
	      if (deltaR(lepton, jet) < 0.5){
		isolated = false;
		break;
	      }
	    }
	  }
	  if (isolated)
	    cand_jets.push_back(jet);
	}
      } 

      unsigned int njet = 0;
      foreach (const Jet& jet,  cand_jets){
	if (jet.momentum().pT()>30*GeV)    
	  njet += 1;
      }

      //W+- 
      FourMomentum  wm;
      FourMomentum  wp;
     
      if(   ((cand_leptons[0].charge() )* (cand_leptons[1].charge()))==1   ){
	vetoEvent;

      }

      //per il leptone 1
      //      for (size_t  j=0; j< 2 ;++j){
      for (size_t  i=0; i< 2 ;++i){
	if   (  (cand_leptons[0].pid()+cand_neutrinos[i].pid())==-1){//caso W-
	  wm=cand_leptons[0].momentum()+cand_neutrinos[i].momentum();
	  numwm++;
	  break;
	}
	if   (  (cand_leptons[0].pid()+cand_neutrinos[i].pid())==1){//caso W+
	  wp=cand_leptons[0].momentum()+cand_neutrinos[i].momentum();
	  numwp++;
	  break;
	  //	}else{
	  //vetoEvent; 
	  }
      }
      // }//for j
  
      
      //per il leptone 2
      for (size_t  i=0; i< 2 ;++i){
	if   (  (cand_leptons[1].pid()+cand_neutrinos[i].pid())==-1){//caso W-
	  wm=cand_leptons[1].momentum()+cand_neutrinos[i].momentum();
	  numwm++;
	  break;
	}
	if   (  (cand_leptons[1].pid()+cand_neutrinos[i].pid())==1){//caso W+
	  wp=cand_leptons[1].momentum()+cand_neutrinos[i].momentum();
	  numwp++;
	  break;
	  //}else{
	  //vetoEvent; 
	  }
      }
        
      //Higgs
      double mh ;     
      double pth;  
      FourMomentum ll(0,0,0,0); 
      FourMomentum nn(0,0,0,0);
      ll= cand_leptons[0].momentum() + cand_leptons[1].momentum();
      nn = cand_neutrinos[0].momentum() + cand_neutrinos[1].momentum();          
      double mll  = ll.mass();         
      double ptll = ll.pT();  
      mh   = (wm+wp).mass();     
      pth=(wm+wp).pT();
      double etah=(wm+wp).eta();
      double mth  = sqrt( 2*ptll*eTmiss*( 1 - cos(fabs(deltaPhi(ll, pTmiss)))));    

      //W+W- mass
      double mwm =wm.mass();
      double mwp =wp.mass();

      //leptons
      double pt1 = cand_leptons[0].momentum().pT();//pt del leading lepton
      double pt2 = cand_leptons[1].momentum().pT();
      double eta_1l= cand_leptons[0].eta();
      double eta_2l= cand_leptons[1].eta();
      double pt_ll=ll.pT();//prima somma vettoriale, poi valuto il pt  
      double deta_l=eta_1l-eta_2l;
      double dphi_l=deltaPhi( eta_1l,  eta_2l);
      

      if(mwp<5. || mwm<5.){
	vetoEvent;        }
      if((nn.pt())<20.){
	vetoEvent;  }

      //per contare gli eventi selezionati (=> eventi in cui ho H->WW)
      numev++; 
      
      //Fill histrogram w/o request
      //Higgs
      _h_mh->fill(mh, weight); 
      _h_mth->fill(mth, weight);
      _h_pth->fill(pth, weight);
 
      if(fabs(etah) <5000000){
	_h_etah->fill(etah, weight); }

      //lepton
      _h_mll->fill(mll, weight);    
      _h_pt_1l->fill(pt1, weight);
      _h_pt_2l->fill(pt2, weight);
      _h_pt_ll->fill(pt_ll, weight);
      _h_eta_1l->fill(eta_1l, weight);
      _h_eta_2l->fill(eta_2l, weight);
      _h_deta_l->fill(deta_l, weight);
      _h_dphi_l->fill(dphi_l, weight);

      //W+-
      _h_mwm->fill(mwm, weight);
      _h_mwp->fill(mwp, weight);
     
      //Njets
      _h_njet->fill(njet, weight);
      if (mh < 130*GeV)
	_h_njet_onPeak->fill(njet, weight);
      else
	_h_njet_offPeak->fill(njet, weight);

      //Fill histrogram with Njet request
      if (njet == 0){
	_h_mll_0j->fill(mll, weight);
	_h_mth_0j->fill(mth, weight);
	_h_mh_0j->fill(mh, weight); 

	_h_mwm_0j->fill(mwm, weight);
	_h_mwp_0j->fill(mwp, weight);

      } else if (njet == 1){
   	_h_mll_1j->fill(mll, weight);
	_h_mth_1j->fill(mth, weight);
	_h_mh_1j->fill(mh, weight);
	double pt_1j= cand_jets[0].momentum().pT();
	double eta_1j= cand_jets[0].momentum().eta();
	_h_pt_1j->fill(pt_1j, weight);
	_h_eta_1j->fill(eta_1j, weight); 
	_h_mwm_1j->fill(mwm, weight);
	_h_mwp_1j->fill(mwp, weight);
      } else {
	_h_mll_2j->fill(mll, weight);
	_h_mth_2j->fill(mth, weight);
	_h_mh_2j->fill(mh, weight);	
	double pt_1j= cand_jets[0].momentum().pT();
	double pt_2j= cand_jets[1].momentum().pT();
	FourMomentum jj = cand_jets[0].momentum() + cand_jets[1].momentum();
	double mjj=jj.mass();
	double eta_1j= cand_jets[0].momentum().eta();
	double eta_2j= cand_jets[1].momentum().eta();
	double deta_j=eta_1j-eta_2j;
	_h_pt_1j->fill(pt_1j, weight);
	_h_pt_2j->fill(pt_2j, weight);
	_h_eta_1j->fill(eta_1j, weight); 
	_h_eta_2j->fill(eta_2j, weight); 
	_h_mjj->fill(mjj, weight);
	_h_deta_j->fill( deta_j,weight);   
    	_h_mwm_2j->fill(mwm, weight);
	_h_mwp_2j->fill(mwp, weight);


      }

      _h_weight->fill(weight);
      
    }//loop void analyze

    
    void finalize() {
      //  su sherpa final cross section = 6.246e-01 +- 2.214e-02 pb
      // double xsec=crossSection();
      // double xsec=sumOfWeights();
      double xsec=1.;
      std::cout << "xsec is: " << xsec << std::endl;
      
      scale(_h_njet,   xsec/sumOfWeights() );
      scale(_h_njet_onPeak, xsec/sumOfWeights());
      scale(_h_njet_offPeak, xsec/sumOfWeights());
      scale(_h_mll, xsec/sumOfWeights());
      scale(_h_mth, xsec/sumOfWeights());
      scale(_h_mh, xsec/sumOfWeights());

      scale(_h_mll_0j, xsec/sumOfWeights());
      scale(_h_mth_0j, xsec/sumOfWeights());
      scale(_h_mh_0j, xsec/sumOfWeights());

      scale(_h_mll_1j, xsec/sumOfWeights());
      scale(_h_mth_1j, xsec/sumOfWeights());
      scale(_h_mh_1j, xsec/sumOfWeights());

      scale(_h_mll_2j, xsec/sumOfWeights());
      scale(_h_mth_2j, xsec/sumOfWeights());
      scale(_h_mh_2j, xsec/sumOfWeights());

      scale(_h_mwm, xsec/sumOfWeights());
      scale(_h_mwp, xsec/sumOfWeights());
      
      scale(_h_mwm_0j, xsec/sumOfWeights());
      scale(_h_mwm_1j, xsec/sumOfWeights());
      scale(_h_mwm_2j, xsec/sumOfWeights());
      scale(_h_mwp_0j, xsec/sumOfWeights());
      scale(_h_mwp_1j, xsec/sumOfWeights());
      scale(_h_mwp_2j, xsec/sumOfWeights());
      

      //new
      scale(_h_pth, xsec/sumOfWeights());
      scale(_h_etah, xsec/sumOfWeights());   

      scale(_h_pt_1l, xsec/sumOfWeights());
      scale(_h_pt_2l, xsec/sumOfWeights());
      scale(_h_pt_ll, xsec/sumOfWeights());
      scale(_h_eta_1l, xsec/sumOfWeights());
      scale(_h_eta_2l, xsec/sumOfWeights());
      scale(_h_deta_l, xsec/sumOfWeights());
      scale(_h_dphi_l, xsec/sumOfWeights());

      scale(_h_mjj, xsec/sumOfWeights());
      scale(_h_pt_1j,xsec/sumOfWeights());
      scale(_h_pt_2j, xsec/sumOfWeights());  
      scale(_h_eta_1j, xsec/sumOfWeights());
      scale(_h_eta_2j, xsec/sumOfWeights());
      scale(_h_deta_j, xsec/sumOfWeights());
      scale(_h_weight, xsec/sumOfWeights());

            
      std::cout<< "num selected event with two W and so H= "<<numev<< std::endl;
      std::cout<< "num W-="<<numwm<<"  num W+= "<<numwp<< std::endl;
     
    }
   
  private:
    
    // Data members like post-cuts event weight counters go here
    
    
  private:
    
    //Higgs
    Histo1DPtr _h_mh;   
    Histo1DPtr _h_mh_0j;
    Histo1DPtr _h_mh_1j;
    Histo1DPtr _h_mh_2j;
    Histo1DPtr _h_mth;   
    Histo1DPtr _h_mth_0j;
    Histo1DPtr _h_mth_1j;
    Histo1DPtr _h_mth_2j;
    Histo1DPtr _h_pth;
    Histo1DPtr _h_etah;

    //leptons
    Histo1DPtr _h_mll;
    Histo1DPtr _h_mll_0j;
    Histo1DPtr _h_mll_1j;
    Histo1DPtr _h_mll_2j;
    Histo1DPtr _h_pt_1l;
    Histo1DPtr _h_pt_2l;
    Histo1DPtr _h_pt_ll;
    Histo1DPtr _h_eta_1l;
    Histo1DPtr _h_eta_2l;
    Histo1DPtr _h_deta_l;
    Histo1DPtr _h_dphi_l;   

    //W+- mass
    Histo1DPtr _h_mwm;
    Histo1DPtr _h_mwp;
    Histo1DPtr _h_mwm_0j;
    Histo1DPtr _h_mwm_1j;    
    Histo1DPtr _h_mwm_2j;
    Histo1DPtr _h_mwp_0j;
    Histo1DPtr _h_mwp_1j;
    Histo1DPtr _h_mwp_2j;

    //Jets
    Histo1DPtr _h_mjj;    
    Histo1DPtr _h_pt_1j;
    Histo1DPtr _h_pt_2j;
    Histo1DPtr _h_eta_1j;
    Histo1DPtr _h_eta_2j;
    Histo1DPtr _h_deta_j;
    Histo1DPtr _h_njet;
    Histo1DPtr _h_njet_onPeak;
    Histo1DPtr _h_njet_offPeak;

    Histo1DPtr _h_weight;

    ofstream myfile;
    double ce;
    double cmu;
    int numev;
    int numlep;
    int numneu;
    int numwp;
    int numwm;

  };
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(WWGen_cut_dl);
   
}


