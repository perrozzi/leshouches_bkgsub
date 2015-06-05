public:
double mass_ll_min_hh, mBB_min_hh, mBB_max_hh;
double METrel_min_hh;
double bjets_central_min_hh;
double m_trans_llMET_min_hh, m_trans_llMET_max_hh;

void initialiseHistos_hh() {
  Histo1DPtr cuts_hh; 
  Histo2DPtr njets_nbjets_before_hh, njets_nbjets_after_hh;
  cuts_hh     = bookHisto1D("cuts_hh", 8,-0.5,7.5);
  njets_nbjets_before_hh     = bookHisto2D("njets_nbjets_before_hh",10,-0.5,9.5, 10,-0.5,9.5);
  njets_nbjets_after_hh     = bookHisto2D("njets_nbjets_after_hh",10,-0.5,9.5, 10,-0.5,9.5);      
}

void analyse_hh(const Event& event) {
  const double weight = event.weight();

  // Objects used: MET, lepton_m, lepton_p, alljets, bjets_central, mass_ll, m_trans_llMET

  // Missing ETrel
  //            *"mismom"   for "delta_phi" >= (0.5*pi)
  //            *"mismom.pT()*sin(delta_phi)"   for "delta_phi" < (0.5*pi)
  FourMomentum mismom;
  double METrel = 0, delta_phi = 0;
  vector<double> vL_MET_angle, vJet_MET_angle;
  mismom = -MET.visibleMomentum();
  vL_MET_angle.push_back(fabs(deltaPhi(lepton_m.momentum(), mismom)));
  vL_MET_angle.push_back(fabs(deltaPhi(lepton_p.momentum(), mismom)));
  foreach (double& lM, vL_MET_angle) if (lM > M_PI) lM = 2*M_PI - lM;
  std::sort(vL_MET_angle.begin(), vL_MET_angle.end());
  if (alljets.size() == 0) delta_phi = vL_MET_angle[0];
  if (alljets.size() > 0) {
    foreach (Jet& vj, alljets) {
      double jet_MET_angle = fabs(deltaPhi(vj.momentum(), mismom));
      if (jet_MET_angle > M_PI) jet_MET_angle = 2*M_PI - jet_MET_angle;
      vJet_MET_angle.push_back(jet_MET_angle);
    }
    std::sort(vJet_MET_angle.begin(), vJet_MET_angle.end());
    if (vL_MET_angle[0] <= vJet_MET_angle[0]) delta_phi = vL_MET_angle[0];
    if (vL_MET_angle[0] > vJet_MET_angle[0]) delta_phi = vJet_MET_angle[0];
  }  
  if (delta_phi >= (0.5*M_PI)) delta_phi = 0.5*M_PI;
  METrel = mismom.pT()*sin(delta_phi);  

  // Selection
  if(METrel <= METrel_min_hh) return;
  njets_nbjets_before_hh->fill(alljets.size(), bjets_central.size());
  cuts_hh->fill(2,weight);

  if(mass_ll <= mass_ll_min_hh) return;
  cuts_hh->fill(3,weight);

  if(bjets_central.size() < bjets_central_min_hh) return;
  cuts_hh->fill(4,weight);

  // 2 leading b-jets invariant mass
  double mBB = (bjets_central[0].momentum() + bjets_central[1].momentum()).mass();

  if(m_trans_llMET < m_trans_llMET_min_hh || m_trans_llMET > m_trans_llMET_max_hh) return;
  cuts_hh->fill(5,weight);

  if(mBB < mBB_min_hh || mBB > mBB_max_hh) return;
  njets_nbjets_after_hh->fill(alljets.size(), bjets_central.size());
  cuts_hh->fill(6,weight);

}
