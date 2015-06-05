public:
double massJJ_min_WBF, deltayJJ_min_WBF;
double m_trans_llMET_min_WBF, m_ll_min_WBF;
double ptlep1_min_WBF,ptlep2_min_WBF,MET_min_WBF;

void initialiseHistos_WBF() {
  Histo1DPtr cuts_WBF;
  cuts_WBF     = bookHisto1D(5,-0.5,4.5);
  njets_before = bookHisto1D(10,-0.5,9.5);
  njets_after  = bookHisto1D(10,-0.5,9.5);
}

void analyse_WBF(const Event& event) {
  const double weight = event.weight();
  njets_before->fill(alljets.size(),weight);

  // cut on 2 jets in opposite hemispheres with minimal mass and rap distance
  if (alljets.size()<2) return;
  double y0(alljets[0]->momentum().eta());
  double y1(alljets[1]->momentum().eta());
  double massJJ((alljets[0]->momentum()+alljets[1]->momentum()).mass())
  if (y0*y1>0. || dabs(y0-y1)<deltaJJ_min_WBF || massJJ<massJJ_min_WBF) return;
  cuts_WBF->fill(2,weight);

  // cuts on 2 lepton MET system
  if (m_trans_llMET<m_trans_llMET_min_WBF || m_ll<m_ll_min_WBF) return;
  double ptm(lepton_m->momentum().pT()),ptp(lepton_p->momentum().pT());
  double ptlep1(ptm>ptp?ptm:ptp), ptlep2(ptm>ptp?ptp:ptm);
  if (ptlep1<ptlep1_min_WBF || ptlep2<ptlep2_min_WBF || MET<MET_min_WBF) return;
  cuts_WBF->fill(3,weight);
  
  njets_after->fill(alljets.size(),weight);

  // veto if tag jets are central (i.e. tagged) bjets
  if (bjets_central.find(alljets[0])!=bjet_central.end() ||
      bjets_central.find(alljets[1])!=bjet_central.end()) return
  cuts_WBF->fill(4,weight);
}
