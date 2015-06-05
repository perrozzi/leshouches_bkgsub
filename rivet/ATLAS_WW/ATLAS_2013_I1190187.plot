# BEGIN PLOT /MC_WHBB/d01-x01-y01
#Title=[Uncomment and insert title for histogram d01-x01-y01 here]
#XLabel=[Uncomment and insert x-axis label for histogram d01-x01-y01 here]
#YLabel=[Uncomment and insert y-axis label for histogram d01-x01-y01 here]
# + any additional plot settings you might like, see make-plots documentation
# END PLOT

# ... add more histograms as you need them ...
#
#
# BEGIN PLOT /MC_WHBB/*
XTwosidedTicks=1
YTwosidedTicks=1
ErrorBands=1
ErrorBars=1
LineStyle=solid
RatioPlotYLabel=ratio
#Rebin=4
# END PLOT

# BEGIN PLOT /MC_WHBB/Z_mass
Title=mass of Z-Bososn
XLabel=$M_Z$ (GeV)
YLabel=$\text{d}\sigma/\text{d}M_Z$ (pb/GeV)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/Z_eta
Title=psudorapidity of Z-Bososn
XLabel=$\eta(Z)$
YLabel=$\text{d}\sigma/\text{d}\eta(Z)$ (pb)
LogY=1
# END PLOT
#
# BEGIN PLOT /MC_WHBB/Z_pT
Title=transverse momentum of Z-Bososn
XLabel=$p_{T}(Z)$ (GeV)
YLabel=$\text{d}\sigma/\text{d}p_{T}(Z)$ (pb/GeV)
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_WHBB/HZ_pT
Title=transverse momentum of H-Z system
XLabel=$p_{T}(HZ)$ (GeV)
YLabel=$\text{d}\sigma/\text{d}p_{T}(HZ)$ (pb/GeV)
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_WHBB/HZ_eta
Title=psudorapidity of H-Z system
XLabel=$\eta(HZ)$
YLabel=$\text{d}\sigma/\text{d}\eta(HZ)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/dijet_mass
Title=mass of di-jet
XLabel=$M_{jj}$ (GeV)
YLabel=$\text{d}\sigma/\text{d}M_{jj}$ (pb/GeV)
LogY=0
# END PLOT

# BEGIN PLOT /MC_WHBB/dijet_eta
Title=psudorapidity of di-jet
XLabel=$\eta(jj)$
YLabel=$\text{d}\sigma/\text{d}\eta(jj)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/dijet_pT
Title=transverse momentum of di-jet
XLabel=$p_{T}(jj)$ (GeV)
YLabel=$\text{d}\sigma/\text{d}p_{T}(jj)$ (pb/GeV)
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_WHBB/h_dau_pT0
Title=transverse momentum of leading higgs jet
XLabel=$p_{T}(h_0)$ (GeV)
YLabel=$\text{d}\sigma/\text{d}p_{T}(h_0)$ (pb/GeV)
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_WHBB/h_dau_pT1
Title=transverse momentum of second higgs jet
XLabel=$p_{T}(h_1)$ (GeV)
YLabel=$\text{d}\sigma/\text{d}p_{T}(h_1)$ (pb/GeV)
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_WHBB/h_dau_eta0
Title=psudorapidity of leading higgs jet
XLabel=$\eta(h_0)$
YLabel=$\text{d}\sigma/\text{d}\eta(h_0)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/dR_jj
Title=distance between the two higgs jets
XLabel=$\Delta_{R}(h_0,h_1)$
YLabel=$\text{d}\sigma/\text{d}\Delta_{R}(h_0,h_1)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/deta_jj
Title=distance in $\eta$ between the two higgs jets
XLabel=$\Delta_{\eta}(h_0,h_1)$
YLabel=$\text{d}\sigma/\text{d}\Delta_{\eta}(h_0,h_1)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/dphi_jj
Title=distance in $\phi$ between the two higgs jets
XLabel=$\Delta_{\phi}(h_0,h_1)$
YLabel=$\text{d}\sigma/\text{d}\Delta_{\phi}(h_0,h_1)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/dR_HZ
Title=distance between di-jet and Z
XLabel=$\Delta_{R}(H,Z)$
YLabel=$\text{d}\sigma/\text{d}\Delta_{R}(H,Z)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/deta_HZ
Title=distance in $\eta$ between di-jet and Z
XLabel=$\Delta_{\eta}(H,Z)$
YLabel=$\text{d}\sigma/\text{d}\Delta_{\eta}(H,Z)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/dphi_HZ
Title=distance in $\phi$ between di-jet and Z
XLabel=$\Delta_{\phi}(h,Z)$
YLabel=$\text{d}\sigma/\text{d}\Delta_{\phi}(H,Z)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/h_dau_eta1
Title=psudorapidity of second higgs jet
XLabel=$\eta(h_1)$ (GeV)
YLabel=$\text{d}\sigma/\text{d}\eta(h_1)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/aJet_pT0
Title=transverse momentum of leading jet
XLabel=$p_{T}(jet_0)$ (GeV)
YLabel=$\text{d}\sigma/\text{d}p_{T}(jet0)$ (pb/GeV)
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_WHBB/aJet_pT1
Title=transverse momentum of second jet
XLabel=$p_{T}(jet_1)$ (GeV)
YLabel=$\text{d}\sigma/\text{d}p_{T}(jet1)$ (pb/GeV)
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_WHBB/aJet_pT2
Title=transverse momentum of third jet
XLabel=$p_{T}(jet_2)$ (GeV)
YLabel=$\text{d}\sigma/\text{d}p_{T}(jet2)$ (pb/GeV)
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_WHBB/aJet_pT3
Title=transverse momentum of fourth jet
XLabel=$p_{T}(jet_3)$ (GeV)
YLabel=$\text{d}\sigma/\text{d}p_{T}(jet3)$ (pb/GeV)
LogY=1
LogX=1
# END PLOT

# BEGIN PLOT /MC_WHBB/aJet_eta0
Title=psudorapidity of leading jet
XLabel=$\eta(jet_0)$
YLabel=$\text{d}\sigma/\text{d}\eta(jet_0)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/aJet_eta1
Title=psudorapidity of second jet
XLabel=$\eta(jet_1)$
YLabel=$\text{d}\sigma/\text{d}\eta(jet_1)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/aJet_eta2
Title=psudorapidity of third jet
XLabel=$\eta(jet_2)$
YLabel=$\text{d}\sigma/\text{d}\eta(jet_2)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/aJet_eta3
Title=psudorapidity of fourth jet
XLabel=$\eta(jet_3)$
YLabel=$\text{d}\sigma/\text{d}\eta(jet_3)$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/najets
Title=number of additional jets
XLabel=$N_{aJets}$
YLabel=$\text{d}\sigma/\text{d}N_{aJets}$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/nbjets
Title=number of b-jets
XLabel=$N_{bJets}$
YLabel=$\text{d}\sigma/\text{d}N_{bJets}$ (pb)
LogY=1
# END PLOT

# BEGIN PLOT /MC_WHBB/log10_d_01
Title=$\log_{10}$($k_\perp$ jet resolution $0 \to 1$ (GeV))
XLabel=$\log_{10}(d_{01}/\text{GeV})$
YLabel=$\text{d}\sigma/\text{d}\log_{10}(d_{01}/\text{GeV})$ (pb)
# END PLOT

# BEGIN PLOT /MC_WHBB/log10_d_12
Title=$\log_{10}$($k_\perp$ jet resolution $1 \to 2$ (GeV))
XLabel=$\log_{10}(d_{12}/\text{GeV})$
YLabel=$\text{d}\sigma/\text{d}\log_{10}(d_{12}/\text{GeV})$ (pb)
# END PLOT

# BEGIN PLOT /MC_WHBB/log10_d_23
Title=$\log_{10}$($k_\perp$ jet resolution $2 \to 3$ (GeV))
XLabel=$\log_{10}(d_{23}/\text{GeV})$
YLabel=$\text{d}\sigma/\text{d}\log_{10}(d_{23}/\text{GeV})$ (pb)
# END PLOT

# BEGIN PLOT /MC_WHBB/log10_d_34
Title=$\log_{10}$($k_\perp$ jet resolution $3 \to 4$ (GeV))
XLabel=$\log_{10}(d_{34}/\text{GeV})$
YLabel=$\text{d}\sigma/\text{d}\log_{10}(d_{34}/\text{GeV})$ (pb)
# END PLOT
