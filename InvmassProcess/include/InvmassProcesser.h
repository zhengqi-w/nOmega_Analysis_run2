// varibles defines here
const TString TreeName_1 = "PiKDCascadeDecayTree";
TTree* fTreeChain;
const Double_t PionMass = 0.13957018;
const Double_t KaonMass = 0.493677;
const Double_t DeuteronMass = 1.875612928;
const Double_t ProtonMass = 0.9382720813;
const Int_t massbinnum = 30;
const Float_t massmin = 2.52;
const Float_t massmax = 2.62;
Int_t xpos=20, ypos=20;
Int_t cansizeX=1500, cansizeY=600;
//Branchs and BranchAddresses Setting
Int_t    T1DecayChannel;
Float_t  T1EventPercentile;
Float_t  T1TVDCA;
Float_t  T1TVCPA;
Float_t  T1TVtoPV;
Float_t  T1SVDCA;
Float_t  T1SVCPA;
Float_t  T1SVtoPV;
Float_t  T1SVtoTV;
Float_t  T1KaonPx;
Float_t  T1KaonPy;
Float_t  T1KaonPz;
Float_t  T1KaonDCAtoPVXY;
Float_t  T1KaonTPCNcls;
Float_t  T1KaonTPCnSigma;
Float_t  T1DeuteronPx;
Float_t  T1DeuteronPy;
Float_t  T1DeuteronPz;
Float_t  T1DeuteronDCAtoPVXY;
Float_t  T1DeuteronTPCNcls;
Float_t  T1DeuteronTPCnSigma;
Float_t  T1PionPx;
Float_t  T1PionPy;
Float_t  T1PionPz;
Float_t  T1PionDCAtoPVXY;
Float_t  T1PionTPCNcls;
Float_t  T1PionTPCnSigma;
//Branches
TBranch* b_T1DecayChannel;
TBranch* b_T1EventPercentile;
TBranch* b_T1TVDCA;
TBranch* b_T1TVCPA;
TBranch* b_T1TVtoPV;
TBranch* b_T1SVDCA;
TBranch* b_T1SVCPA;
TBranch* b_T1SVtoPV;
TBranch* b_T1SVtoTV;
TBranch* b_T1KaonPx;
TBranch* b_T1KaonPy;
TBranch* b_T1KaonPz;
TBranch* b_T1KaonDCAtoPVXY;
TBranch* b_T1KaonTPCNcls;
TBranch* b_T1KaonTPCnSigma;
TBranch* b_T1DeuteronPx;
TBranch* b_T1DeuteronPy;
TBranch* b_T1DeuteronPz;
TBranch* b_T1DeuteronDCAtoPVXY;
TBranch* b_T1DeuteronTPCNcls;
TBranch* b_T1DeuteronTPCnSigma;
TBranch* b_T1PionPx;
TBranch* b_T1PionPy;
TBranch* b_T1PionPz;
TBranch* b_T1PionDCAtoPVXY;
TBranch* b_T1PionTPCNcls;
TBranch* b_T1PionTPCnSigma;
//vectors for branchaddress
std::vector<Int_t*> PiKDCascadeDecayTree_BranchAddress_int = 
{
  &T1DecayChannel,
};

std::vector<Float_t*> PiKDCascadeDecayTree_BranchAddress_float = 
{
  &T1EventPercentile,
  &T1TVDCA,
  &T1TVCPA,
  &T1TVtoPV,
  &T1SVDCA,
  &T1SVCPA,
  &T1SVtoPV,
  &T1SVtoTV,
  &T1KaonPx,
  &T1KaonPy,
  &T1KaonPz,
  &T1KaonDCAtoPVXY,
  &T1KaonTPCNcls,
  &T1KaonTPCnSigma,
  &T1DeuteronPx,
  &T1DeuteronPy,
  &T1DeuteronPz,
  &T1DeuteronDCAtoPVXY,
  &T1DeuteronTPCNcls,
  &T1DeuteronTPCnSigma,
  &T1PionPx,
  &T1PionPy,
  &T1PionPz,
  &T1PionDCAtoPVXY,
  &T1PionTPCNcls,
  &T1PionTPCnSigma
};

std::vector<Double_t*> PiKDCascadeDecayTree_BranchAddress_double = 
{

};
//vector for branches
std::vector<TBranch*> PiKDCascadeDecayTree_Branches = 
{
  b_T1DecayChannel,
  b_T1EventPercentile,
  b_T1TVDCA,
  b_T1TVCPA,
  b_T1TVtoPV,
  b_T1SVDCA,
  b_T1SVCPA,
  b_T1SVtoPV,
  b_T1SVtoTV,
  b_T1KaonPx,
  b_T1KaonPy,
  b_T1KaonPz,
  b_T1KaonDCAtoPVXY,
  b_T1KaonTPCNcls,
  b_T1KaonTPCnSigma,
  b_T1DeuteronPx,
  b_T1DeuteronPy,
  b_T1DeuteronPz,
  b_T1DeuteronDCAtoPVXY,
  b_T1DeuteronTPCNcls,
  b_T1DeuteronTPCnSigma,
  b_T1PionPx,
  b_T1PionPy,
  b_T1PionPz,
  b_T1PionDCAtoPVXY,
  b_T1PionTPCNcls,
  b_T1PionTPCnSigma
};
//functions
void DrawSignalOnly(TH1* h_signal,TH1* h_signal_anti,TString dataset,TString target,TString beoraf);
void DrawSignalvsBKG(TH1* h_signal,TH1* h_signal_anti,TH1* h_bkg,TH1* h_bkg_anti,TString dataset,TString target, TString beofaf);
void Draw_Split_Invmass_Spectrum(std::vector<std::vector<TH1F*>> h_signal, std::vector<std::vector<TH1F*>> h_bkg, int xbins, int ybins, TString dataset, TString target, TString type);

