#include "AliAnalysisTaskHyperFinder3Body.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TMath.h"
#include <TString.h>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliESDcascade.h"
#include "AliTrackerBase.h"
#include <vector>
#include <unordered_set>
#include <TLorentzVector.h>
#include "AliVVertex.h"

class AliAnalysisTaskHyperFinder3Body;

using namespace std;

ClassImp(AliAnalysisTaskHyperFinder3Body);

AliAnalysisTaskHyperFinder3Body::AliAnalysisTaskHyperFinder3Body() : AliAnalysisTaskSE(),
fInputHandler(),
fESDevent(),
fOutputList(),
fESDtrackCuts(),
feventCut(),
fPIDResponse(),
nRotation(2),
TrackDCAtoPVcut(0.05),
TrackTPCRefit(kTRUE),
TrackAcceptKinkDaughters(kFALSE),
kOpenTrackImpacCutV0(kFALSE),
kRotatePion(kFALSE),
kRotateKaon(kFALSE),
kRotateDeuteron(kTRUE),
TrackTPCNClscut(50),
TrackLengthcut(80),
TrackNCrossedRowscut(60),
TrackEtacut(0.8),
TrackChi2PerClusterTPCcut(7),
V0Radiusmin(0.9),
V0Radiusmax(200),
V0DCAdaughtercut(1.2),//1.5
V0CPAcut(0.995),//0.98
V0DCAtoPVcut(0.1),
CascadeDCAcut(1.2),//1.5
CascadeV0CPAcut(0.995),//0.98
CascadeCPAcut(0.995),//0.98
CascadeRadiusmin(0.9),
CascadeRadiusmax(100),
Cascadeptmin(0.3),
Cascadeptmax(100),
CascadeV0Invmassmax(2.0553),
CascadeInvmassmin(2.52),
CascadeInvmassmax(2.62),
CascadeInvmassmin2(3.13),
CascadeInvmassmax2(3.23),
V0Invmassmin(2.02),
V0Invmassmax(2.14),
V0Invmassmin2(2.07),
V0Invmassmax2(2.13),
fhEventCount(),
fhcentralitybeforecut(),
fhcentralityaftercut(),
fhCausalitycheck(),
fhTrackDCAtoPVxycheck(),
fhPiKDChargeCheckC1(),
fhPiKDChargeCheckC2(),
fhPiKDChargeCheckC3(),
fhPiKDChargeCheckC4(),
fhPiKDChargeCheckC5(),
fhPiKDChargeCheckC6(),
fhPiKDChargeCheckC7(),
fhPiKDChargeCheckC8(),
fhPiDChargeCheckC1(),
fhPiDChargeCheckC2(),
fhPiDChargeCheckC3(),
fhPiDChargeCheckC4(),
fhPropagationcheck(),
fhCascadePropagationPcheck(),
fhDeuteronEventCheck(),
fhTrigger(),
fhBBPion(),
fhBBDeuteron(),
fhBBKaon(),
PiKDCascadeDecayTree(),
T1DecayChannel(0),
T1EventPercentile(0),
T1TVDCA(0),
T1TVCPA(0),
T1TVtoPV(0),
T1SVDCA(0),
T1SVCPA(0),
T1SVtoPV(0),
T1SVtoTV(0),
T1KaonPx(0),
T1KaonPy(0),
T1KaonPz(0),
T1KaonDCAtoPVXY(0),
T1KaonTPCNcls(0),
T1KaonTPCnSigma(0),
T1DeuteronPx(0),
T1DeuteronPy(0),
T1DeuteronPz(0),
T1DeuteronDCAtoPVXY(0),
T1DeuteronTPCNcls(0),
T1DeuteronTPCnSigma(0),
T1PionPx(0),
T1PionPy(0),
T1PionPz(0),
T1PionDCAtoPVXY(0),
T1PionTPCNcls(0),
T1PionTPCnSigma(0),
PiDDecayTree(),
T2DecayChannel(0),
T2EventPercentile(0),
T2SVDCA(0),
T2SVCPA(0),
T2SVtoPV(0),
T2DeuteronPx(0),
T2DeuteronPy(0),
T2DeuteronPz(0),
T2DeuteronDCAtoPVXY(0),
T2DeuteronDCAtoPVZ(0),
T2DeuteronDCAtoSVXY(0),
T2DeuteronTPCNcls(0),
T2DeuteronTPCnSigma(0),
T2PionPx(0),
T2PionPy(0),
T2PionPz(0),
T2PionDCAtoPVXY(0),
T2PionDCAtoPVZ(0),
T2PionDCAtoSVXY(0),
T2PionTPCNcls(0),
T2PionTPCnSigma(0)
{

} 
AliAnalysisTaskHyperFinder3Body::AliAnalysisTaskHyperFinder3Body(const char* name) : AliAnalysisTaskSE(name),
fInputHandler(),
fESDevent(),
fOutputList(),
fESDtrackCuts(),
feventCut(),
fPIDResponse(),
nRotation(2),
TrackDCAtoPVcut(0.05),
TrackTPCRefit(kTRUE),
TrackAcceptKinkDaughters(kFALSE),
kOpenTrackImpacCutV0(kFALSE),
kRotatePion(kFALSE),
kRotateKaon(kFALSE),
kRotateDeuteron(kTRUE),
TrackTPCNClscut(50),
TrackLengthcut(80),
TrackNCrossedRowscut(60),
TrackEtacut(0.8),
TrackChi2PerClusterTPCcut(7),
V0Radiusmin(0.9),
V0Radiusmax(200),
V0DCAdaughtercut(1.2),//1.5
V0CPAcut(0.995),//0.98
V0DCAtoPVcut(0.1),
CascadeDCAcut(1.2),//1.5
CascadeV0CPAcut(0.995),//0.98
CascadeCPAcut(0.995),//0.98
CascadeRadiusmin(0.9),
CascadeRadiusmax(100),
Cascadeptmin(0.3),
Cascadeptmax(100),
CascadeV0Invmassmax(2.0553),
CascadeInvmassmin(2.52),
CascadeInvmassmax(2.62),
CascadeInvmassmin2(3.13),
CascadeInvmassmax2(3.23),
V0Invmassmin(2.02),
V0Invmassmax(2.14),
V0Invmassmin2(2.07),
V0Invmassmax2(2.13),
fhEventCount(),
fhcentralitybeforecut(),
fhcentralityaftercut(),
fhCausalitycheck(),
fhTrackDCAtoPVxycheck(),
fhPiKDChargeCheckC1(),
fhPiKDChargeCheckC2(),
fhPiKDChargeCheckC3(),
fhPiKDChargeCheckC4(),
fhPiKDChargeCheckC5(),
fhPiKDChargeCheckC6(),
fhPiKDChargeCheckC7(),
fhPiKDChargeCheckC8(),
fhPiDChargeCheckC1(),
fhPiDChargeCheckC2(),
fhPiDChargeCheckC3(),
fhPiDChargeCheckC4(),
fhPropagationcheck(),
fhCascadePropagationPcheck(),
fhDeuteronEventCheck(),
fhTrigger(),
fhBBPion(),
fhBBDeuteron(),
fhBBKaon(),
PiKDCascadeDecayTree(),
T1DecayChannel(0),
T1EventPercentile(0),
T1TVDCA(0),
T1TVCPA(0),
T1TVtoPV(0),
T1SVDCA(0),
T1SVCPA(0),
T1SVtoPV(0),
T1SVtoTV(0),
T1KaonPx(0),
T1KaonPy(0),
T1KaonPz(0),
T1KaonDCAtoPVXY(0),
T1KaonTPCNcls(0),
T1KaonTPCnSigma(0),
T1DeuteronPx(0),
T1DeuteronPy(0),
T1DeuteronPz(0),
T1DeuteronDCAtoPVXY(0),
T1DeuteronTPCNcls(0),
T1DeuteronTPCnSigma(0),
T1PionPx(0),
T1PionPy(0),
T1PionPz(0),
T1PionDCAtoPVXY(0),
T1PionTPCNcls(0),
T1PionTPCnSigma(0),
PiDDecayTree(),
T2DecayChannel(0),
T2EventPercentile(0),
T2SVDCA(0),
T2SVCPA(0),
T2SVtoPV(0),
T2DeuteronPx(0),
T2DeuteronPy(0),
T2DeuteronPz(0),
T2DeuteronDCAtoPVXY(0),
T2DeuteronDCAtoPVZ(0),
T2DeuteronDCAtoSVXY(0),
T2DeuteronTPCNcls(0),
T2DeuteronTPCnSigma(0),
T2PionPx(0),
T2PionPy(0),
T2PionPz(0),
T2PionDCAtoPVXY(0),
T2PionDCAtoPVZ(0),
T2PionDCAtoSVXY(0),
T2PionTPCNcls(0),
T2PionTPCnSigma(0)
{
  fESDtrackCuts = new AliESDtrackCuts("fESDtrackCuts");
  fESDtrackCuts->SetRequireTPCRefit(TrackTPCRefit);
  fESDtrackCuts->SetAcceptKinkDaughters(TrackAcceptKinkDaughters);
  fESDtrackCuts->SetMinNClustersTPC(TrackTPCNClscut);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(TrackChi2PerClusterTPCcut);
  fESDtrackCuts->SetEtaRange(-1.*TrackEtacut, TrackEtacut);
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}
AliAnalysisTaskHyperFinder3Body::~AliAnalysisTaskHyperFinder3Body()
{
  if(fInputHandler)delete fInputHandler;
  if(fESDevent)delete fESDevent;
  if(fOutputList)delete fOutputList;
  if(fESDtrackCuts)delete fESDtrackCuts;
  if(fPIDResponse)delete fPIDResponse;
  if(fhEventCount)delete fhEventCount;
  if(fhcentralitybeforecut)delete fhcentralitybeforecut;
  if(fhcentralityaftercut)delete fhcentralityaftercut;
  if(fhCausalitycheck)delete fhCausalitycheck;
  if(fhTrackDCAtoPVxycheck)delete fhTrackDCAtoPVxycheck;
  if(fhPiKDChargeCheckC1)delete fhPiKDChargeCheckC1;
  if(fhPiKDChargeCheckC2)delete fhPiKDChargeCheckC2;
  if(fhPiKDChargeCheckC3)delete fhPiKDChargeCheckC3;
  if(fhPiKDChargeCheckC4)delete fhPiKDChargeCheckC4;
  if(fhPiKDChargeCheckC5)delete fhPiKDChargeCheckC5;
  if(fhPiKDChargeCheckC6)delete fhPiKDChargeCheckC6;
  if(fhPiKDChargeCheckC7)delete fhPiKDChargeCheckC7;
  if(fhPiKDChargeCheckC8)delete fhPiKDChargeCheckC8;
  if(fhPiDChargeCheckC1)delete fhPiDChargeCheckC1;
  if(fhPiDChargeCheckC2)delete fhPiDChargeCheckC2;
  if(fhPiDChargeCheckC3)delete fhPiDChargeCheckC3;
  if(fhPiDChargeCheckC4)delete fhPiDChargeCheckC4;
  if(fhPropagationcheck)delete fhPropagationcheck;
  if(fhCascadePropagationPcheck)delete fhCascadePropagationPcheck;
  if(fhDeuteronEventCheck)delete fhDeuteronEventCheck;
  if(fhTrigger)delete fhTrigger;
  if(fhBBPion)delete fhBBPion;
  if(fhBBDeuteron)delete fhBBDeuteron;
  if(fhBBKaon)delete fhBBKaon;
  if(PiKDCascadeDecayTree)delete PiKDCascadeDecayTree;
  if(PiDDecayTree)delete PiDDecayTree;
}
void AliAnalysisTaskHyperFinder3Body::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  if(!fhEventCount)
  {
    fhEventCount = new TH1I("fhEventCount","event counter",3,0,3);
    fhEventCount->GetXaxis()->SetBinLabel(1,"beginning");
    fhEventCount->GetXaxis()->SetBinLabel(2,"read in");
    fhEventCount->GetXaxis()->SetBinLabel(3,"after cuts");
    fOutputList->Add(fhEventCount);
  }
  if(!fhcentralityaftercut)
  {
    fhcentralitybeforecut = new TH1F("fhcentralitybeforecut","centrality distribution before event cuts",100,0,100);
    fOutputList->Add(fhcentralitybeforecut);
  }
  if(!fhcentralityaftercut)
  {
    fhcentralityaftercut = new TH1F("fhcentralityaftercut","centrality distribution after even cuts",100,0,100);
    fOutputList->Add(fhcentralityaftercut);
  }
  if(!fhCausalitycheck)
  {
    fhCausalitycheck = new TH1I("fhCausalitycheck"," before and after casality cut",2,0,2);
    fOutputList->Add(fhCausalitycheck);
  }
  if(!fhTrackDCAtoPVxycheck)
  {
    fhTrackDCAtoPVxycheck = new TH1D("fhTrackDCAtoPVxycheck","checking for trackdcaxy (impact - d)",100,-2,2);
    fOutputList->Add(fhTrackDCAtoPVxycheck);
  }
  if(!fhPiKDChargeCheckC1)
  {
    fhPiKDChargeCheckC1 = new TH2I("fhPiKDChargeCheckC1","PiKDChargeCheckC1",3,0,3,4,-2,2);
    fhPiKDChargeCheckC1->GetXaxis()->SetTitle("particle type");
    fhPiKDChargeCheckC1->GetXaxis()->SetBinLabel(1,"pion");
    fhPiKDChargeCheckC1->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiKDChargeCheckC1->GetXaxis()->SetBinLabel(3,"kaon");
    fhPiKDChargeCheckC1->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiKDChargeCheckC1);
  }
  if(!fhPiKDChargeCheckC2)
  {
    fhPiKDChargeCheckC2 = new TH2I("fhPiKDChargeCheckC2","PiKDChargeCheckC2",3,0,3,4,-2,2);
    fhPiKDChargeCheckC2->GetXaxis()->SetTitle("particle type");
    fhPiKDChargeCheckC2->GetXaxis()->SetBinLabel(1,"pion");
    fhPiKDChargeCheckC2->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiKDChargeCheckC2->GetXaxis()->SetBinLabel(3,"kaon");
    fhPiKDChargeCheckC2->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiKDChargeCheckC2);
  }
  if(!fhPiKDChargeCheckC3)
  {
    fhPiKDChargeCheckC3 = new TH2I("fhPiKDChargeCheckC3","PiKDChargeCheckC3",3,0,3,4,-2,2);
    fhPiKDChargeCheckC3->GetXaxis()->SetTitle("particle type");
    fhPiKDChargeCheckC3->GetXaxis()->SetBinLabel(1,"pion");
    fhPiKDChargeCheckC3->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiKDChargeCheckC3->GetXaxis()->SetBinLabel(3,"kaon");
    fhPiKDChargeCheckC3->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiKDChargeCheckC3);
  }
  if(!fhPiKDChargeCheckC4)
  {
    fhPiKDChargeCheckC4 = new TH2I("fhPiKDChargeCheckC4","PiKDChargeCheckC4",3,0,3,4,-2,2);
    fhPiKDChargeCheckC4->GetXaxis()->SetTitle("particle type");
    fhPiKDChargeCheckC4->GetXaxis()->SetBinLabel(1,"pion");
    fhPiKDChargeCheckC4->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiKDChargeCheckC4->GetXaxis()->SetBinLabel(3,"kaon");
    fhPiKDChargeCheckC4->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiKDChargeCheckC4);
  }
  if(!fhPiKDChargeCheckC5)
  {
    fhPiKDChargeCheckC5 = new TH2I("fhPiKDChargeCheckC5","PiKDChargeCheckC5",3,0,3,4,-2,2);
    fhPiKDChargeCheckC5->GetXaxis()->SetTitle("particle type");
    fhPiKDChargeCheckC5->GetXaxis()->SetBinLabel(1,"pion");
    fhPiKDChargeCheckC5->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiKDChargeCheckC5->GetXaxis()->SetBinLabel(3,"kaon");
    fhPiKDChargeCheckC5->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiKDChargeCheckC5);
  }
  if(!fhPiKDChargeCheckC6)
  {
    fhPiKDChargeCheckC6 = new TH2I("fhPiKDChargeCheckC6","PiKDChargeCheckC6",3,0,3,4,-2,2);
    fhPiKDChargeCheckC6->GetXaxis()->SetTitle("particle type");
    fhPiKDChargeCheckC6->GetXaxis()->SetBinLabel(1,"pion");
    fhPiKDChargeCheckC6->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiKDChargeCheckC6->GetXaxis()->SetBinLabel(3,"kaon");
    fhPiKDChargeCheckC6->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiKDChargeCheckC6);
  }
  if(!fhPiKDChargeCheckC7)
  {
    fhPiKDChargeCheckC7 = new TH2I("fhPiKDChargeCheckC7","PiKDChargeCheckC7",3,0,3,4,-2,2);
    fhPiKDChargeCheckC7->GetXaxis()->SetTitle("particle type");
    fhPiKDChargeCheckC7->GetXaxis()->SetBinLabel(1,"pion");
    fhPiKDChargeCheckC7->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiKDChargeCheckC7->GetXaxis()->SetBinLabel(3,"kaon");
    fhPiKDChargeCheckC7->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiKDChargeCheckC7);
  }
  if(!fhPiKDChargeCheckC8)
  {
    fhPiKDChargeCheckC8 = new TH2I("fhPiKDChargeCheckC8","PiKDChargeCheckC8",3,0,3,4,-2,2);
    fhPiKDChargeCheckC8->GetXaxis()->SetTitle("particle type");
    fhPiKDChargeCheckC8->GetXaxis()->SetBinLabel(1,"pion");
    fhPiKDChargeCheckC8->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiKDChargeCheckC8->GetXaxis()->SetBinLabel(3,"kaon");
    fhPiKDChargeCheckC8->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiKDChargeCheckC8);
  }
  if(!fhPiDChargeCheckC1)
  {
    fhPiDChargeCheckC1 = new TH2I("fhPiDChargeCheckC1","PiDChargeCheckC1",2,0,3,4,-2,2);
    fhPiDChargeCheckC1->GetXaxis()->SetTitle("particle type");
    fhPiDChargeCheckC1->GetXaxis()->SetBinLabel(1,"pion");
    fhPiDChargeCheckC1->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiDChargeCheckC1->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiDChargeCheckC1);
  }
  if(!fhPiDChargeCheckC2)
  {
    fhPiDChargeCheckC2 = new TH2I("fhPiDChargeCheckC2","PiDChargeCheckC2",2,0,3,4,-2,2);
    fhPiDChargeCheckC2->GetXaxis()->SetTitle("particle type");
    fhPiDChargeCheckC2->GetXaxis()->SetBinLabel(1,"pion");
    fhPiDChargeCheckC2->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiDChargeCheckC2->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiDChargeCheckC2);
  }
  if(!fhPiDChargeCheckC3)
  {
    fhPiDChargeCheckC3 = new TH2I("fhPiDChargeCheckC3","PiDChargeCheckC3",2,0,3,4,-2,2);
    fhPiDChargeCheckC3->GetXaxis()->SetTitle("particle type");
    fhPiDChargeCheckC3->GetXaxis()->SetBinLabel(1,"pion");
    fhPiDChargeCheckC3->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiDChargeCheckC3->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiDChargeCheckC3);
  }
  if(!fhPiDChargeCheckC4)
  {
    fhPiDChargeCheckC4 = new TH2I("fhPiDChargeCheckC4","PiDChargeCheckC4",2,0,3,4,-2,2);
    fhPiDChargeCheckC4->GetXaxis()->SetTitle("particle type");
    fhPiDChargeCheckC4->GetXaxis()->SetBinLabel(1,"pion");
    fhPiDChargeCheckC4->GetXaxis()->SetBinLabel(2,"deuteron");
    fhPiDChargeCheckC4->GetYaxis()->SetTitle("Charge");
    fOutputList->Add(fhPiDChargeCheckC4);
  }
  
  if(!fhPropagationcheck)
  {
    fhPropagationcheck = new TH1I("fhPropagationcheck","fhPropagationcheck",2,0,2);
    fhPropagationcheck->GetXaxis()->SetBinLabel(1, "failure");
    fhPropagationcheck->GetXaxis()->SetBinLabel(2, "success");
    fOutputList->Add(fhPropagationcheck);
  }
  if(!fhCascadePropagationPcheck)
  {
    fhCascadePropagationPcheck = new TH2D("fhCascadePropagationPcheck","fhCascadePropagationPcheck",3,0,3,1000,-0.1,0.1);
    fhCascadePropagationPcheck->GetXaxis()->SetTitle("p type");
    fhCascadePropagationPcheck->GetXaxis()->SetBinLabel(1,"px");
    fhCascadePropagationPcheck->GetXaxis()->SetBinLabel(2,"py");
    fhCascadePropagationPcheck->GetXaxis()->SetBinLabel(3,"pz");
    fhCascadePropagationPcheck->GetYaxis()->SetTitle("deviation between track->GetPAt and trackin->P()");
    fOutputList->Add(fhCascadePropagationPcheck);
  }
  if(!fhDeuteronEventCheck)
  {
    fhDeuteronEventCheck = new TH1I("fhDeuteronEventCheck","Events has vs has no deuteron",2,0,2);
    fhDeuteronEventCheck->GetXaxis()->SetBinLabel(1,"with D");
    fhDeuteronEventCheck->GetXaxis()->SetBinLabel(2,"without D");
    fOutputList->Add(fhDeuteronEventCheck);
  }
  if(!fhTrigger)
  {
    fhTrigger = new TH1F("fhTrigger", "Trigger", 8, 0, 8);
    fhTrigger->GetXaxis()->SetBinLabel(1, "other");
    fhTrigger->GetXaxis()->SetBinLabel(2, "kINT7");
    fhTrigger->GetXaxis()->SetBinLabel(3, "kHighMultV0");
    fhTrigger->GetXaxis()->SetBinLabel(4, "kHighMultSPD");
    fhTrigger->GetXaxis()->SetBinLabel(5, "HNU");
    fhTrigger->GetXaxis()->SetBinLabel(6, "HQU");
    fhTrigger->GetXaxis()->SetBinLabel(7, "kCentral");
    fhTrigger->GetXaxis()->SetBinLabel(8, "kSemiCentral");
    fOutputList->Add(fhTrigger);
  }
  if(!fhBBPion)
  {
    fhBBPion = new TH2D("fhBBPion","Pion dE/dX",240,-6,6,300,0,1500);
    fhBBPion->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhBBPion->GetYaxis()->SetTitle("TPC Signal");
    fOutputList->Add(fhBBPion);
  }
  if(!fhBBDeuteron)
  {
    fhBBDeuteron = new TH2D("fhBBDeuteron","Deuteron dE/dx",240,-6,6,300,0,1500);
    fhBBDeuteron->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhBBDeuteron->GetYaxis()->SetTitle("TPC Signal");
    fOutputList->Add(fhBBDeuteron);
  }
  if(!fhBBKaon)
  {
    fhBBKaon = new TH2D("fhBBKaon","Kaon dE/dX",240,-6,6,300,0,1500);
    fhBBKaon->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhBBKaon->GetYaxis()->SetTitle("TPC Signal");
    fOutputList->Add(fhBBKaon);
  }
  feventCut.AddQAplotsToList(fOutputList);
  if(!PiKDCascadeDecayTree)
  {
    PiKDCascadeDecayTree = new TTree("PiKDCascadeDecayTree","PiKDCascadeDecayTree");
    PiKDCascadeDecayTree->Branch("T1DecayChannel",&T1DecayChannel,"T1DecayChannel/I");
    PiKDCascadeDecayTree->Branch("T1EventPercentile",&T1EventPercentile,"T1EventPercentile/F");
    PiKDCascadeDecayTree->Branch("T1TVDCA",&T1TVDCA,"T1TVDCA/F");
    PiKDCascadeDecayTree->Branch("T1TVCPA",&T1TVCPA,"T1TVCPA/F");
    PiKDCascadeDecayTree->Branch("T1TVtoPV",&T1TVtoPV,"T1TVtoPV/F");
    PiKDCascadeDecayTree->Branch("T1SVDCA",&T1SVDCA,"T1SVDCA/F");
    PiKDCascadeDecayTree->Branch("T1SVCPA",&T1SVCPA,"T1SVCPA/F");
    PiKDCascadeDecayTree->Branch("T1SVtoPV",&T1SVtoPV,"T1SVtoPV/F");
    PiKDCascadeDecayTree->Branch("T1SVtoTV",&T1SVtoTV,"T1SVtoTV/F");
    PiKDCascadeDecayTree->Branch("T1KaonPx",&T1KaonPx,"T1KaonPx/F");
    PiKDCascadeDecayTree->Branch("T1KaonPy",&T1KaonPy,"T1KaonPy/F");
    PiKDCascadeDecayTree->Branch("T1KaonPz",&T1KaonPz,"T1KaonPz/F");
    PiKDCascadeDecayTree->Branch("T1KaonDCAtoPVXY",&T1KaonDCAtoPVXY,"T1KaonDCAtoPVXY/F");
    PiKDCascadeDecayTree->Branch("T1KaonTPCNcls",&T1KaonTPCNcls,"T1KaonTPCNcls/F");
    PiKDCascadeDecayTree->Branch("T1KaonTPCnSigma",&T1KaonTPCnSigma,"T1KaonTPCnSigma/F");
    PiKDCascadeDecayTree->Branch("T1DeuteronPx",&T1DeuteronPx,"T1DeuteronPx/F");
    PiKDCascadeDecayTree->Branch("T1DeuteronPy",&T1DeuteronPy,"T1DeuteronPy/F");
    PiKDCascadeDecayTree->Branch("T1DeuteronPz",&T1DeuteronPz,"T1DeuteronPz/F");
    PiKDCascadeDecayTree->Branch("T1DeuteronDCAtoPVXY",&T1DeuteronDCAtoPVXY,"T1DeuteronDCAtoPVXY/F");
    PiKDCascadeDecayTree->Branch("T1DeuteronTPCNcls",&T1DeuteronTPCNcls,"T1DeuteronTPCNcls/F");
    PiKDCascadeDecayTree->Branch("T1DeuteronTPCnSigma",&T1DeuteronTPCnSigma,"T1DeuteronTPCnSigma/F");
    PiKDCascadeDecayTree->Branch("T1PionPx",&T1PionPx,"T1PionPx/F");
    PiKDCascadeDecayTree->Branch("T1PionPy",&T1PionPy,"T1PionPy/F");
    PiKDCascadeDecayTree->Branch("T1PionPz",&T1PionPz,"T1PionPz/F");
    PiKDCascadeDecayTree->Branch("T1PionDCAtoPVXY",&T1PionDCAtoPVXY,"T1PionDCAtoPVXY/F");
    PiKDCascadeDecayTree->Branch("T1PionTPCNcls",&T1PionTPCNcls,"T1PionTPCNcls/F");
    PiKDCascadeDecayTree->Branch("T1PionTPCnSigma",&T1PionTPCnSigma,"T1PionTPCnSigma/F");
  }
  if(!PiDDecayTree)
  {
    PiDDecayTree = new TTree("PiDDecayTree","PiDDecayTree");
    PiDDecayTree->Branch("T2DecayChannel",&T2DecayChannel,"T2DecayChannel/I");
    PiDDecayTree->Branch("T2EventPercentile",&T2EventPercentile,"T2EventPercentile/F");
    PiDDecayTree->Branch("T2SVDCA",&T2SVDCA,"T2SVDCA/F");
    PiDDecayTree->Branch("T2SVCPA",&T2SVCPA,"T2SVCPA/F");
    PiDDecayTree->Branch("T2SVtoPV",&T2SVtoPV,"T2SVtoPV/F");
    PiDDecayTree->Branch("T2DeuteronPx",&T2DeuteronPx,"T2DeuteronPx/F");
    PiDDecayTree->Branch("T2DeuteronPy",&T2DeuteronPy,"T2DeuteronPy/F");
    PiDDecayTree->Branch("T2DeuteronPz",&T2DeuteronPz,"T2DeuteronPz/F");
    PiDDecayTree->Branch("T2DeuteronDCAtoPVXY",&T2DeuteronDCAtoPVXY,"T2DeuteronDCAtoPVXY/F");
    PiDDecayTree->Branch("T2DeuteronDCAtoPVZ",&T2DeuteronDCAtoPVZ,"T2DeuteronDCAtoPVZ/F");
    PiDDecayTree->Branch("T2DeuteronDCAtoSVXY",&T2DeuteronDCAtoSVXY,"T2DeuteronDCAtoSVXY/F");
    PiDDecayTree->Branch("T2DeuteronTPCNcls",&T2DeuteronTPCNcls,"T2DeuteronTPCNcls/F");
    PiDDecayTree->Branch("T2DeuteronTPCnSigma",&T2DeuteronTPCnSigma,"T2DeuteronTPCnSigma/F");
    PiDDecayTree->Branch("T2PionPx",&T2PionPx,"T2PionPx/F");
    PiDDecayTree->Branch("T2PionPy",&T2PionPy,"T2PionPy/F");
    PiDDecayTree->Branch("T2PionPz",&T2PionPz,"T2PionPz/F");
    PiDDecayTree->Branch("T2PionDCAtoPVXY",&T2PionDCAtoPVXY,"T2PionDCAtoPVXY/F");
    PiDDecayTree->Branch("T2PionDCAtoPVZ",&T2PionDCAtoPVZ,"T2PionDCAtoPVZ/F");
    PiDDecayTree->Branch("T2PionDCAtoSVXY",&T2PionDCAtoSVXY,"T2PionDCAtoSVXY/F");
    PiDDecayTree->Branch("T2PionTPCNcls",&T2PionTPCNcls,"T2PionTPCNcls/F");
    PiDDecayTree->Branch("T2PionTPCnSigma",&T2PionTPCnSigma,"T2PionTPCnSigma/F");
  }
  PostData(1,fOutputList);
  PostData(2,PiKDCascadeDecayTree);
  PostData(3,PiDDecayTree);
}
void AliAnalysisTaskHyperFinder3Body::UserExec(Option_t *)
{
  fhEventCount->Fill(0.5);
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { AliError(Form("%s: Could not get Analysis Manager", GetName())); return;}
  fInputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
  if (!fInputHandler) { AliError(Form("%s: Could not get Input Handler", GetName())); return;}
  fPIDResponse=fInputHandler->GetPIDResponse(); 
  if (!fPIDResponse) { AliError(Form("%s: Could not get PIDResponse", GetName())); return;}
  fESDevent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESDevent) { Printf("ERROR: Could not retrieve event"); return;}
  fhEventCount->Fill(1.5);
  AliMultSelection *fMultSelection = (AliMultSelection *) fESDevent->FindListObject("MultSelection");
  Float_t percentile;
  percentile = fMultSelection->GetMultiplicityPercentile("V0M");
  fhcentralitybeforecut->Fill(percentile);
  Int_t frunnumber = fESDevent->GetRunNumber();
  if (frunnumber == 297219 || frunnumber == 297194 || frunnumber == 297029
		  || frunnumber == 296890 || frunnumber == 296849 || frunnumber == 296750
		  || frunnumber == 296749 || frunnumber == 297481) feventCut.UseTimeRangeCut();
  feventCut.fRequireTrackVertex = kTRUE;
  feventCut.OverrideAutomaticTriggerSelection(AliVEvent::kAny);
  if(!feventCut.AcceptEvent(fESDevent))return;
  fhEventCount->Fill(2.5);
  fhcentralityaftercut->Fill(percentile);
  Bool_t MB = kFALSE;
  Bool_t HMV0 = kFALSE;
  Bool_t HMSPD = kFALSE;
  Bool_t HNU = kFALSE;
  Bool_t HQU = kFALSE;
  Bool_t Central = kFALSE;
  Bool_t SemiCentral = kFALSE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)        MB = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)  HMV0 = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD) HMSPD = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)     Central = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral) SemiCentral = kTRUE;
  TString classes = fESDevent->GetFiredTriggerClasses();
  if (classes.Contains("HNU")) HNU = kTRUE;
  if (classes.Contains("HQU")) HQU = kTRUE;
  fhTrigger->Fill(0);
  if (MB)          fhTrigger->Fill(1);
  if (HMV0)        fhTrigger->Fill(2);
  if (HMSPD)       fhTrigger->Fill(3);
  if (HNU)         fhTrigger->Fill(4);
  if (HQU)         fhTrigger->Fill(5);
  if (Central)     fhTrigger->Fill(6);
  if (SemiCentral) fhTrigger->Fill(7);
  Double_t PVposition[3];
  const AliESDVertex *vtTRC = fESDevent->GetPrimaryVertex();
  const AliVVertex *primaryVertex = InputEvent()->GetPrimaryVertex();
  vtTRC->GetXYZ(PVposition);
  Double_t lMagneticField = fESDevent->GetMagneticField();
  Int_t tracknumber;
  tracknumber = fESDevent->GetNumberOfTracks();
  vector<Int_t> Pionid;
  vector<Int_t> Deuteronid;
  vector<Int_t> Kaonid;
  for (Int_t i=0; i<tracknumber; i++)
  {
    AliESDtrack* track = 0x0;
    track = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(i));
    if(!track) continue;
    ULong_t status = (ULong_t)track->GetStatus();
    Bool_t isTPC = (((status) & AliVTrack::kTPCin) !=0);
    if(!isTPC) continue;
    Double_t pinTPC = track->GetTPCmomentum();
    Double_t d = track->GetD(PVposition[0],PVposition[1],lMagneticField);
    Float_t Trackimpactxy;
    Float_t Trackimpactz;
    track->GetImpactParameters(Trackimpactxy,Trackimpactz);
    fhTrackDCAtoPVxycheck->Fill(Trackimpactxy-d);
    Double_t TrackLength = track->GetLengthInActiveZone(1, 2.0, 220.0, lMagneticField);
    Double_t TrackNCrossedRows = track->GetTPCClusterInfo(2,1);
    Double_t TPCSignal = track->GetTPCsignal();
    Double_t nSigmaPion = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
    Double_t nSigmaPionTOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
    Double_t nSigmaDeuteron = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kDeuteron);
    Double_t nSigmaDeuteronTOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kDeuteron);
    Double_t nSigmaKaon = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
    Double_t nSigmaKaonTOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
    Bool_t trackflag = kFALSE;
    Bool_t trackflagPion = kTRUE;
    trackflag = fESDtrackCuts->AcceptTrack(track);
    //if(d <= TrackDCAtoPVcut || TrackLength <= TrackLengthcut || TrackNCrossedRows <= TrackNCrossedRowscut)trackflag = kFALSE;
    if(TMath::Abs(nSigmaPion) < 4 && d >= TrackDCAtoPVcut && trackflag && track->Pt() < 2)
    {
      //if(pinTPC > 0.5 && TMath::Abs(nSigmaPionTOF) >= 3)trackflagPion = kFALSE;
      if(trackflagPion)
      {
        fhBBPion->Fill(pinTPC*track->Charge(),TPCSignal);
        Pionid.push_back(i);
      }
    }
    if(TMath::Abs(nSigmaDeuteron) < 4 && d >= TrackDCAtoPVcut && trackflag)
    {
      fhBBDeuteron->Fill(pinTPC*track->Charge(),TPCSignal);
      Deuteronid.push_back(i);
    }
    if(TMath::Abs(nSigmaKaon) < 4 /*&& d >= TrackDCAtoPVcut*/ && trackflag)
    {
      fhBBKaon->Fill(pinTPC*track->Charge(),TPCSignal);
      Kaonid.push_back(i);
    }
  }
  if(Deuteronid.size() == 0) {fhDeuteronEventCheck->Fill(1.5); return;}
  else {fhDeuteronEventCheck->Fill(0.5);}
  for(Int_t i=0; i<Pionid.size(); i++)
  {
    Int_t Pionidx = Pionid[i];
    AliESDtrack* PionTrack = 0x0;
    PionTrack = (fESDevent->GetTrack(Pionidx));
    AliExternalTrackParam trackInPion(*PionTrack);
    Float_t PionImpactXY =0;
    Float_t PionImpactZ =0;
    // rotate PionTrack
    // rotate over
    PionTrack->GetImpactParameters(PionImpactXY,PionImpactZ);
    if(PionTrack->Charge() < 0)
    {
      for(Int_t j=0; j<Deuteronid.size(); j++)
      {
        Int_t Deuteronidx = Deuteronid[j];
        if(Pionidx == Deuteronidx)continue;
        AliESDtrack* DeuteronTrack = 0x0;
        DeuteronTrack = (fESDevent->GetTrack(Deuteronidx));
        if(DeuteronTrack->Charge() < 0)continue;
        AliExternalTrackParam trackInDeuteron(*DeuteronTrack);
        //rotate DeuteronTrack
        //rotate over
        Float_t DeuteronImpactXY;
        Float_t DeuteronImpactZ;
        DeuteronTrack->GetImpactParameters(DeuteronImpactXY,DeuteronImpactZ);
        Double_t xn, xp;
        Double_t Dcadaughters = 0;
        Dcadaughters = PionTrack->GetDCA(DeuteronTrack,lMagneticField,xn,xp);
        if((xn+xp) < 2*V0Radiusmin || (xn+xp) > 2*V0Radiusmax)continue;
        Bool_t corrected=kFALSE;
        if((trackInPion.GetX() > 3.) && (xn < 3.)) {corrected=kTRUE;}
        if((trackInDeuteron.GetX() > 3.) && (xp < 3.)) {corrected=kTRUE;}
        if(corrected && PionTrack->Charge()<0) {Dcadaughters = trackInPion.GetDCA(&trackInDeuteron,lMagneticField,xn,xp);}
        if(Dcadaughters > V0DCAdaughtercut)continue;
        trackInPion.PropagateTo(xn,lMagneticField);
        trackInDeuteron.PropagateTo(xp,lMagneticField);
        AliESDv0 vertex(trackInPion,Pionidx,trackInDeuteron,Deuteronidx);
        Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
        if(SVtoPV < V0DCAtoPVcut)continue;
        Double_t SVposition[3] = {0., 0., 0.};
        SVposition[0] = vertex.Xv();
        SVposition[1] = vertex.Yv();
        SVposition[2] = vertex.Zv();
        Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
        if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax)continue;
        Double_t pPionAt[3] = {0., 0., 0.};
        Double_t pDeuteronAt[3] = {0., 0., 0.};
        PionTrack->GetPxPyPzAt(xn,lMagneticField,pPionAt);
        DeuteronTrack->GetPxPyPzAt(xp,lMagneticField,pDeuteronAt);
        Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
        if(CPAV0 < CascadeV0CPAcut)continue;
        TLorentzVector vDeuteron,vPion,vsumV0;
        vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
        vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
        vsumV0 = vDeuteron + vPion;
        Double_t V0Invmass = vsumV0.M();
        //--------leviate the cut----------
        if(V0Invmass > CascadeV0Invmassmax)continue;
        if(DeuteronTrack->Charge()>0)
        {
          for(Int_t k=0; k<Kaonid.size(); k++)
          {
            Int_t Kaonidx = Kaonid[k];
            AliESDtrack* KaonTrack = 0x0;
            KaonTrack = (fESDevent->GetTrack(Kaonidx));
            AliExternalTrackParam trackInKaon(*KaonTrack);
            //rotate KaonTrack
            //rotate over
            if(KaonTrack->Charge() > 0)continue;
            if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
            fhPiKDChargeCheckC1->Fill(0.5,PionTrack->Charge());
            fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrack->Charge());
            fhPiKDChargeCheckC1->Fill(2.5,KaonTrack->Charge());
            Float_t KaonImpactXY;
            Float_t KaonImpactZ;
            KaonTrack->GetImpactParameters(KaonImpactXY,KaonImpactZ);
            Double_t xkaon = 0;
            Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&trackInKaon,lMagneticField,&xkaon);
            if(DCACascade > CascadeDCAcut)continue;
            AliESDcascade cascade(vertex,trackInKaon,Kaonidx);
            Double_t Cascadeposition[3] = {0., 0., 0.};
            cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
            Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
            if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
            Double_t pV0[3] = {0., 0., 0.};
            vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
            fhCausalitycheck->Fill(0.5);
            if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
            fhCausalitycheck->Fill(1.5);
            if(r2cascade > r2)continue;
            Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
            if(CascadeCPA < CascadeCPAcut)continue;
            Double_t pCascade[3] = {0., 0., 0.};
            cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
            Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
            if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
            Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
            Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
            Double_t pKaonAt[3] = {0., 0., 0.};
            KaonTrack->GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
            Double_t pKaon[3] = {0., 0., 0.};
            pKaon[0] = trackInKaon.Px();
            pKaon[1] = trackInKaon.Py();
            pKaon[2] = trackInKaon.Pz();
            fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
            fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
            fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
            TLorentzVector vKaon, vsumCascade;
            vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
            vsumCascade = vsumV0 + vKaon;
            Float_t InvmassCascade =  (Float_t)vsumCascade.M();
            //----------leviate the cut----------
            if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
            T1DecayChannel = (Int_t)1;
            T1EventPercentile = (Float_t)percentile;
            T1TVDCA = (Float_t)DCACascade;
            T1TVCPA = (Float_t)CascadeCPA;
            T1TVtoPV = (Float_t)DCACascadetoPV;
            T1SVDCA = (Float_t)Dcadaughters;
            T1SVCPA = (Float_t)CPAV0;
            T1SVtoPV = (Float_t)SVtoPV;
            T1SVtoTV = (Float_t)DCACascadetoSV;
            T1KaonPx = (Float_t)pKaonAt[0];
            T1KaonPy = (Float_t)pKaonAt[1];
            T1KaonPz = (Float_t)pKaonAt[2];
            T1KaonDCAtoPVXY = (Float_t)KaonImpactXY;
            T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
            T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
            T1DeuteronPx = (Float_t)pDeuteronAt[0];
            T1DeuteronPy = (Float_t)pDeuteronAt[1];
            T1DeuteronPz = (Float_t)pDeuteronAt[2];
            T1DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
            T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
            T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
            T1PionPx = (Float_t)pPionAt[0];
            T1PionPy = (Float_t)pPionAt[1];
            T1PionPz = (Float_t)pPionAt[2];
            T1PionDCAtoPVXY = (Float_t)PionImpactXY;
            T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
            T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
            PiKDCascadeDecayTree->Fill();
          }
        }
        else continue;
      }
    }
    else
    {
      for(Int_t j=0; j<Deuteronid.size(); j++)
      {
        Int_t Deuteronidx = Deuteronid[j];
        if(Pionidx == Deuteronidx)continue;
        AliESDtrack* DeuteronTrack = 0x0;
        DeuteronTrack = (fESDevent->GetTrack(Deuteronidx));
        AliExternalTrackParam trackInDeuteron(*DeuteronTrack);
        //rotate Deuterontrack
        //rotate over
        Float_t DeuteronImpactXY;
        Float_t DeuteronImpactZ;
        DeuteronTrack->GetImpactParameters(DeuteronImpactXY,DeuteronImpactZ);
        Double_t xn, xp;
        Double_t Dcadaughters = 0;
        Dcadaughters = DeuteronTrack->GetDCA(PionTrack,lMagneticField,xn,xp);
        if((xn+xp) < 2*V0Radiusmin || (xn+xp) > 2*V0Radiusmax)continue;
        Bool_t corrected=kFALSE;
        if((trackInPion.GetX() > 3.) && (xp < 3.)) {corrected=kTRUE;}
        if((trackInDeuteron.GetX() > 3.) && (xn < 3.)) {corrected=kTRUE;}
        if(corrected && PionTrack->Charge()>0) {Dcadaughters = trackInDeuteron.GetDCA(&trackInPion,lMagneticField,xn,xp);}
        if(Dcadaughters > V0DCAdaughtercut)continue;
        trackInPion.PropagateTo(xp,lMagneticField);
        trackInDeuteron.PropagateTo(xn,lMagneticField);
        AliESDv0 vertex(trackInDeuteron,Deuteronidx,trackInPion,Pionidx);
        Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
        if(SVtoPV < V0DCAtoPVcut)continue;
        Double_t SVposition[3] = {0., 0., 0.};
        SVposition[0] = vertex.Xv();
        SVposition[1] = vertex.Yv();
        SVposition[2] = vertex.Zv();
        Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
        if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax)continue;
        Double_t pPionAt[3] = {0., 0., 0.};
        Double_t pDeuteronAt[3] = {0., 0., 0.};
        PionTrack->GetPxPyPzAt(xp,lMagneticField,pPionAt);
        DeuteronTrack->GetPxPyPzAt(xn,lMagneticField,pDeuteronAt);
        Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
        if(CPAV0 < CascadeV0CPAcut)continue;
        TLorentzVector vDeuteron,vPion,vsumV0;
        vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
        vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
        vsumV0 = vDeuteron + vPion;
        Double_t V0Invmass = vsumV0.M();
        //-----------leviate the cut--------------
        if(V0Invmass > CascadeV0Invmassmax)continue;
        if(DeuteronTrack->Charge()>0)continue;
        else if(DeuteronTrack->Charge()<0)
        {
          for(Int_t k=0; k<Kaonid.size(); k++)
          {
            Int_t Kaonidx = Kaonid[k];
            AliESDtrack* KaonTrack = 0x0;
            KaonTrack = (fESDevent->GetTrack(Kaonidx));
            AliExternalTrackParam trackInKaon(*KaonTrack);
            //rotate KaonTrack
            //rotate over
            if(KaonTrack->Charge() < 0)continue;
            if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
            fhPiKDChargeCheckC2->Fill(0.5,PionTrack->Charge());
            fhPiKDChargeCheckC2->Fill(1.5,DeuteronTrack->Charge());
            fhPiKDChargeCheckC2->Fill(2.5,KaonTrack->Charge());
            Float_t KaonImpactXY;
            Float_t KaonImpactZ;
            KaonTrack->GetImpactParameters(KaonImpactXY,KaonImpactZ);
            Double_t xkaon = 0;
            Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&trackInKaon,lMagneticField,&xkaon);
            if(DCACascade > CascadeDCAcut)continue;
            AliESDcascade cascade(vertex,trackInKaon,Kaonidx);
            Double_t Cascadeposition[3] = {0., 0., 0.};
            cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
            Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
            if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
            Double_t pV0[3] = {0., 0., 0.};
            vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
            fhCausalitycheck->Fill(0.5);
            if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
            fhCausalitycheck->Fill(1.5);
            if(r2cascade > r2)continue;
            Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
            if(CascadeCPA < CascadeCPAcut)continue;
            Double_t pCascade[3] = {0., 0., 0.};
            cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
            Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
            if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
            Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
            Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
            Double_t pKaonAt[3] = {0., 0., 0.};
            KaonTrack->GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
            TLorentzVector vKaon, vsumCascade;
            vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
            vsumCascade = vsumV0 + vKaon;
            Float_t InvmassCascade =  (Float_t)vsumCascade.M();
            if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
            T1DecayChannel = (Int_t)2;
            T1EventPercentile = (Float_t)percentile;
            T1TVDCA = (Float_t)DCACascade;
            T1TVCPA = (Float_t)CascadeCPA;
            T1TVtoPV = (Float_t)DCACascadetoPV;
            T1SVDCA = (Float_t)Dcadaughters;
            T1SVCPA = (Float_t)CPAV0;
            T1SVtoPV = (Float_t)SVtoPV;
            T1SVtoTV = (Float_t)DCACascadetoSV;
            T1KaonPx = (Float_t)pKaonAt[0];
            T1KaonPy = (Float_t)pKaonAt[1];
            T1KaonPz = (Float_t)pKaonAt[2];
            T1KaonDCAtoPVXY = (Float_t)KaonImpactXY;
            T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
            T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
            T1DeuteronPx = (Float_t)pDeuteronAt[0];
            T1DeuteronPy = (Float_t)pDeuteronAt[1];
            T1DeuteronPz = (Float_t)pDeuteronAt[2];
            T1DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
            T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
            T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
            T1PionPx = (Float_t)pPionAt[0];
            T1PionPy = (Float_t)pPionAt[1];
            T1PionPz = (Float_t)pPionAt[2];
            T1PionDCAtoPVXY = (Float_t)PionImpactXY;
            T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
            T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
            PiKDCascadeDecayTree->Fill();
          }
        }
      }
    }
  }
//=====================================================================================================================================//
//===================================================Rotation Track Loop===============================================================//
//=====================================================================================================================================//
  for(Int_t i=0; i<Pionid.size(); i++)
  {
    Int_t Pionidx = Pionid[i];
    AliESDtrack* PionTrack = 0x0;
    PionTrack = (fESDevent->GetTrack(Pionidx));
    AliExternalTrackParam trackInPion(*PionTrack);
    if(kRotatePion)
    {
      Double_t pseudoXPion[3], pseudoPPion[3], CovPseudoPion[21];
      Double_t pPionRot[3];
      PionTrack->GetCovarianceXYZPxPyPz(CovPseudoPion);
      PionTrack->GetXYZ(pseudoXPion);
      PionTrack->GetPxPyPz(pseudoPPion);
      Short_t signPion = PionTrack->Charge();
      Double_t fRot = nRotation;
      Double_t fAngle = 0.0872; // 6 gradi
      if(nRotation == 20)
      {
        fAngle = 0.1047;// 5 gradi
      }
      for(Int_t r =0; r < fRot ; r++)
      {
        AliExternalTrackParam PionTrackRotated;
        if(nRotation == 20)
        {
          pPionRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoPPion[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoPPion[1];
          pPionRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoPPion[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoPPion[1];
          pPionRot[2] = pseudoPPion[2];
          PionTrackRotated = AliExternalTrackParam(pseudoXPion,pPionRot,CovPseudoPion,signPion);
        }
        else
        {
          pPionRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoPPion[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoPPion[1];
          pPionRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoPPion[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoPPion[1];
          pPionRot[2] = pseudoPPion[2];
          PionTrackRotated = AliExternalTrackParam(pseudoXPion,pPionRot,CovPseudoPion,signPion);
        }
        Double_t Piond0Rot[2], Pioncovd0Rot[3];
        PionTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Piond0Rot, Pioncovd0Rot);
        if(PionTrackRotated.Charge() < 0)
        {
          for(Int_t j=0; j<Deuteronid.size(); j++)
          {
            Int_t Deuteronidx = Deuteronid[j];
            if(Pionidx == Deuteronidx)continue;
            AliESDtrack* DeuteronTrack = 0x0;
            DeuteronTrack = (fESDevent->GetTrack(Deuteronidx));
            if(DeuteronTrack->Charge() < 0)continue;
            AliExternalTrackParam trackInDeuteron(*DeuteronTrack);
            if(kRotateDeuteron)
            {
              Double_t pseudoXDeuteron[3], pseudoPDeuteron[3], CovPseudoDeuteron[21];
              Double_t pDeuteronRot[3];
              DeuteronTrack->GetCovarianceXYZPxPyPz(CovPseudoDeuteron);
              DeuteronTrack->GetXYZ(pseudoXDeuteron);
              DeuteronTrack->GetPxPyPz(pseudoPDeuteron);
              Short_t signDeuteron = DeuteronTrack->Charge();
              for(Int_t r_Deuteron=0; r_Deuteron<fRot; r_Deuteron++)
              {
                AliExternalTrackParam DeuteronTrackRotated;
                if(nRotation == 20)
                {
                  pDeuteronRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoPDeuteron[1];
                  pDeuteronRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoPDeuteron[1];
                  pDeuteronRot[2] = pseudoPDeuteron[2];
                  DeuteronTrackRotated = AliExternalTrackParam(pseudoXDeuteron,pDeuteronRot,CovPseudoDeuteron,signDeuteron);
                }
                else
                {
                  pDeuteronRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoPDeuteron[1];
                  pDeuteronRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoPDeuteron[1];
                  pDeuteronRot[2] = pseudoPDeuteron[2];
                  DeuteronTrackRotated = AliExternalTrackParam(pseudoXDeuteron,pDeuteronRot,CovPseudoDeuteron,signDeuteron);
                }
                Double_t Deuterond0Rot[2], Deuteroncovd0Rot[3];
                DeuteronTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Deuterond0Rot, Deuteroncovd0Rot);
                Double_t xn, xp;
                Double_t Dcadaughters = 0;
                Dcadaughters = PionTrackRotated.GetDCA(&DeuteronTrackRotated,lMagneticField,xn,xp);
                if(Dcadaughters > V0DCAdaughtercut)continue;
                PionTrackRotated.PropagateTo(xn,lMagneticField);
                DeuteronTrackRotated.PropagateTo(xp,lMagneticField);
                AliESDv0 vertex(PionTrackRotated,Pionidx,DeuteronTrackRotated,Deuteronidx);
                Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
                if(SVtoPV < V0DCAtoPVcut)continue;
                Double_t SVposition[3] = {0., 0., 0.};
                SVposition[0] = vertex.Xv();
                SVposition[1] = vertex.Yv();
                SVposition[2] = vertex.Zv();
                Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
                if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax)continue;
                Double_t pPionAt[3] = {0., 0., 0.};
                Double_t pDeuteronAt[3] = {0., 0., 0.};
                PionTrackRotated.GetPxPyPzAt(xn,lMagneticField,pPionAt);
                DeuteronTrackRotated.GetPxPyPzAt(xp,lMagneticField,pDeuteronAt);
                Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                if(CPAV0 < CascadeV0CPAcut)continue;
                TLorentzVector vDeuteron,vPion,vsumV0;
                vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
                vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
                vsumV0 = vDeuteron + vPion;
                Double_t V0Invmass = vsumV0.M();
                if(V0Invmass > CascadeV0Invmassmax)continue;
                if(DeuteronTrackRotated.Charge() > 0)
                {
                  for(Int_t k=0; k<Kaonid.size(); k++)
                  {
                    Int_t Kaonidx = Kaonid[k];
                    AliESDtrack* KaonTrack = 0x0;
                    KaonTrack = (fESDevent->GetTrack(Kaonidx));
                    AliExternalTrackParam trackInKaon(*KaonTrack);
                    if(kRotateKaon)
                    {
                      Double_t pseudoXKaon[3], pseudoPKaon[3], CovPseudoKaon[21];
                      Double_t pKaonRot[3];
                      KaonTrack->GetCovarianceXYZPxPyPz(CovPseudoKaon);
                      KaonTrack->GetXYZ(pseudoXKaon);
                      KaonTrack->GetPxPyPz(pseudoPKaon);
                      Short_t signKaon = KaonTrack->Charge();
                      for(Int_t r_Kaon=0; r_Kaon<fRot; r_Kaon++)
                      {
                        AliExternalTrackParam KaonTrackRotated;
                        if(nRotation == 20)
                        {
                          pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                          pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                          pKaonRot[2] = pseudoPKaon[2];
                          KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                        }
                        else
                        {
                          pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                          pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                          pKaonRot[2] = pseudoPKaon[2];
                          KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                        }
                        if(KaonTrackRotated.Charge() > 0)continue;
                        if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                        fhPiKDChargeCheckC1->Fill(0.5,PionTrackRotated.Charge());
                        fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrackRotated.Charge());
                        fhPiKDChargeCheckC1->Fill(2.5,KaonTrackRotated.Charge());
                        Double_t Kaond0Rot[2], Kaoncovd0Rot[3];
                        KaonTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Kaond0Rot, Kaoncovd0Rot);
                        Double_t xkaon = 0;
                        Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&KaonTrackRotated,lMagneticField,&xkaon);
                        if(DCACascade > CascadeDCAcut)continue;
                        AliESDcascade cascade(vertex,KaonTrackRotated,Kaonidx);
                        Double_t Cascadeposition[3] = {0., 0., 0.};
                        cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                        Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                        if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                        Double_t pV0[3] = {0., 0., 0.};
                        vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                        fhCausalitycheck->Fill(0.5);
                        if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                        fhCausalitycheck->Fill(1.5);
                        if(r2cascade > r2)continue;
                        Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                        if(CascadeCPA < CascadeCPAcut)continue;
                        Double_t pCascade[3] = {0., 0., 0.};
                        cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                        Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                        if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                        Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                        Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                        Double_t pKaonAt[3] = {0., 0., 0.};
                        KaonTrackRotated.GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                        Double_t pKaon[3] = {0., 0., 0.};
                        pKaon[0] = KaonTrackRotated.Px();
                        pKaon[1] = KaonTrackRotated.Py();
                        pKaon[2] = KaonTrackRotated.Pz();
                        fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                        fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                        fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                        TLorentzVector vKaon, vsumCascade;
                        vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                        vsumCascade = vsumV0 + vKaon;
                        Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                        if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                        T1DecayChannel = (Int_t)3;
                        T1EventPercentile = (Float_t)percentile;
                        T1TVDCA = (Float_t)DCACascade;
                        T1TVCPA = (Float_t)CascadeCPA;
                        T1TVtoPV = (Float_t)DCACascadetoPV;
                        T1SVDCA = (Float_t)Dcadaughters;
                        T1SVCPA = (Float_t)CPAV0;
                        T1SVtoPV = (Float_t)SVtoPV;
                        T1SVtoTV = (Float_t)DCACascadetoSV;
                        T1KaonPx = (Float_t)pKaonAt[0];
                        T1KaonPy = (Float_t)pKaonAt[1];
                        T1KaonPz = (Float_t)pKaonAt[2];
                        T1KaonDCAtoPVXY = (Float_t)Kaond0Rot[0];
                        T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                        T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                        T1DeuteronPx = (Float_t)pDeuteronAt[0];
                        T1DeuteronPy = (Float_t)pDeuteronAt[1];
                        T1DeuteronPz = (Float_t)pDeuteronAt[2];
                        T1DeuteronDCAtoPVXY = (Float_t)Deuterond0Rot[0];
                        T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                        T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                        T1PionPx = (Float_t)pPionAt[0];
                        T1PionPy = (Float_t)pPionAt[1];
                        T1PionPz = (Float_t)pPionAt[2];
                        T1PionDCAtoPVXY = (Float_t)Piond0Rot[0];
                        T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                        T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                        PiKDCascadeDecayTree->Fill();
                      }
                    }
                    else
                    {
                      if(KaonTrack->Charge() > 0)continue;
                      if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                      fhPiKDChargeCheckC1->Fill(0.5,PionTrackRotated.Charge());
                      fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrackRotated.Charge());
                      fhPiKDChargeCheckC1->Fill(2.5,KaonTrack->Charge());
                      Float_t KaonImpactXY;
                      Float_t KaonImpactZ;
                      KaonTrack->GetImpactParameters(KaonImpactXY,KaonImpactZ);
                      Double_t xkaon = 0;
                      Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&trackInKaon,lMagneticField,&xkaon);
                      if(DCACascade > CascadeDCAcut)continue;
                      AliESDcascade cascade(vertex,trackInKaon,Kaonidx);
                      Double_t Cascadeposition[3] = {0., 0., 0.};
                      cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                      Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                      if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                      Double_t pV0[3] = {0., 0., 0.};
                      vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                      fhCausalitycheck->Fill(0.5);
                      if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                      fhCausalitycheck->Fill(1.5);
                      if(r2cascade > r2)continue;
                      Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                      if(CascadeCPA < CascadeCPAcut)continue;
                      Double_t pCascade[3] = {0., 0., 0.};
                      cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                      Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                      if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                      Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                      Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                      Double_t pKaonAt[3] = {0., 0., 0.};
                      KaonTrack->GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                      Double_t pKaon[3] = {0., 0., 0.};
                      pKaon[0] = trackInKaon.Px();
                      pKaon[1] = trackInKaon.Py();
                      pKaon[2] = trackInKaon.Pz();
                      fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                      fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                      fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                      TLorentzVector vKaon, vsumCascade;
                      vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                      vsumCascade = vsumV0 + vKaon;
                      Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                      if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                      T1DecayChannel = (Int_t)3;
                      T1EventPercentile = (Float_t)percentile;
                      T1TVDCA = (Float_t)DCACascade;
                      T1TVCPA = (Float_t)CascadeCPA;
                      T1TVtoPV = (Float_t)DCACascadetoPV;
                      T1SVDCA = (Float_t)Dcadaughters;
                      T1SVCPA = (Float_t)CPAV0;
                      T1SVtoPV = (Float_t)SVtoPV;
                      T1SVtoTV = (Float_t)DCACascadetoSV;
                      T1KaonPx = (Float_t)pKaonAt[0];
                      T1KaonPy = (Float_t)pKaonAt[1];
                      T1KaonPz = (Float_t)pKaonAt[2];
                      T1KaonDCAtoPVXY = (Float_t)KaonImpactXY;
                      T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                      T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                      T1DeuteronPx = (Float_t)pDeuteronAt[0];
                      T1DeuteronPy = (Float_t)pDeuteronAt[1];
                      T1DeuteronPz = (Float_t)pDeuteronAt[2];
                      T1DeuteronDCAtoPVXY = (Float_t)Deuterond0Rot[0];
                      T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                      T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                      T1PionPx = (Float_t)pPionAt[0];
                      T1PionPy = (Float_t)pPionAt[1];
                      T1PionPz = (Float_t)pPionAt[2];
                      T1PionDCAtoPVXY = (Float_t)Piond0Rot[0];
                      T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                      T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                      PiKDCascadeDecayTree->Fill();
                    }
                  }
                }
              }
            }
            else
            {
              Float_t DeuteronImpactXY;
              Float_t DeuteronImpactZ;
              DeuteronTrack->GetImpactParameters(DeuteronImpactXY,DeuteronImpactZ);
              Double_t xn, xp;
              Double_t Dcadaughters = 0;
              Dcadaughters = PionTrackRotated.GetDCA(&trackInDeuteron,lMagneticField,xn,xp);
              if((xn+xp) < 2*V0Radiusmin || (xn+xp) > 2*V0Radiusmax)continue;
              if(Dcadaughters > V0DCAdaughtercut)continue;
              PionTrackRotated.PropagateTo(xn,lMagneticField);
              trackInDeuteron.PropagateTo(xp,lMagneticField);
              AliESDv0 vertex(PionTrackRotated,Pionidx,trackInDeuteron,Deuteronidx);
              Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
              if(SVtoPV < V0DCAtoPVcut)continue;
              Double_t SVposition[3] = {0., 0., 0.};
              SVposition[0] = vertex.Xv();
              SVposition[1] = vertex.Yv();
              SVposition[2] = vertex.Zv();
              Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
              if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax)continue;
              Double_t pPionAt[3] = {0., 0., 0.};
              Double_t pDeuteronAt[3] = {0., 0., 0.};
              PionTrackRotated.GetPxPyPzAt(xn,lMagneticField,pPionAt);
              DeuteronTrack->GetPxPyPzAt(xp,lMagneticField,pDeuteronAt);
              Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
              if(CPAV0 < CascadeV0CPAcut)continue;
              TLorentzVector vDeuteron,vPion,vsumV0;
              vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
              vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
              vsumV0 = vDeuteron + vPion;
              Double_t V0Invmass = vsumV0.M();
              if(V0Invmass > CascadeV0Invmassmax)continue;
              if(DeuteronTrack->Charge()>0)
              {
                for(Int_t k=0; k<Kaonid.size(); k++)
                {
                  Int_t Kaonidx = Kaonid[k];
                  AliESDtrack* KaonTrack = 0x0;
                  KaonTrack = (fESDevent->GetTrack(Kaonidx));
                  AliExternalTrackParam trackInKaon(*KaonTrack);
                  if(KaonTrack->Charge() > 0)continue;
                  if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                  fhPiKDChargeCheckC1->Fill(0.5,PionTrackRotated.Charge());
                  fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrack->Charge());
                  fhPiKDChargeCheckC1->Fill(2.5,KaonTrack->Charge());
                  if(kRotateKaon)
                  {
                    Double_t pseudoXKaon[3], pseudoPKaon[3], CovPseudoKaon[21];
                    Double_t pKaonRot[3];
                    KaonTrack->GetCovarianceXYZPxPyPz(CovPseudoKaon);
                    KaonTrack->GetXYZ(pseudoXKaon);
                    KaonTrack->GetPxPyPz(pseudoPKaon);
                    Short_t signKaon = KaonTrack->Charge();
                    for(Int_t r_Kaon=0; r_Kaon<fRot; r_Kaon++)
                    {
                      AliExternalTrackParam KaonTrackRotated;
                      if(nRotation == 20)
                      {
                        pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[2] = pseudoPKaon[2];
                        KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                      }
                      else
                      {
                        pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[2] = pseudoPKaon[2];
                        KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                      }
                      if(KaonTrackRotated.Charge() > 0)continue;
                      if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                      fhPiKDChargeCheckC1->Fill(0.5,PionTrackRotated.Charge());
                      fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrack->Charge());
                      fhPiKDChargeCheckC1->Fill(2.5,KaonTrackRotated.Charge());
                      Double_t Kaond0Rot[2], Kaoncovd0Rot[3];
                      KaonTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Kaond0Rot, Kaoncovd0Rot);
                      Double_t xkaon = 0;
                      Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&KaonTrackRotated,lMagneticField,&xkaon);
                      if(DCACascade > CascadeDCAcut)continue;
                      AliESDcascade cascade(vertex,KaonTrackRotated,Kaonidx);
                      Double_t Cascadeposition[3] = {0., 0., 0.};
                      cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                      Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                      if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                      Double_t pV0[3] = {0., 0., 0.};
                      vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                      fhCausalitycheck->Fill(0.5);
                      if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                      fhCausalitycheck->Fill(1.5);
                      if(r2cascade > r2)continue;
                      Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                      if(CascadeCPA < CascadeCPAcut)continue;
                      Double_t pCascade[3] = {0., 0., 0.};
                      cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                      Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                      if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                      Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                      Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                      Double_t pKaonAt[3] = {0., 0., 0.};
                      KaonTrackRotated.GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                      Double_t pKaon[3] = {0., 0., 0.};
                      pKaon[0] = KaonTrackRotated.Px();
                      pKaon[1] = KaonTrackRotated.Py();
                      pKaon[2] = KaonTrackRotated.Pz();
                      fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                      fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                      fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                      TLorentzVector vKaon, vsumCascade;
                      vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                      vsumCascade = vsumV0 + vKaon;
                      Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                      if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                      T1DecayChannel = (Int_t)3;
                      T1EventPercentile = (Float_t)percentile;
                      T1TVDCA = (Float_t)DCACascade;
                      T1TVCPA = (Float_t)CascadeCPA;
                      T1TVtoPV = (Float_t)DCACascadetoPV;
                      T1SVDCA = (Float_t)Dcadaughters;
                      T1SVCPA = (Float_t)CPAV0;
                      T1SVtoPV = (Float_t)SVtoPV;
                      T1SVtoTV = (Float_t)DCACascadetoSV;
                      T1KaonPx = (Float_t)pKaonAt[0];
                      T1KaonPy = (Float_t)pKaonAt[1];
                      T1KaonPz = (Float_t)pKaonAt[2];
                      T1KaonDCAtoPVXY = (Float_t)Kaond0Rot[0];
                      T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                      T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                      T1DeuteronPx = (Float_t)pDeuteronAt[0];
                      T1DeuteronPy = (Float_t)pDeuteronAt[1];
                      T1DeuteronPz = (Float_t)pDeuteronAt[2];
                      T1DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
                      T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                      T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                      T1PionPx = (Float_t)pPionAt[0];
                      T1PionPy = (Float_t)pPionAt[1];
                      T1PionPz = (Float_t)pPionAt[2];
                      T1PionDCAtoPVXY = (Float_t)Piond0Rot[0];
                      T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                      T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                      PiKDCascadeDecayTree->Fill();
                    }
                  }
                  else
                  {
                    Float_t KaonImpactXY;
                    Float_t KaonImpactZ;
                    KaonTrack->GetImpactParameters(KaonImpactXY,KaonImpactZ);
                    Double_t xkaon = 0;
                    Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&trackInKaon,lMagneticField,&xkaon);
                    if(DCACascade > CascadeDCAcut)continue;
                    AliESDcascade cascade(vertex,trackInKaon,Kaonidx);
                    Double_t Cascadeposition[3] = {0., 0., 0.};
                    cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                    Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                    if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                    Double_t pV0[3] = {0., 0., 0.};
                    vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                    fhCausalitycheck->Fill(0.5);
                    if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                    fhCausalitycheck->Fill(1.5);
                    if(r2cascade > r2)continue;
                    Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                    if(CascadeCPA < CascadeCPAcut)continue;
                    Double_t pCascade[3] = {0., 0., 0.};
                    cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                    Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                    if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                    Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                    Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                    Double_t pKaonAt[3] = {0., 0., 0.};
                    KaonTrack->GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                    Double_t pKaon[3] = {0., 0., 0.};
                    pKaon[0] = trackInKaon.Px();
                    pKaon[1] = trackInKaon.Py();
                    pKaon[2] = trackInKaon.Pz();
                    fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                    fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                    fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                    TLorentzVector vKaon, vsumCascade;
                    vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                    vsumCascade = vsumV0 + vKaon;
                    Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                    //----------leviate the cut----------
                    if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                    T1DecayChannel = (Int_t)3;
                    T1EventPercentile = (Float_t)percentile;
                    T1TVDCA = (Float_t)DCACascade;
                    T1TVCPA = (Float_t)CascadeCPA;
                    T1TVtoPV = (Float_t)DCACascadetoPV;
                    T1SVDCA = (Float_t)Dcadaughters;
                    T1SVCPA = (Float_t)CPAV0;
                    T1SVtoPV = (Float_t)SVtoPV;
                    T1SVtoTV = (Float_t)DCACascadetoSV;
                    T1KaonPx = (Float_t)pKaonAt[0];
                    T1KaonPy = (Float_t)pKaonAt[1];
                    T1KaonPz = (Float_t)pKaonAt[2];
                    T1KaonDCAtoPVXY = (Float_t)KaonImpactXY;
                    T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                    T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                    T1DeuteronPx = (Float_t)pDeuteronAt[0];
                    T1DeuteronPy = (Float_t)pDeuteronAt[1];
                    T1DeuteronPz = (Float_t)pDeuteronAt[2];
                    T1DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
                    T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                    T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                    T1PionPx = (Float_t)pPionAt[0];
                    T1PionPy = (Float_t)pPionAt[1];
                    T1PionPz = (Float_t)pPionAt[2];
                    T1PionDCAtoPVXY = (Float_t)Piond0Rot[0];
                    T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                    T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                    PiKDCascadeDecayTree->Fill();
                  }
                }
              }
            }
          }
        }
        else
        {
          for(Int_t j=0; j<Deuteronid.size(); j++)
          {
            Int_t Deuteronidx = Deuteronid[j];
            if(Pionidx == Deuteronidx)continue;
            AliESDtrack* DeuteronTrack = 0x0;
            DeuteronTrack = (fESDevent->GetTrack(Deuteronidx));
            if(DeuteronTrack->Charge() > 0)continue;
            AliExternalTrackParam trackInDeuteron(*DeuteronTrack);
            if(kRotateDeuteron)
            {
              Double_t pseudoXDeuteron[3], pseudoPDeuteron[3], CovPseudoDeuteron[21];
              Double_t pDeuteronRot[3];
              DeuteronTrack->GetCovarianceXYZPxPyPz(CovPseudoDeuteron);
              DeuteronTrack->GetXYZ(pseudoXDeuteron);
              DeuteronTrack->GetPxPyPz(pseudoPDeuteron);
              Short_t signDeuteron = DeuteronTrack->Charge();
              for(Int_t r_Deuteron=0; r_Deuteron<fRot; r_Deuteron++)
              {
                AliExternalTrackParam DeuteronTrackRotated;
                if(nRotation == 20)
                {
                  pDeuteronRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                  pDeuteronRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                  pDeuteronRot[2] = pseudoPDeuteron[2];
                  DeuteronTrackRotated = AliExternalTrackParam(pseudoXDeuteron,pDeuteronRot,CovPseudoDeuteron,signDeuteron);
                }
                else
                {
                  pDeuteronRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                  pDeuteronRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                  pDeuteronRot[2] = pseudoPDeuteron[2];
                  DeuteronTrackRotated = AliExternalTrackParam(pseudoXDeuteron,pDeuteronRot,CovPseudoDeuteron,signDeuteron);
                }
                Double_t Deuterond0Rot[2], Deuteroncovd0Rot[3];
                DeuteronTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Deuterond0Rot, Deuteroncovd0Rot);
                Double_t xn, xp;
                Double_t Dcadaughters = 0;
                Dcadaughters = DeuteronTrackRotated.GetDCA(&PionTrackRotated,lMagneticField,xn,xp);
                if(Dcadaughters > V0DCAdaughtercut)continue;
                PionTrackRotated.PropagateTo(xp,lMagneticField);
                DeuteronTrackRotated.PropagateTo(xn,lMagneticField);
                AliESDv0 vertex(DeuteronTrackRotated,Deuteronidx,PionTrackRotated,Pionidx);
                Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
                if(SVtoPV < V0DCAtoPVcut)continue;
                Double_t SVposition[3] = {0., 0., 0.};
                SVposition[0] = vertex.Xv();
                SVposition[1] = vertex.Yv();
                SVposition[2] = vertex.Zv();
                Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
                if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax)continue;
                Double_t pPionAt[3] = {0., 0., 0.};
                Double_t pDeuteronAt[3] = {0., 0., 0.};
                PionTrackRotated.GetPxPyPzAt(xp,lMagneticField,pPionAt);
                DeuteronTrackRotated.GetPxPyPzAt(xn,lMagneticField,pDeuteronAt);
                Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                if(CPAV0 < CascadeV0CPAcut)continue;
                TLorentzVector vDeuteron,vPion,vsumV0;
                vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
                vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
                vsumV0 = vDeuteron + vPion;
                Double_t V0Invmass = vsumV0.M();
                if(V0Invmass > CascadeV0Invmassmax)continue;
                if(DeuteronTrackRotated.Charge() < 0)
                {
                  for(Int_t k=0; k<Kaonid.size(); k++)
                  {
                    Int_t Kaonidx = Kaonid[k];
                    AliESDtrack* KaonTrack = 0x0;
                    KaonTrack = (fESDevent->GetTrack(Kaonidx));
                    AliExternalTrackParam trackInKaon(*KaonTrack);
                    if(kRotateKaon)
                    {
                      Double_t pseudoXKaon[3], pseudoPKaon[3], CovPseudoKaon[21];
                      Double_t pKaonRot[3];
                      KaonTrack->GetCovarianceXYZPxPyPz(CovPseudoKaon);
                      KaonTrack->GetXYZ(pseudoXKaon);
                      KaonTrack->GetPxPyPz(pseudoPKaon);
                      Short_t signKaon = KaonTrack->Charge();
                      for(Int_t r_Kaon=0; r_Kaon<fRot; r_Kaon++)
                      {
                        AliExternalTrackParam KaonTrackRotated;
                        if(nRotation == 20)
                        {
                          pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                          pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                          pKaonRot[2] = pseudoPKaon[2];
                          KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                        }
                        else
                        {
                          pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                          pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                          pKaonRot[2] = pseudoPKaon[2];
                          KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                        }
                        if(KaonTrackRotated.Charge() < 0)continue;
                        if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                        fhPiKDChargeCheckC2->Fill(0.5,PionTrackRotated.Charge());
                        fhPiKDChargeCheckC2->Fill(1.5,DeuteronTrackRotated.Charge());
                        fhPiKDChargeCheckC2->Fill(2.5,KaonTrackRotated.Charge());
                        Double_t Kaond0Rot[2], Kaoncovd0Rot[3];
                        KaonTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Kaond0Rot, Kaoncovd0Rot);
                        Double_t xkaon = 0;
                        Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&KaonTrackRotated,lMagneticField,&xkaon);
                        if(DCACascade > CascadeDCAcut)continue;
                        AliESDcascade cascade(vertex,KaonTrackRotated,Kaonidx);
                        Double_t Cascadeposition[3] = {0., 0., 0.};
                        cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                        Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                        if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                        Double_t pV0[3] = {0., 0., 0.};
                        vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                        fhCausalitycheck->Fill(0.5);
                        if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                        fhCausalitycheck->Fill(1.5);
                        if(r2cascade > r2)continue;
                        Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                        if(CascadeCPA < CascadeCPAcut)continue;
                        Double_t pCascade[3] = {0., 0., 0.};
                        cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                        Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                        if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                        Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                        Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                        Double_t pKaonAt[3] = {0., 0., 0.};
                        KaonTrackRotated.GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                        Double_t pKaon[3] = {0., 0., 0.};
                        pKaon[0] = KaonTrackRotated.Px();
                        pKaon[1] = KaonTrackRotated.Py();
                        pKaon[2] = KaonTrackRotated.Pz();
                        fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                        fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                        fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                        TLorentzVector vKaon, vsumCascade;
                        vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                        vsumCascade = vsumV0 + vKaon;
                        Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                        if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                        T1DecayChannel = (Int_t)4;
                        T1EventPercentile = (Float_t)percentile;
                        T1TVDCA = (Float_t)DCACascade;
                        T1TVCPA = (Float_t)CascadeCPA;
                        T1TVtoPV = (Float_t)DCACascadetoPV;
                        T1SVDCA = (Float_t)Dcadaughters;
                        T1SVCPA = (Float_t)CPAV0;
                        T1SVtoPV = (Float_t)SVtoPV;
                        T1SVtoTV = (Float_t)DCACascadetoSV;
                        T1KaonPx = (Float_t)pKaonAt[0];
                        T1KaonPy = (Float_t)pKaonAt[1];
                        T1KaonPz = (Float_t)pKaonAt[2];
                        T1KaonDCAtoPVXY = (Float_t)Kaond0Rot[0];
                        T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                        T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                        T1DeuteronPx = (Float_t)pDeuteronAt[0];
                        T1DeuteronPy = (Float_t)pDeuteronAt[1];
                        T1DeuteronPz = (Float_t)pDeuteronAt[2];
                        T1DeuteronDCAtoPVXY = (Float_t)Deuterond0Rot[0];
                        T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                        T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                        T1PionPx = (Float_t)pPionAt[0];
                        T1PionPy = (Float_t)pPionAt[1];
                        T1PionPz = (Float_t)pPionAt[2];
                        T1PionDCAtoPVXY = (Float_t)Piond0Rot[0];
                        T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                        T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                        PiKDCascadeDecayTree->Fill();
                      }
                    }
                    else
                    {
                      if(KaonTrack->Charge() < 0)continue;
                      if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                      fhPiKDChargeCheckC2->Fill(0.5,PionTrackRotated.Charge());
                      fhPiKDChargeCheckC2->Fill(1.5,DeuteronTrackRotated.Charge());
                      fhPiKDChargeCheckC2->Fill(2.5,KaonTrack->Charge());
                      Float_t KaonImpactXY;
                      Float_t KaonImpactZ;
                      KaonTrack->GetImpactParameters(KaonImpactXY,KaonImpactZ);
                      Double_t xkaon = 0;
                      Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&trackInKaon,lMagneticField,&xkaon);
                      if(DCACascade > CascadeDCAcut)continue;
                      AliESDcascade cascade(vertex,trackInKaon,Kaonidx);
                      Double_t Cascadeposition[3] = {0., 0., 0.};
                      cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                      Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                      if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                      Double_t pV0[3] = {0., 0., 0.};
                      vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                      fhCausalitycheck->Fill(0.5);
                      if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                      fhCausalitycheck->Fill(1.5);
                      if(r2cascade > r2)continue;
                      Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                      if(CascadeCPA < CascadeCPAcut)continue;
                      Double_t pCascade[3] = {0., 0., 0.};
                      cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                      Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                      if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                      Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                      Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                      Double_t pKaonAt[3] = {0., 0., 0.};
                      KaonTrack->GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                      Double_t pKaon[3] = {0., 0., 0.};
                      pKaon[0] = trackInKaon.Px();
                      pKaon[1] = trackInKaon.Py();
                      pKaon[2] = trackInKaon.Pz();
                      fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                      fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                      fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                      TLorentzVector vKaon, vsumCascade;
                      vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                      vsumCascade = vsumV0 + vKaon;
                      Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                      if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                      T1DecayChannel = (Int_t)4;
                      T1EventPercentile = (Float_t)percentile;
                      T1TVDCA = (Float_t)DCACascade;
                      T1TVCPA = (Float_t)CascadeCPA;
                      T1TVtoPV = (Float_t)DCACascadetoPV;
                      T1SVDCA = (Float_t)Dcadaughters;
                      T1SVCPA = (Float_t)CPAV0;
                      T1SVtoPV = (Float_t)SVtoPV;
                      T1SVtoTV = (Float_t)DCACascadetoSV;
                      T1KaonPx = (Float_t)pKaonAt[0];
                      T1KaonPy = (Float_t)pKaonAt[1];
                      T1KaonPz = (Float_t)pKaonAt[2];
                      T1KaonDCAtoPVXY = (Float_t)KaonImpactXY;
                      T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                      T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                      T1DeuteronPx = (Float_t)pDeuteronAt[0];
                      T1DeuteronPy = (Float_t)pDeuteronAt[1];
                      T1DeuteronPz = (Float_t)pDeuteronAt[2];
                      T1DeuteronDCAtoPVXY = (Float_t)Deuterond0Rot[0];
                      T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                      T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                      T1PionPx = (Float_t)pPionAt[0];
                      T1PionPy = (Float_t)pPionAt[1];
                      T1PionPz = (Float_t)pPionAt[2];
                      T1PionDCAtoPVXY = (Float_t)Piond0Rot[0];
                      T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                      T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                      PiKDCascadeDecayTree->Fill();
                    }
                  }
                }
              }
            }
            else
            {
              Float_t DeuteronImpactXY;
              Float_t DeuteronImpactZ;
              DeuteronTrack->GetImpactParameters(DeuteronImpactXY,DeuteronImpactZ);
              Double_t xn, xp;
              Double_t Dcadaughters = 0;
              Dcadaughters = trackInDeuteron.GetDCA(&PionTrackRotated,lMagneticField,xn,xp);
              if((xn+xp) < 2*V0Radiusmin || (xn+xp) > 2*V0Radiusmax)continue;
              if(Dcadaughters > V0DCAdaughtercut)continue;
              PionTrackRotated.PropagateTo(xp,lMagneticField);
              trackInDeuteron.PropagateTo(xn,lMagneticField);
              AliESDv0 vertex(trackInDeuteron,Deuteronidx,PionTrackRotated,Pionidx);
              Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
              if(SVtoPV < V0DCAtoPVcut)continue;
              Double_t SVposition[3] = {0., 0., 0.};
              SVposition[0] = vertex.Xv();
              SVposition[1] = vertex.Yv();
              SVposition[2] = vertex.Zv();
              Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
              if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax)continue;
              Double_t pPionAt[3] = {0., 0., 0.};
              Double_t pDeuteronAt[3] = {0., 0., 0.};
              PionTrackRotated.GetPxPyPzAt(xp,lMagneticField,pPionAt);
              DeuteronTrack->GetPxPyPzAt(xn,lMagneticField,pDeuteronAt);
              Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
              if(CPAV0 < CascadeV0CPAcut)continue;
              TLorentzVector vDeuteron,vPion,vsumV0;
              vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
              vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
              vsumV0 = vDeuteron + vPion;
              Double_t V0Invmass = vsumV0.M();
              if(V0Invmass > CascadeV0Invmassmax)continue;
              if(DeuteronTrack->Charge()<0)
              {
                for(Int_t k=0; k<Kaonid.size(); k++)
                {
                  Int_t Kaonidx = Kaonid[k];
                  AliESDtrack* KaonTrack = 0x0;
                  KaonTrack = (fESDevent->GetTrack(Kaonidx));
                  AliExternalTrackParam trackInKaon(*KaonTrack);
                  if(KaonTrack->Charge() < 0)continue;
                  if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                  fhPiKDChargeCheckC2->Fill(0.5,PionTrackRotated.Charge());
                  fhPiKDChargeCheckC2->Fill(1.5,DeuteronTrack->Charge());
                  fhPiKDChargeCheckC2->Fill(2.5,KaonTrack->Charge());
                  if(kRotateKaon)
                  {
                    Double_t pseudoXKaon[3], pseudoPKaon[3], CovPseudoKaon[21];
                    Double_t pKaonRot[3];
                    KaonTrack->GetCovarianceXYZPxPyPz(CovPseudoKaon);
                    KaonTrack->GetXYZ(pseudoXKaon);
                    KaonTrack->GetPxPyPz(pseudoPKaon);
                    Short_t signKaon = KaonTrack->Charge();
                    for(Int_t r_Kaon=0; r_Kaon<fRot; r_Kaon++)
                    {
                      AliExternalTrackParam KaonTrackRotated;
                      if(nRotation == 20)
                      {
                        pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[2] = pseudoPKaon[2];
                        KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                      }
                      else
                      {
                        pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[2] = pseudoPKaon[2];
                        KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                      }
                      if(KaonTrackRotated.Charge() < 0)continue;
                      if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                      fhPiKDChargeCheckC1->Fill(0.5,PionTrackRotated.Charge());
                      fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrack->Charge());
                      fhPiKDChargeCheckC1->Fill(2.5,KaonTrackRotated.Charge());
                      Double_t Kaond0Rot[2], Kaoncovd0Rot[3];
                      KaonTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Kaond0Rot, Kaoncovd0Rot);
                      Double_t xkaon = 0;
                      Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&KaonTrackRotated,lMagneticField,&xkaon);
                      if(DCACascade > CascadeDCAcut)continue;
                      AliESDcascade cascade(vertex,KaonTrackRotated,Kaonidx);
                      Double_t Cascadeposition[3] = {0., 0., 0.};
                      cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                      Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                      if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                      Double_t pV0[3] = {0., 0., 0.};
                      vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                      fhCausalitycheck->Fill(0.5);
                      if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                      fhCausalitycheck->Fill(1.5);
                      if(r2cascade > r2)continue;
                      Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                      if(CascadeCPA < CascadeCPAcut)continue;
                      Double_t pCascade[3] = {0., 0., 0.};
                      cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                      Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                      if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                      Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                      Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                      Double_t pKaonAt[3] = {0., 0., 0.};
                      KaonTrackRotated.GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                      Double_t pKaon[3] = {0., 0., 0.};
                      pKaon[0] = KaonTrackRotated.Px();
                      pKaon[1] = KaonTrackRotated.Py();
                      pKaon[2] = KaonTrackRotated.Pz();
                      fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                      fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                      fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                      TLorentzVector vKaon, vsumCascade;
                      vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                      vsumCascade = vsumV0 + vKaon;
                      Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                      if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                      T1DecayChannel = (Int_t)4;
                      T1EventPercentile = (Float_t)percentile;
                      T1TVDCA = (Float_t)DCACascade;
                      T1TVCPA = (Float_t)CascadeCPA;
                      T1TVtoPV = (Float_t)DCACascadetoPV;
                      T1SVDCA = (Float_t)Dcadaughters;
                      T1SVCPA = (Float_t)CPAV0;
                      T1SVtoPV = (Float_t)SVtoPV;
                      T1SVtoTV = (Float_t)DCACascadetoSV;
                      T1KaonPx = (Float_t)pKaonAt[0];
                      T1KaonPy = (Float_t)pKaonAt[1];
                      T1KaonPz = (Float_t)pKaonAt[2];
                      T1KaonDCAtoPVXY = (Float_t)Kaond0Rot[0];
                      T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                      T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                      T1DeuteronPx = (Float_t)pDeuteronAt[0];
                      T1DeuteronPy = (Float_t)pDeuteronAt[1];
                      T1DeuteronPz = (Float_t)pDeuteronAt[2];
                      T1DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
                      T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                      T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                      T1PionPx = (Float_t)pPionAt[0];
                      T1PionPy = (Float_t)pPionAt[1];
                      T1PionPz = (Float_t)pPionAt[2];
                      T1PionDCAtoPVXY = (Float_t)Piond0Rot[0];
                      T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                      T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                      PiKDCascadeDecayTree->Fill();
                    }
                  }
                  else
                  {
                    Float_t KaonImpactXY;
                    Float_t KaonImpactZ;
                    KaonTrack->GetImpactParameters(KaonImpactXY,KaonImpactZ);
                    Double_t xkaon = 0;
                    Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&trackInKaon,lMagneticField,&xkaon);
                    if(DCACascade > CascadeDCAcut)continue;
                    AliESDcascade cascade(vertex,trackInKaon,Kaonidx);
                    Double_t Cascadeposition[3] = {0., 0., 0.};
                    cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                    Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                    if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                    Double_t pV0[3] = {0., 0., 0.};
                    vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                    fhCausalitycheck->Fill(0.5);
                    if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                    fhCausalitycheck->Fill(1.5);
                    if(r2cascade > r2)continue;
                    Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                    if(CascadeCPA < CascadeCPAcut)continue;
                    Double_t pCascade[3] = {0., 0., 0.};
                    cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                    Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                    if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                    Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                    Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                    Double_t pKaonAt[3] = {0., 0., 0.};
                    KaonTrack->GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                    Double_t pKaon[3] = {0., 0., 0.};
                    pKaon[0] = trackInKaon.Px();
                    pKaon[1] = trackInKaon.Py();
                    pKaon[2] = trackInKaon.Pz();
                    fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                    fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                    fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                    TLorentzVector vKaon, vsumCascade;
                    vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                    vsumCascade = vsumV0 + vKaon;
                    Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                    //----------leviate the cut----------
                    if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                    T1DecayChannel = (Int_t)4;
                    T1EventPercentile = (Float_t)percentile;
                    T1TVDCA = (Float_t)DCACascade;
                    T1TVCPA = (Float_t)CascadeCPA;
                    T1TVtoPV = (Float_t)DCACascadetoPV;
                    T1SVDCA = (Float_t)Dcadaughters;
                    T1SVCPA = (Float_t)CPAV0;
                    T1SVtoPV = (Float_t)SVtoPV;
                    T1SVtoTV = (Float_t)DCACascadetoSV;
                    T1KaonPx = (Float_t)pKaonAt[0];
                    T1KaonPy = (Float_t)pKaonAt[1];
                    T1KaonPz = (Float_t)pKaonAt[2];
                    T1KaonDCAtoPVXY = (Float_t)KaonImpactXY;
                    T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                    T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                    T1DeuteronPx = (Float_t)pDeuteronAt[0];
                    T1DeuteronPy = (Float_t)pDeuteronAt[1];
                    T1DeuteronPz = (Float_t)pDeuteronAt[2];
                    T1DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
                    T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                    T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                    T1PionPx = (Float_t)pPionAt[0];
                    T1PionPy = (Float_t)pPionAt[1];
                    T1PionPz = (Float_t)pPionAt[2];
                    T1PionDCAtoPVXY = (Float_t)Piond0Rot[0];
                    T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                    T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                    PiKDCascadeDecayTree->Fill();
                  }
                }
              }
            }
          }
        }
      }
    }
    else
    {
      Float_t PionImpactXY =0;
      Float_t PionImpactZ =0;
      PionTrack->GetImpactParameters(PionImpactXY,PionImpactZ);
      Double_t fRot = nRotation;
      Double_t fAngle = 0.0872; // 6 gradi
      if(nRotation == 20)
      {
        fAngle = 0.1047;// 5 gradi
      }
      if(PionTrack->Charge() < 0)
      {
        for(Int_t j=0; j<Deuteronid.size(); j++)
        {
          Int_t Deuteronidx = Deuteronid[j];
          if(Pionidx == Deuteronidx)continue;
          AliESDtrack* DeuteronTrack = 0x0;
          DeuteronTrack = (fESDevent->GetTrack(Deuteronidx));
          if(DeuteronTrack->Charge() < 0)continue;
          AliExternalTrackParam trackInDeuteron(*DeuteronTrack);
          if(kRotateDeuteron)
          {
            Double_t pseudoXDeuteron[3], pseudoPDeuteron[3], CovPseudoDeuteron[21];
            Double_t pDeuteronRot[3];
            DeuteronTrack->GetCovarianceXYZPxPyPz(CovPseudoDeuteron);
            DeuteronTrack->GetXYZ(pseudoXDeuteron);
            DeuteronTrack->GetPxPyPz(pseudoPDeuteron);
            Short_t signDeuteron = DeuteronTrack->Charge();
            for(Int_t r_Deuteron=0; r_Deuteron<fRot; r_Deuteron++)
            {
              AliExternalTrackParam DeuteronTrackRotated;
              if(nRotation == 20)
              {
                pDeuteronRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                pDeuteronRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                pDeuteronRot[2] = pseudoPDeuteron[2];
                DeuteronTrackRotated = AliExternalTrackParam(pseudoXDeuteron,pDeuteronRot,CovPseudoDeuteron,signDeuteron);
              }
              else
              {
                pDeuteronRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                pDeuteronRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                pDeuteronRot[2] = pseudoPDeuteron[2];
                DeuteronTrackRotated = AliExternalTrackParam(pseudoXDeuteron,pDeuteronRot,CovPseudoDeuteron,signDeuteron);
              }
              Double_t Deuterond0Rot[2], Deuteroncovd0Rot[3];
              DeuteronTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Deuterond0Rot, Deuteroncovd0Rot);
              Double_t xn, xp;
              Double_t Dcadaughters = 0;
              Dcadaughters = PionTrack->GetDCA(&DeuteronTrackRotated,lMagneticField,xn,xp);
              if(Dcadaughters > V0DCAdaughtercut)continue;
              PionTrack->PropagateTo(xn,lMagneticField);
              DeuteronTrackRotated.PropagateTo(xp,lMagneticField);
              AliESDv0 vertex(trackInPion,Pionidx,DeuteronTrackRotated,Deuteronidx);
              Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
              if(SVtoPV < V0DCAtoPVcut)continue;
              Double_t SVposition[3] = {0., 0., 0.};
              SVposition[0] = vertex.Xv();
              SVposition[1] = vertex.Yv();
              SVposition[2] = vertex.Zv();
              Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
              if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax)continue;
              Double_t pPionAt[3] = {0., 0., 0.};
              Double_t pDeuteronAt[3] = {0., 0., 0.};
              PionTrack->GetPxPyPzAt(xn,lMagneticField,pPionAt);
              DeuteronTrackRotated.GetPxPyPzAt(xp,lMagneticField,pDeuteronAt);
              Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
              if(CPAV0 < CascadeV0CPAcut)continue;
              TLorentzVector vDeuteron,vPion,vsumV0;
              vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
              vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
              vsumV0 = vDeuteron + vPion;
              Double_t V0Invmass = vsumV0.M();
              if(V0Invmass > CascadeV0Invmassmax)continue;
              if(DeuteronTrackRotated.Charge() > 0)
              {
                for(Int_t k=0; k<Kaonid.size(); k++)
                {
                  Int_t Kaonidx = Kaonid[k];
                  AliESDtrack* KaonTrack = 0x0;
                  KaonTrack = (fESDevent->GetTrack(Kaonidx));
                  AliExternalTrackParam trackInKaon(*KaonTrack);
                  if(kRotateKaon)
                  {
                    Double_t pseudoXKaon[3], pseudoPKaon[3], CovPseudoKaon[21];
                    Double_t pKaonRot[3];
                    KaonTrack->GetCovarianceXYZPxPyPz(CovPseudoKaon);
                    KaonTrack->GetXYZ(pseudoXKaon);
                    KaonTrack->GetPxPyPz(pseudoPKaon);
                    Short_t signKaon = KaonTrack->Charge();
                    for(Int_t r_Kaon=0; r_Kaon<fRot; r_Kaon++)
                    {
                      AliExternalTrackParam KaonTrackRotated;
                      if(nRotation == 20)
                      {
                        pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[2] = pseudoPKaon[2];
                        KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                      }
                      else
                      {
                        pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[2] = pseudoPKaon[2];
                        KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                      }
                      if(KaonTrackRotated.Charge() > 0)continue;
                      if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                      fhPiKDChargeCheckC1->Fill(0.5,PionTrack->Charge());
                      fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrackRotated.Charge());
                      fhPiKDChargeCheckC1->Fill(2.5,KaonTrackRotated.Charge());
                      Double_t Kaond0Rot[2], Kaoncovd0Rot[3];
                      KaonTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Kaond0Rot, Kaoncovd0Rot);
                      Double_t xkaon = 0;
                      Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&KaonTrackRotated,lMagneticField,&xkaon);
                      if(DCACascade > CascadeDCAcut)continue;
                      AliESDcascade cascade(vertex,KaonTrackRotated,Kaonidx);
                      Double_t Cascadeposition[3] = {0., 0., 0.};
                      cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                      Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                      if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                      Double_t pV0[3] = {0., 0., 0.};
                      vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                      fhCausalitycheck->Fill(0.5);
                      if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                      fhCausalitycheck->Fill(1.5);
                      if(r2cascade > r2)continue;
                      Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                      if(CascadeCPA < CascadeCPAcut)continue;
                      Double_t pCascade[3] = {0., 0., 0.};
                      cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                      Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                      if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                      Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                      Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                      Double_t pKaonAt[3] = {0., 0., 0.};
                      KaonTrackRotated.GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                      Double_t pKaon[3] = {0., 0., 0.};
                      pKaon[0] = KaonTrackRotated.Px();
                      pKaon[1] = KaonTrackRotated.Py();
                      pKaon[2] = KaonTrackRotated.Pz();
                      fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                      fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                      fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                      TLorentzVector vKaon, vsumCascade;
                      vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                      vsumCascade = vsumV0 + vKaon;
                      Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                      if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                      T1DecayChannel = (Int_t)3;
                      T1EventPercentile = (Float_t)percentile;
                      T1TVDCA = (Float_t)DCACascade;
                      T1TVCPA = (Float_t)CascadeCPA;
                      T1TVtoPV = (Float_t)DCACascadetoPV;
                      T1SVDCA = (Float_t)Dcadaughters;
                      T1SVCPA = (Float_t)CPAV0;
                      T1SVtoPV = (Float_t)SVtoPV;
                      T1SVtoTV = (Float_t)DCACascadetoSV;
                      T1KaonPx = (Float_t)pKaonAt[0];
                      T1KaonPy = (Float_t)pKaonAt[1];
                      T1KaonPz = (Float_t)pKaonAt[2];
                      T1KaonDCAtoPVXY = (Float_t)Kaond0Rot[0];
                      T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                      T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                      T1DeuteronPx = (Float_t)pDeuteronAt[0];
                      T1DeuteronPy = (Float_t)pDeuteronAt[1];
                      T1DeuteronPz = (Float_t)pDeuteronAt[2];
                      T1DeuteronDCAtoPVXY = (Float_t)Deuterond0Rot[0];
                      T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                      T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                      T1PionPx = (Float_t)pPionAt[0];
                      T1PionPy = (Float_t)pPionAt[1];
                      T1PionPz = (Float_t)pPionAt[2];
                      T1PionDCAtoPVXY = (Float_t)PionImpactXY;
                      T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                      T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                      PiKDCascadeDecayTree->Fill();
                    }
                  }
                  else
                  {
                    if(KaonTrack->Charge() > 0)continue;
                    if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                    fhPiKDChargeCheckC1->Fill(0.5,PionTrack->Charge());
                    fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrackRotated.Charge());
                    fhPiKDChargeCheckC1->Fill(2.5,KaonTrack->Charge());
                    Float_t KaonImpactXY;
                    Float_t KaonImpactZ;
                    KaonTrack->GetImpactParameters(KaonImpactXY,KaonImpactZ);
                    Double_t xkaon = 0;
                    Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&trackInKaon,lMagneticField,&xkaon);
                    if(DCACascade > CascadeDCAcut)continue;
                    AliESDcascade cascade(vertex,trackInKaon,Kaonidx);
                    Double_t Cascadeposition[3] = {0., 0., 0.};
                    cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                    Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                    if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                    Double_t pV0[3] = {0., 0., 0.};
                    vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                    fhCausalitycheck->Fill(0.5);
                    if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                    fhCausalitycheck->Fill(1.5);
                    if(r2cascade > r2)continue;
                    Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                    if(CascadeCPA < CascadeCPAcut)continue;
                    Double_t pCascade[3] = {0., 0., 0.};
                    cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                    Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                    if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                    Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                    Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                    Double_t pKaonAt[3] = {0., 0., 0.};
                    KaonTrack->GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                    Double_t pKaon[3] = {0., 0., 0.};
                    pKaon[0] = trackInKaon.Px();
                    pKaon[1] = trackInKaon.Py();
                    pKaon[2] = trackInKaon.Pz();
                    fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                    fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                    fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                    TLorentzVector vKaon, vsumCascade;
                    vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                    vsumCascade = vsumV0 + vKaon;
                    Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                    if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                    T1DecayChannel = (Int_t)3;
                    T1EventPercentile = (Float_t)percentile;
                    T1TVDCA = (Float_t)DCACascade;
                    T1TVCPA = (Float_t)CascadeCPA;
                    T1TVtoPV = (Float_t)DCACascadetoPV;
                    T1SVDCA = (Float_t)Dcadaughters;
                    T1SVCPA = (Float_t)CPAV0;
                    T1SVtoPV = (Float_t)SVtoPV;
                    T1SVtoTV = (Float_t)DCACascadetoSV;
                    T1KaonPx = (Float_t)pKaonAt[0];
                    T1KaonPy = (Float_t)pKaonAt[1];
                    T1KaonPz = (Float_t)pKaonAt[2];
                    T1KaonDCAtoPVXY = (Float_t)KaonImpactXY;
                    T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                    T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                    T1DeuteronPx = (Float_t)pDeuteronAt[0];
                    T1DeuteronPy = (Float_t)pDeuteronAt[1];
                    T1DeuteronPz = (Float_t)pDeuteronAt[2];
                    T1DeuteronDCAtoPVXY = (Float_t)Deuterond0Rot[0];
                    T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                    T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                    T1PionPx = (Float_t)pPionAt[0];
                    T1PionPy = (Float_t)pPionAt[1];
                    T1PionPz = (Float_t)pPionAt[2];
                    T1PionDCAtoPVXY = (Float_t)PionImpactXY;
                    T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                    T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                    PiKDCascadeDecayTree->Fill();
                  }
                }
              }
            }
          }
          else
          {
            Float_t DeuteronImpactXY;
            Float_t DeuteronImpactZ;
            DeuteronTrack->GetImpactParameters(DeuteronImpactXY,DeuteronImpactZ);
            Double_t xn, xp;
            Double_t Dcadaughters = 0;
            Dcadaughters = PionTrack->GetDCA(&trackInDeuteron,lMagneticField,xn,xp);
            if((xn+xp) < 2*V0Radiusmin || (xn+xp) > 2*V0Radiusmax)continue;
            if(Dcadaughters > V0DCAdaughtercut)continue;
            PionTrack->PropagateTo(xn,lMagneticField);
            trackInDeuteron.PropagateTo(xp,lMagneticField);
            AliESDv0 vertex(trackInPion,Pionidx,trackInDeuteron,Deuteronidx);
            Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
            if(SVtoPV < V0DCAtoPVcut)continue;
            Double_t SVposition[3] = {0., 0., 0.};
            SVposition[0] = vertex.Xv();
            SVposition[1] = vertex.Yv();
            SVposition[2] = vertex.Zv();
            Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
            if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax)continue;
            Double_t pPionAt[3] = {0., 0., 0.};
            Double_t pDeuteronAt[3] = {0., 0., 0.};
            PionTrack->GetPxPyPzAt(xn,lMagneticField,pPionAt);
            DeuteronTrack->GetPxPyPzAt(xp,lMagneticField,pDeuteronAt);
            Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
            if(CPAV0 < CascadeV0CPAcut)continue;
            TLorentzVector vDeuteron,vPion,vsumV0;
            vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
            vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
            vsumV0 = vDeuteron + vPion;
            Double_t V0Invmass = vsumV0.M();
            if(V0Invmass > CascadeV0Invmassmax)continue;
            if(DeuteronTrack->Charge()>0)
            {
              for(Int_t k=0; k<Kaonid.size(); k++)
              {
                Int_t Kaonidx = Kaonid[k];
                AliESDtrack* KaonTrack = 0x0;
                KaonTrack = (fESDevent->GetTrack(Kaonidx));
                AliExternalTrackParam trackInKaon(*KaonTrack);
                if(KaonTrack->Charge() > 0)continue;
                if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                fhPiKDChargeCheckC1->Fill(0.5,PionTrack->Charge());
                fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrack->Charge());
                fhPiKDChargeCheckC1->Fill(2.5,KaonTrack->Charge());
                if(kRotateKaon)
                {
                  Double_t pseudoXKaon[3], pseudoPKaon[3], CovPseudoKaon[21];
                  Double_t pKaonRot[3];
                  KaonTrack->GetCovarianceXYZPxPyPz(CovPseudoKaon);
                  KaonTrack->GetXYZ(pseudoXKaon);
                  KaonTrack->GetPxPyPz(pseudoPKaon);
                  Short_t signKaon = KaonTrack->Charge();
                  for(Int_t r_Kaon=0; r_Kaon<fRot; r_Kaon++)
                  {
                    AliExternalTrackParam KaonTrackRotated;
                    if(nRotation == 20)
                    {
                      pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                      pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                      pKaonRot[2] = pseudoPKaon[2];
                      KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                    }
                    else
                    {
                      pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                      pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                      pKaonRot[2] = pseudoPKaon[2];
                      KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                    }
                    if(KaonTrackRotated.Charge() > 0)continue;
                    if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                    fhPiKDChargeCheckC1->Fill(0.5,PionTrack->Charge());
                    fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrack->Charge());
                    fhPiKDChargeCheckC1->Fill(2.5,KaonTrackRotated.Charge());
                    Double_t Kaond0Rot[2], Kaoncovd0Rot[3];
                    KaonTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Kaond0Rot, Kaoncovd0Rot);
                    Double_t xkaon = 0;
                    Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&KaonTrackRotated,lMagneticField,&xkaon);
                    if(DCACascade > CascadeDCAcut)continue;
                    AliESDcascade cascade(vertex,KaonTrackRotated,Kaonidx);
                    Double_t Cascadeposition[3] = {0., 0., 0.};
                    cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                    Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                    if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                    Double_t pV0[3] = {0., 0., 0.};
                    vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                    fhCausalitycheck->Fill(0.5);
                    if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                    fhCausalitycheck->Fill(1.5);
                    if(r2cascade > r2)continue;
                    Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                    if(CascadeCPA < CascadeCPAcut)continue;
                    Double_t pCascade[3] = {0., 0., 0.};
                    cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                    Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                    if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                    Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                    Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                    Double_t pKaonAt[3] = {0., 0., 0.};
                    KaonTrackRotated.GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                    Double_t pKaon[3] = {0., 0., 0.};
                    pKaon[0] = KaonTrackRotated.Px();
                    pKaon[1] = KaonTrackRotated.Py();
                    pKaon[2] = KaonTrackRotated.Pz();
                    fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                    fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                    fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                    TLorentzVector vKaon, vsumCascade;
                    vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                    vsumCascade = vsumV0 + vKaon;
                    Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                    if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                    T1DecayChannel = (Int_t)3;
                    T1EventPercentile = (Float_t)percentile;
                    T1TVDCA = (Float_t)DCACascade;
                    T1TVCPA = (Float_t)CascadeCPA;
                    T1TVtoPV = (Float_t)DCACascadetoPV;
                    T1SVDCA = (Float_t)Dcadaughters;
                    T1SVCPA = (Float_t)CPAV0;
                    T1SVtoPV = (Float_t)SVtoPV;
                    T1SVtoTV = (Float_t)DCACascadetoSV;
                    T1KaonPx = (Float_t)pKaonAt[0];
                    T1KaonPy = (Float_t)pKaonAt[1];
                    T1KaonPz = (Float_t)pKaonAt[2];
                    T1KaonDCAtoPVXY = (Float_t)Kaond0Rot[0];
                    T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                    T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                    T1DeuteronPx = (Float_t)pDeuteronAt[0];
                    T1DeuteronPy = (Float_t)pDeuteronAt[1];
                    T1DeuteronPz = (Float_t)pDeuteronAt[2];
                    T1DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
                    T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                    T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                    T1PionPx = (Float_t)pPionAt[0];
                    T1PionPy = (Float_t)pPionAt[1];
                    T1PionPz = (Float_t)pPionAt[2];
                    T1PionDCAtoPVXY = (Float_t)PionImpactXY;
                    T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                    T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                    PiKDCascadeDecayTree->Fill();
                  }
                }
                else
                {
                  Float_t KaonImpactXY;
                  Float_t KaonImpactZ;
                  KaonTrack->GetImpactParameters(KaonImpactXY,KaonImpactZ);
                  Double_t xkaon = 0;
                  Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&trackInKaon,lMagneticField,&xkaon);
                  if(DCACascade > CascadeDCAcut)continue;
                  AliESDcascade cascade(vertex,trackInKaon,Kaonidx);
                  Double_t Cascadeposition[3] = {0., 0., 0.};
                  cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                  Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                  if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                  Double_t pV0[3] = {0., 0., 0.};
                  vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                  fhCausalitycheck->Fill(0.5);
                  if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                  fhCausalitycheck->Fill(1.5);
                  if(r2cascade > r2)continue;
                  Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                  if(CascadeCPA < CascadeCPAcut)continue;
                  Double_t pCascade[3] = {0., 0., 0.};
                  cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                  Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                  if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                  Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                  Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                  Double_t pKaonAt[3] = {0., 0., 0.};
                  KaonTrack->GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                  Double_t pKaon[3] = {0., 0., 0.};
                  pKaon[0] = trackInKaon.Px();
                  pKaon[1] = trackInKaon.Py();
                  pKaon[2] = trackInKaon.Pz();
                  fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                  fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                  fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                  TLorentzVector vKaon, vsumCascade;
                  vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                  vsumCascade = vsumV0 + vKaon;
                  Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                  //----------leviate the cut----------
                  if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                  T1DecayChannel = (Int_t)3;
                  T1EventPercentile = (Float_t)percentile;
                  T1TVDCA = (Float_t)DCACascade;
                  T1TVCPA = (Float_t)CascadeCPA;
                  T1TVtoPV = (Float_t)DCACascadetoPV;
                  T1SVDCA = (Float_t)Dcadaughters;
                  T1SVCPA = (Float_t)CPAV0;
                  T1SVtoPV = (Float_t)SVtoPV;
                  T1SVtoTV = (Float_t)DCACascadetoSV;
                  T1KaonPx = (Float_t)pKaonAt[0];
                  T1KaonPy = (Float_t)pKaonAt[1];
                  T1KaonPz = (Float_t)pKaonAt[2];
                  T1KaonDCAtoPVXY = (Float_t)KaonImpactXY;
                  T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                  T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                  T1DeuteronPx = (Float_t)pDeuteronAt[0];
                  T1DeuteronPy = (Float_t)pDeuteronAt[1];
                  T1DeuteronPz = (Float_t)pDeuteronAt[2];
                  T1DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
                  T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                  T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                  T1PionPx = (Float_t)pPionAt[0];
                  T1PionPy = (Float_t)pPionAt[1];
                  T1PionPz = (Float_t)pPionAt[2];
                  T1PionDCAtoPVXY = (Float_t)PionImpactXY;
                  T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                  T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                  PiKDCascadeDecayTree->Fill();
                }
              }
            }
          }
        }
      }
      else
      {
        for(Int_t j=0; j<Deuteronid.size(); j++)
        {
          Int_t Deuteronidx = Deuteronid[j];
          if(Pionidx == Deuteronidx)continue;
          AliESDtrack* DeuteronTrack = 0x0;
          DeuteronTrack = (fESDevent->GetTrack(Deuteronidx));
          if(DeuteronTrack->Charge() > 0)continue;
          AliExternalTrackParam trackInDeuteron(*DeuteronTrack);
          if(kRotateDeuteron)
          {
            Double_t pseudoXDeuteron[3], pseudoPDeuteron[3], CovPseudoDeuteron[21];
            Double_t pDeuteronRot[3];
            DeuteronTrack->GetCovarianceXYZPxPyPz(CovPseudoDeuteron);
            DeuteronTrack->GetXYZ(pseudoXDeuteron);
            DeuteronTrack->GetPxPyPz(pseudoPDeuteron);
            Short_t signDeuteron = DeuteronTrack->Charge();
            for(Int_t r_Deuteron=0; r_Deuteron<fRot; r_Deuteron++)
            {
              AliExternalTrackParam DeuteronTrackRotated;
              if(nRotation == 20)
              {
                pDeuteronRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                pDeuteronRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                pDeuteronRot[2] = pseudoPDeuteron[2];
                DeuteronTrackRotated = AliExternalTrackParam(pseudoXDeuteron,pDeuteronRot,CovPseudoDeuteron,signDeuteron);
              }
              else
              {
                pDeuteronRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                pDeuteronRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Deuteron*fAngle)*pseudoPDeuteron[1];
                pDeuteronRot[2] = pseudoPDeuteron[2];
                DeuteronTrackRotated = AliExternalTrackParam(pseudoXDeuteron,pDeuteronRot,CovPseudoDeuteron,signDeuteron);
              }
              Double_t Deuterond0Rot[2], Deuteroncovd0Rot[3];
              DeuteronTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Deuterond0Rot, Deuteroncovd0Rot);
              Double_t xn, xp;
              Double_t Dcadaughters = 0;
              Dcadaughters = DeuteronTrackRotated.GetDCA(PionTrack,lMagneticField,xn,xp);
              if(Dcadaughters > V0DCAdaughtercut)continue;
              PionTrack->PropagateTo(xp,lMagneticField);
              DeuteronTrackRotated.PropagateTo(xn,lMagneticField);
              AliESDv0 vertex(DeuteronTrackRotated,Deuteronidx,trackInPion,Pionidx);
              Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
              if(SVtoPV < V0DCAtoPVcut)continue;
              Double_t SVposition[3] = {0., 0., 0.};
              SVposition[0] = vertex.Xv();
              SVposition[1] = vertex.Yv();
              SVposition[2] = vertex.Zv();
              Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
              if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax)continue;
              Double_t pPionAt[3] = {0., 0., 0.};
              Double_t pDeuteronAt[3] = {0., 0., 0.};
              PionTrack->GetPxPyPzAt(xp,lMagneticField,pPionAt);
              DeuteronTrackRotated.GetPxPyPzAt(xn,lMagneticField,pDeuteronAt);
              Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
              if(CPAV0 < CascadeV0CPAcut)continue;
              TLorentzVector vDeuteron,vPion,vsumV0;
              vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
              vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
              vsumV0 = vDeuteron + vPion;
              Double_t V0Invmass = vsumV0.M();
              if(V0Invmass > CascadeV0Invmassmax)continue;
              if(DeuteronTrackRotated.Charge() < 0)
              {
                for(Int_t k=0; k<Kaonid.size(); k++)
                {
                  Int_t Kaonidx = Kaonid[k];
                  AliESDtrack* KaonTrack = 0x0;
                  KaonTrack = (fESDevent->GetTrack(Kaonidx));
                  AliExternalTrackParam trackInKaon(*KaonTrack);
                  if(kRotateKaon)
                  {
                    Double_t pseudoXKaon[3], pseudoPKaon[3], CovPseudoKaon[21];
                    Double_t pKaonRot[3];
                    KaonTrack->GetCovarianceXYZPxPyPz(CovPseudoKaon);
                    KaonTrack->GetXYZ(pseudoXKaon);
                    KaonTrack->GetPxPyPz(pseudoPKaon);
                    Short_t signKaon = KaonTrack->Charge();
                    for(Int_t r_Kaon=0; r_Kaon<fRot; r_Kaon++)
                    {
                      AliExternalTrackParam KaonTrackRotated;
                      if(nRotation == 20)
                      {
                        pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[2] = pseudoPKaon[2];
                        KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                      }
                      else
                      {
                        pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                        pKaonRot[2] = pseudoPKaon[2];
                        KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                      }
                      if(KaonTrackRotated.Charge() < 0)continue;
                      if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                      fhPiKDChargeCheckC2->Fill(0.5,PionTrack->Charge());
                      fhPiKDChargeCheckC2->Fill(1.5,DeuteronTrackRotated.Charge());
                      fhPiKDChargeCheckC2->Fill(2.5,KaonTrackRotated.Charge());
                      Double_t Kaond0Rot[2], Kaoncovd0Rot[3];
                      KaonTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Kaond0Rot, Kaoncovd0Rot);
                      Double_t xkaon = 0;
                      Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&KaonTrackRotated,lMagneticField,&xkaon);
                      if(DCACascade > CascadeDCAcut)continue;
                      AliESDcascade cascade(vertex,KaonTrackRotated,Kaonidx);
                      Double_t Cascadeposition[3] = {0., 0., 0.};
                      cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                      Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                      if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                      Double_t pV0[3] = {0., 0., 0.};
                      vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                      fhCausalitycheck->Fill(0.5);
                      if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                      fhCausalitycheck->Fill(1.5);
                      if(r2cascade > r2)continue;
                      Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                      if(CascadeCPA < CascadeCPAcut)continue;
                      Double_t pCascade[3] = {0., 0., 0.};
                      cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                      Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                      if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                      Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                      Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                      Double_t pKaonAt[3] = {0., 0., 0.};
                      KaonTrackRotated.GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                      Double_t pKaon[3] = {0., 0., 0.};
                      pKaon[0] = KaonTrackRotated.Px();
                      pKaon[1] = KaonTrackRotated.Py();
                      pKaon[2] = KaonTrackRotated.Pz();
                      fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                      fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                      fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                      TLorentzVector vKaon, vsumCascade;
                      vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                      vsumCascade = vsumV0 + vKaon;
                      Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                      if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                      T1DecayChannel = (Int_t)4;
                      T1EventPercentile = (Float_t)percentile;
                      T1TVDCA = (Float_t)DCACascade;
                      T1TVCPA = (Float_t)CascadeCPA;
                      T1TVtoPV = (Float_t)DCACascadetoPV;
                      T1SVDCA = (Float_t)Dcadaughters;
                      T1SVCPA = (Float_t)CPAV0;
                      T1SVtoPV = (Float_t)SVtoPV;
                      T1SVtoTV = (Float_t)DCACascadetoSV;
                      T1KaonPx = (Float_t)pKaonAt[0];
                      T1KaonPy = (Float_t)pKaonAt[1];
                      T1KaonPz = (Float_t)pKaonAt[2];
                      T1KaonDCAtoPVXY = (Float_t)Kaond0Rot[0];
                      T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                      T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                      T1DeuteronPx = (Float_t)pDeuteronAt[0];
                      T1DeuteronPy = (Float_t)pDeuteronAt[1];
                      T1DeuteronPz = (Float_t)pDeuteronAt[2];
                      T1DeuteronDCAtoPVXY = (Float_t)Deuterond0Rot[0];
                      T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                      T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                      T1PionPx = (Float_t)pPionAt[0];
                      T1PionPy = (Float_t)pPionAt[1];
                      T1PionPz = (Float_t)pPionAt[2];
                      T1PionDCAtoPVXY = (Float_t)PionImpactXY;
                      T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                      T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                      PiKDCascadeDecayTree->Fill();
                    }
                  }
                  else
                  {
                    if(KaonTrack->Charge() < 0)continue;
                    if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                    fhPiKDChargeCheckC2->Fill(0.5,PionTrack->Charge());
                    fhPiKDChargeCheckC2->Fill(1.5,DeuteronTrackRotated.Charge());
                    fhPiKDChargeCheckC2->Fill(2.5,KaonTrack->Charge());
                    Float_t KaonImpactXY;
                    Float_t KaonImpactZ;
                    KaonTrack->GetImpactParameters(KaonImpactXY,KaonImpactZ);
                    Double_t xkaon = 0;
                    Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&trackInKaon,lMagneticField,&xkaon);
                    if(DCACascade > CascadeDCAcut)continue;
                    AliESDcascade cascade(vertex,trackInKaon,Kaonidx);
                    Double_t Cascadeposition[3] = {0., 0., 0.};
                    cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                    Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                    if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                    Double_t pV0[3] = {0., 0., 0.};
                    vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                    fhCausalitycheck->Fill(0.5);
                    if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                    fhCausalitycheck->Fill(1.5);
                    if(r2cascade > r2)continue;
                    Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                    if(CascadeCPA < CascadeCPAcut)continue;
                    Double_t pCascade[3] = {0., 0., 0.};
                    cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                    Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                    if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                    Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                    Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                    Double_t pKaonAt[3] = {0., 0., 0.};
                    KaonTrack->GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                    Double_t pKaon[3] = {0., 0., 0.};
                    pKaon[0] = trackInKaon.Px();
                    pKaon[1] = trackInKaon.Py();
                    pKaon[2] = trackInKaon.Pz();
                    fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                    fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                    fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                    TLorentzVector vKaon, vsumCascade;
                    vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                    vsumCascade = vsumV0 + vKaon;
                    Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                    if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                    T1DecayChannel = (Int_t)4;
                    T1EventPercentile = (Float_t)percentile;
                    T1TVDCA = (Float_t)DCACascade;
                    T1TVCPA = (Float_t)CascadeCPA;
                    T1TVtoPV = (Float_t)DCACascadetoPV;
                    T1SVDCA = (Float_t)Dcadaughters;
                    T1SVCPA = (Float_t)CPAV0;
                    T1SVtoPV = (Float_t)SVtoPV;
                    T1SVtoTV = (Float_t)DCACascadetoSV;
                    T1KaonPx = (Float_t)pKaonAt[0];
                    T1KaonPy = (Float_t)pKaonAt[1];
                    T1KaonPz = (Float_t)pKaonAt[2];
                    T1KaonDCAtoPVXY = (Float_t)KaonImpactXY;
                    T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                    T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                    T1DeuteronPx = (Float_t)pDeuteronAt[0];
                    T1DeuteronPy = (Float_t)pDeuteronAt[1];
                    T1DeuteronPz = (Float_t)pDeuteronAt[2];
                    T1DeuteronDCAtoPVXY = (Float_t)Deuterond0Rot[0];
                    T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                    T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                    T1PionPx = (Float_t)pPionAt[0];
                    T1PionPy = (Float_t)pPionAt[1];
                    T1PionPz = (Float_t)pPionAt[2];
                    T1PionDCAtoPVXY = (Float_t)PionImpactXY;
                    T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                    T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                    PiKDCascadeDecayTree->Fill();
                  }
                }
              }
            }
          }
          else
          {
            Float_t DeuteronImpactXY;
            Float_t DeuteronImpactZ;
            DeuteronTrack->GetImpactParameters(DeuteronImpactXY,DeuteronImpactZ);
            Double_t xn, xp;
            Double_t Dcadaughters = 0;
            Dcadaughters = trackInDeuteron.GetDCA(PionTrack,lMagneticField,xn,xp);
            if((xn+xp) < 2*V0Radiusmin || (xn+xp) > 2*V0Radiusmax)continue;
            if(Dcadaughters > V0DCAdaughtercut)continue;
            PionTrack->PropagateTo(xp,lMagneticField);
            trackInDeuteron.PropagateTo(xn,lMagneticField);
            AliESDv0 vertex(trackInDeuteron,Deuteronidx,trackInPion,Pionidx);
            Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
            if(SVtoPV < V0DCAtoPVcut)continue;
            Double_t SVposition[3] = {0., 0., 0.};
            SVposition[0] = vertex.Xv();
            SVposition[1] = vertex.Yv();
            SVposition[2] = vertex.Zv();
            Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
            if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax)continue;
            Double_t pPionAt[3] = {0., 0., 0.};
            Double_t pDeuteronAt[3] = {0., 0., 0.};
            PionTrack->GetPxPyPzAt(xp,lMagneticField,pPionAt);
            DeuteronTrack->GetPxPyPzAt(xn,lMagneticField,pDeuteronAt);
            Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
            if(CPAV0 < CascadeV0CPAcut)continue;
            TLorentzVector vDeuteron,vPion,vsumV0;
            vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
            vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
            vsumV0 = vDeuteron + vPion;
            Double_t V0Invmass = vsumV0.M();
            if(V0Invmass > CascadeV0Invmassmax)continue;
            if(DeuteronTrack->Charge()<0)
            {
              for(Int_t k=0; k<Kaonid.size(); k++)
              {
                Int_t Kaonidx = Kaonid[k];
                AliESDtrack* KaonTrack = 0x0;
                KaonTrack = (fESDevent->GetTrack(Kaonidx));
                AliExternalTrackParam trackInKaon(*KaonTrack);
                if(KaonTrack->Charge() < 0)continue;
                if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                fhPiKDChargeCheckC2->Fill(0.5,PionTrack->Charge());
                fhPiKDChargeCheckC2->Fill(1.5,DeuteronTrack->Charge());
                fhPiKDChargeCheckC2->Fill(2.5,KaonTrack->Charge());
                if(kRotateKaon)
                {
                  Double_t pseudoXKaon[3], pseudoPKaon[3], CovPseudoKaon[21];
                  Double_t pKaonRot[3];
                  KaonTrack->GetCovarianceXYZPxPyPz(CovPseudoKaon);
                  KaonTrack->GetXYZ(pseudoXKaon);
                  KaonTrack->GetPxPyPz(pseudoPKaon);
                  Short_t signKaon = KaonTrack->Charge();
                  for(Int_t r_Kaon=0; r_Kaon<fRot; r_Kaon++)
                  {
                    AliExternalTrackParam KaonTrackRotated;
                    if(nRotation == 20)
                    {
                      pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                      pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                      pKaonRot[2] = pseudoPKaon[2];
                      KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                    }
                    else
                    {
                      pKaonRot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                      pKaonRot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r_Kaon*fAngle)*pseudoPKaon[1];
                      pKaonRot[2] = pseudoPKaon[2];
                      KaonTrackRotated = AliExternalTrackParam(pseudoXKaon,pKaonRot,CovPseudoKaon,signKaon);
                    }
                    if(KaonTrackRotated.Charge() < 0)continue;
                    if(Kaonidx == Pionidx || Kaonidx == Deuteronidx)continue;
                    fhPiKDChargeCheckC1->Fill(0.5,PionTrack->Charge());
                    fhPiKDChargeCheckC1->Fill(1.5,DeuteronTrack->Charge());
                    fhPiKDChargeCheckC1->Fill(2.5,KaonTrackRotated.Charge());
                    Double_t Kaond0Rot[2], Kaoncovd0Rot[3];
                    KaonTrackRotated.PropagateToDCA(primaryVertex, lMagneticField, 100., Kaond0Rot, Kaoncovd0Rot);
                    Double_t xkaon = 0;
                    Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&KaonTrackRotated,lMagneticField,&xkaon);
                    if(DCACascade > CascadeDCAcut)continue;
                    AliESDcascade cascade(vertex,KaonTrackRotated,Kaonidx);
                    Double_t Cascadeposition[3] = {0., 0., 0.};
                    cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                    Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                    if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                    Double_t pV0[3] = {0., 0., 0.};
                    vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                    fhCausalitycheck->Fill(0.5);
                    if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                    fhCausalitycheck->Fill(1.5);
                    if(r2cascade > r2)continue;
                    Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                    if(CascadeCPA < CascadeCPAcut)continue;
                    Double_t pCascade[3] = {0., 0., 0.};
                    cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                    Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                    if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                    Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                    Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                    Double_t pKaonAt[3] = {0., 0., 0.};
                    KaonTrackRotated.GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                    Double_t pKaon[3] = {0., 0., 0.};
                    pKaon[0] = KaonTrackRotated.Px();
                    pKaon[1] = KaonTrackRotated.Py();
                    pKaon[2] = KaonTrackRotated.Pz();
                    fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                    fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                    fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                    TLorentzVector vKaon, vsumCascade;
                    vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                    vsumCascade = vsumV0 + vKaon;
                    Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                    if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                    T1DecayChannel = (Int_t)4;
                    T1EventPercentile = (Float_t)percentile;
                    T1TVDCA = (Float_t)DCACascade;
                    T1TVCPA = (Float_t)CascadeCPA;
                    T1TVtoPV = (Float_t)DCACascadetoPV;
                    T1SVDCA = (Float_t)Dcadaughters;
                    T1SVCPA = (Float_t)CPAV0;
                    T1SVtoPV = (Float_t)SVtoPV;
                    T1SVtoTV = (Float_t)DCACascadetoSV;
                    T1KaonPx = (Float_t)pKaonAt[0];
                    T1KaonPy = (Float_t)pKaonAt[1];
                    T1KaonPz = (Float_t)pKaonAt[2];
                    T1KaonDCAtoPVXY = (Float_t)Kaond0Rot[0];
                    T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                    T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                    T1DeuteronPx = (Float_t)pDeuteronAt[0];
                    T1DeuteronPy = (Float_t)pDeuteronAt[1];
                    T1DeuteronPz = (Float_t)pDeuteronAt[2];
                    T1DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
                    T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                    T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                    T1PionPx = (Float_t)pPionAt[0];
                    T1PionPy = (Float_t)pPionAt[1];
                    T1PionPz = (Float_t)pPionAt[2];
                    T1PionDCAtoPVXY = (Float_t)PionImpactXY;
                    T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                    T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                    PiKDCascadeDecayTree->Fill();
                  }
                }
                else
                {
                  Float_t KaonImpactXY;
                  Float_t KaonImpactZ;
                  KaonTrack->GetImpactParameters(KaonImpactXY,KaonImpactZ);
                  Double_t xkaon = 0;
                  Double_t DCACascade = AliAnalysisTaskHyperFinder3Body::PropagateToDCA(&vertex,&trackInKaon,lMagneticField,&xkaon);
                  if(DCACascade > CascadeDCAcut)continue;
                  AliESDcascade cascade(vertex,trackInKaon,Kaonidx);
                  Double_t Cascadeposition[3] = {0., 0., 0.};
                  cascade.GetXYZcascade(Cascadeposition[0],Cascadeposition[1],Cascadeposition[2]);
                  Double_t r2cascade = Cascadeposition[0]*Cascadeposition[0]+Cascadeposition[1]*Cascadeposition[1];
                  if(r2cascade < CascadeRadiusmin*CascadeRadiusmin || r2cascade > CascadeRadiusmax*CascadeRadiusmax)continue;
                  Double_t pV0[3] = {0., 0., 0.};
                  vertex.GetPxPyPz(pV0[0],pV0[1],pV0[2]);
                  fhCausalitycheck->Fill(0.5);
                  if(Cascadeposition[0]*pV0[0]+Cascadeposition[1]*pV0[1]+Cascadeposition[2]*pV0[2] < 0)continue;
                  fhCausalitycheck->Fill(1.5);
                  if(r2cascade > r2)continue;
                  Double_t CascadeCPA = cascade.GetCascadeCosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
                  if(CascadeCPA < CascadeCPAcut)continue;
                  Double_t pCascade[3] = {0., 0., 0.};
                  cascade.GetPxPyPz(pCascade[0],pCascade[1],pCascade[2]);
                  Double_t ptCascade = TMath::Sqrt(pCascade[0]*pCascade[0]+pCascade[1]*pCascade[1]);
                  if(ptCascade < Cascadeptmin || ptCascade > Cascadeptmax)continue;
                  Double_t DCACascadetoPV = cascade.GetDcascade(PVposition[0],PVposition[1],PVposition[2]);
                  Double_t DCACascadetoSV = cascade.GetDcascade(SVposition[0],SVposition[1],SVposition[2]);
                  Double_t pKaonAt[3] = {0., 0., 0.};
                  KaonTrack->GetPxPyPzAt(xkaon,lMagneticField,pKaonAt);
                  Double_t pKaon[3] = {0., 0., 0.};
                  pKaon[0] = trackInKaon.Px();
                  pKaon[1] = trackInKaon.Py();
                  pKaon[2] = trackInKaon.Pz();
                  fhCascadePropagationPcheck->Fill(0.5,pKaon[0]-pKaonAt[0]);
                  fhCascadePropagationPcheck->Fill(1.5,pKaon[1]-pKaonAt[1]);
                  fhCascadePropagationPcheck->Fill(2.5,pKaon[2]-pKaonAt[2]);
                  TLorentzVector vKaon, vsumCascade;
                  vKaon.SetXYZM(pKaonAt[0],pKaonAt[1],pKaonAt[2],KaonMass);
                  vsumCascade = vsumV0 + vKaon;
                  Float_t InvmassCascade =  (Float_t)vsumCascade.M();
                  //----------leviate the cut----------
                  if(InvmassCascade < CascadeInvmassmin || InvmassCascade > CascadeInvmassmax)continue;
                  T1DecayChannel = (Int_t)4;
                  T1EventPercentile = (Float_t)percentile;
                  T1TVDCA = (Float_t)DCACascade;
                  T1TVCPA = (Float_t)CascadeCPA;
                  T1TVtoPV = (Float_t)DCACascadetoPV;
                  T1SVDCA = (Float_t)Dcadaughters;
                  T1SVCPA = (Float_t)CPAV0;
                  T1SVtoPV = (Float_t)SVtoPV;
                  T1SVtoTV = (Float_t)DCACascadetoSV;
                  T1KaonPx = (Float_t)pKaonAt[0];
                  T1KaonPy = (Float_t)pKaonAt[1];
                  T1KaonPz = (Float_t)pKaonAt[2];
                  T1KaonDCAtoPVXY = (Float_t)KaonImpactXY;
                  T1KaonTPCNcls = (Float_t)KaonTrack->GetTPCNcls();
                  T1KaonTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(KaonTrack,AliPID::kKaon);
                  T1DeuteronPx = (Float_t)pDeuteronAt[0];
                  T1DeuteronPy = (Float_t)pDeuteronAt[1];
                  T1DeuteronPz = (Float_t)pDeuteronAt[2];
                  T1DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
                  T1DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
                  T1DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
                  T1PionPx = (Float_t)pPionAt[0];
                  T1PionPy = (Float_t)pPionAt[1];
                  T1PionPz = (Float_t)pPionAt[2];
                  T1PionDCAtoPVXY = (Float_t)PionImpactXY;
                  T1PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
                  T1PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
                  PiKDCascadeDecayTree->Fill();
                }
              }
            }
          }
        }
      }
    }
  }

  /*for(Int_t i=0; i<Pionid.size(); i++)
  {
    AliESDtrack *PionTrack = 0x0;
    Int_t Pionidx = Pionid[i];
    PionTrack = fESDevent->GetTrack(Pionidx);
    if(kOpenTrackImpacCutV0 && PionTrack->GetD(PVposition[0],PVposition[1],lMagneticField) < TrackDCAtoPVcut) continue; 
    AliExternalTrackParam trackInPion(*PionTrack);
    Float_t PionImpactXY =0;
    Float_t PionImpactZ =0;
    PionTrack->GetImpactParameters(PionImpactXY,PionImpactZ);
    if(PionTrack->Charge() < 0)
    {
      for(Int_t j=0; j<Deuteronid.size(); j++)
      {
        AliESDtrack *DeuteronTrack = 0x0;
        Int_t Deuteronidx = Deuteronid[j];
        if(Pionidx == Deuteronidx)continue;
        DeuteronTrack = fESDevent->GetTrack(Deuteronidx);
        if(kOpenTrackImpacCutV0 && DeuteronTrack->GetD(PVposition[0],PVposition[1],lMagneticField) < TrackDCAtoPVcut) continue; 
        AliExternalTrackParam trackInDeuteron(*DeuteronTrack);
        Float_t DeuteronImpactXY;
        Float_t DeuteronImpactZ;
        DeuteronTrack->GetImpactParameters(DeuteronImpactXY,DeuteronImpactZ);
        Double_t xn, xp;
        Double_t Dcadaughters = 0;
        Dcadaughters = PionTrack->GetDCA(DeuteronTrack,lMagneticField,xn,xp);
        if((xn+xp) < 2*V0Radiusmin || (xn+xp) > 2*V0Radiusmax) continue;
        Bool_t corrected=kFALSE;
        if((trackInPion.GetX() > 3.) && (xn < 3.)) {corrected=kTRUE;}
        if((trackInDeuteron.GetX() > 3.) && (xp < 3.)) {corrected=kTRUE;}
        if(corrected && PionTrack->Charge()<0) {Dcadaughters = trackInPion.GetDCA(&trackInDeuteron,lMagneticField,xn,xp);}
        if(Dcadaughters > V0DCAdaughtercut) continue;
        trackInPion.PropagateTo(xn,lMagneticField);
        trackInDeuteron.PropagateTo(xp,lMagneticField);
        AliESDv0 vertex(trackInPion,Pionidx,trackInDeuteron,Deuteronidx);
        Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
        //if(SVtoPV < V0DCAtoPVcut)continue;
        Double_t SVposition[3] = {0., 0., 0.};
        SVposition[0] = vertex.Xv();
        SVposition[1] = vertex.Yv();
        SVposition[2] = vertex.Zv();
        Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
        if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax) continue;
        Double_t pPionAt[3] = {0., 0., 0.};
        Double_t pDeuteronAt[3] = {0., 0., 0.};
        PionTrack->GetPxPyPzAt(xn,lMagneticField,pPionAt);
        DeuteronTrack->GetPxPyPzAt(xp,lMagneticField,pDeuteronAt); 
        Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
        if(CPAV0 < V0CPAcut)continue;
        TLorentzVector vDeuteron,vPion,vsumV0;
        vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
        vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
        vsumV0 = vDeuteron + vPion;
        Float_t InvmassV0 = (Float_t)vsumV0.M();
        if(DeuteronTrack->Charge()>0 && kEnablePiDChannel1)
        {
          fhPiDChargeCheckC1->Fill(0.5,PionTrack->Charge());
          fhPiDChargeCheckC1->Fill(1.5,DeuteronTrack->Charge());
          if(InvmassV0 < V0Invmassmin || InvmassV0 > V0Invmassmax)continue;
          T2DecayChannel = (Int_t)1;
          T2EventPercentile = (Float_t)percentile;
          T2SVDCA = (Float_t)Dcadaughters;
          T2SVCPA = (Float_t)CPAV0;
          T2SVtoPV = (Float_t)SVtoPV;
          T2DeuteronPx = (Float_t)pDeuteronAt[0];
          T2DeuteronPy = (Float_t)pDeuteronAt[1];
          T2DeuteronPz = (Float_t)pDeuteronAt[2];
          T2DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
          T2DeuteronDCAtoPVZ = (Float_t)DeuteronImpactZ;
          T2DeuteronDCAtoSVXY = (Float_t)DeuteronTrack->GetD(SVposition[0],SVposition[1],lMagneticField);
          T2DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
          T2DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
          T2PionPx = (Float_t)pPionAt[0];
          T2PionPy = (Float_t)pPionAt[1];
          T2PionPz = (Float_t)pPionAt[2];
          T2PionDCAtoPVXY = (Float_t)PionImpactXY;
          T2PionDCAtoPVZ = (Float_t)PionImpactZ;
          T2PionDCAtoSVXY = (Float_t)PionTrack->GetD(SVposition[0],SVposition[1],lMagneticField);
          T2PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
          T2PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
          PiDDecayTree->Fill();
        }
        else if(DeuteronTrack->Charge() <0 && kEnablePiDChannel4)
        {
          fhPiDChargeCheckC4->Fill(0.5,PionTrack->Charge());
          fhPiDChargeCheckC4->Fill(1.5,DeuteronTrack->Charge());
          if(InvmassV0 < V0Invmassmin2 || InvmassV0 > V0Invmassmax2)continue;
          T2DecayChannel = (Int_t)4;
          T2EventPercentile = (Float_t)percentile;
          T2SVDCA = (Float_t)Dcadaughters;
          T2SVCPA = (Float_t)CPAV0;
          T2SVtoPV = (Float_t)SVtoPV;
          T2DeuteronPx = (Float_t)pDeuteronAt[0];
          T2DeuteronPy = (Float_t)pDeuteronAt[1];
          T2DeuteronPz = (Float_t)pDeuteronAt[2];
          T2DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
          T2DeuteronDCAtoPVZ = (Float_t)DeuteronImpactZ;
          T2DeuteronDCAtoSVXY = (Float_t)DeuteronTrack->GetD(SVposition[0],SVposition[1],lMagneticField);
          T2DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
          T2DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
          T2PionPx = (Float_t)pPionAt[0];
          T2PionPy = (Float_t)pPionAt[1];
          T2PionPz = (Float_t)pPionAt[2];
          T2PionDCAtoPVXY = (Float_t)PionImpactXY;
          T2PionDCAtoPVZ = (Float_t)PionImpactZ;
          T2PionDCAtoSVXY = (Float_t)PionTrack->GetD(SVposition[0],SVposition[1],lMagneticField);
          T2PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
          T2PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
          PiDDecayTree->Fill();
        }
      }
    }
    else
    {
      for(Int_t j=0; j<Deuteronid.size(); j++)
      {
        AliESDtrack *DeuteronTrack = 0x0;
        Int_t Deuteronidx = Deuteronid[j];
        if(Pionidx == Deuteronidx)continue;
        DeuteronTrack = fESDevent->GetTrack(Deuteronidx);
        if(kOpenTrackImpacCutV0 && DeuteronTrack->GetD(PVposition[0],PVposition[1],lMagneticField) < TrackDCAtoPVcut) continue; 
        AliExternalTrackParam trackInDeuteron(*DeuteronTrack);
        Float_t DeuteronImpactXY;
        Float_t DeuteronImpactZ;
        DeuteronTrack->GetImpactParameters(DeuteronImpactXY,DeuteronImpactZ);
        Double_t xn, xp;
        Double_t Dcadaughters = 0;
        Dcadaughters = DeuteronTrack->GetDCA(PionTrack,lMagneticField,xn,xp);
        if((xn+xp) < 2*V0Radiusmin || (xn+xp) > 2*V0Radiusmax) continue;
        Bool_t corrected=kFALSE;
        if((trackInPion.GetX() > 3.) && (xp < 3.)) {corrected=kTRUE;}
        if((trackInDeuteron.GetX() > 3.) && (xn < 3.)) {corrected=kTRUE;}
        if(corrected && PionTrack->Charge()>0) {Dcadaughters = trackInPion.GetDCA(&trackInDeuteron,lMagneticField,xn,xp);}
        if(Dcadaughters > V0DCAdaughtercut) continue;
        trackInPion.PropagateTo(xp,lMagneticField);
        trackInDeuteron.PropagateTo(xn,lMagneticField);
        AliESDv0 vertex(trackInPion,Pionidx,trackInDeuteron,Deuteronidx);
        Double_t SVtoPV = vertex.GetD(PVposition[0],PVposition[1],PVposition[2]);
        //if(SVtoPV < V0DCAtoPVcut)continue;
        Double_t SVposition[3] = {0., 0., 0.};
        SVposition[0] = vertex.Xv();
        SVposition[1] = vertex.Yv();
        SVposition[2] = vertex.Zv();
        Double_t r2  = SVposition[0]*SVposition[0] + SVposition[1]*SVposition[1];
        if(r2 < V0Radiusmin*V0Radiusmin || r2 > V0Radiusmax*V0Radiusmax) continue;
        Double_t pPionAt[3] = {0., 0., 0.};
        Double_t pDeuteronAt[3] = {0., 0., 0.};
        PionTrack->GetPxPyPzAt(xp,lMagneticField,pPionAt);
        DeuteronTrack->GetPxPyPzAt(xn,lMagneticField,pDeuteronAt); 
        Float_t CPAV0 = vertex.GetV0CosineOfPointingAngle(PVposition[0],PVposition[1],PVposition[2]);
        if(CPAV0 < V0CPAcut)continue;
        TLorentzVector vDeuteron,vPion,vsumV0;
        vDeuteron.SetXYZM(pDeuteronAt[0],pDeuteronAt[1],pDeuteronAt[2],DeuteronMass);
        vPion.SetXYZM(pPionAt[0],pPionAt[1],pPionAt[2],PionMass);
        vsumV0 = vDeuteron + vPion;
        Float_t InvmassV0 = (Float_t)vsumV0.M();
        if(DeuteronTrack->Charge()>0 && kEnablePiDChannel3)
        {
          fhPiDChargeCheckC3->Fill(0.5,PionTrack->Charge());
          fhPiDChargeCheckC3->Fill(1.5,DeuteronTrack->Charge());
          if(InvmassV0 < V0Invmassmin2 || InvmassV0 > V0Invmassmax2)continue;
          T2DecayChannel = (Int_t)3;
          T2EventPercentile = (Float_t)percentile;
          T2SVDCA = (Float_t)Dcadaughters;
          T2SVCPA = (Float_t)CPAV0;
          T2SVtoPV = (Float_t)SVtoPV;
          T2DeuteronPx = (Float_t)pDeuteronAt[0];
          T2DeuteronPy = (Float_t)pDeuteronAt[1];
          T2DeuteronPz = (Float_t)pDeuteronAt[2];
          T2DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
          T2DeuteronDCAtoPVZ = (Float_t)DeuteronImpactZ;
          T2DeuteronDCAtoSVXY = (Float_t)DeuteronTrack->GetD(SVposition[0],SVposition[1],lMagneticField);
          T2DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
          T2DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
          T2PionPx = (Float_t)pPionAt[0];
          T2PionPy = (Float_t)pPionAt[1];
          T2PionPz = (Float_t)pPionAt[2];
          T2PionDCAtoPVXY = (Float_t)PionImpactXY;
          T2PionDCAtoPVZ = (Float_t)PionImpactZ;
          T2PionDCAtoSVXY = (Float_t)PionTrack->GetD(SVposition[0],SVposition[1],lMagneticField);
          T2PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
          T2PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
          PiDDecayTree->Fill();
        }
        else if(DeuteronTrack->Charge() <0 && kEnablePiDChannel2)
        {
          fhPiDChargeCheckC2->Fill(0.5,PionTrack->Charge());
          fhPiDChargeCheckC2->Fill(1.5,DeuteronTrack->Charge());
          if(InvmassV0 < V0Invmassmin || InvmassV0 > V0Invmassmax)continue;
          T2DecayChannel = (Int_t)2;
          T2EventPercentile = (Float_t)percentile;
          T2SVDCA = (Float_t)Dcadaughters;
          T2SVCPA = (Float_t)CPAV0;
          T2SVtoPV = (Float_t)SVtoPV;
          T2DeuteronPx = (Float_t)pDeuteronAt[0];
          T2DeuteronPy = (Float_t)pDeuteronAt[1];
          T2DeuteronPz = (Float_t)pDeuteronAt[2];
          T2DeuteronDCAtoPVXY = (Float_t)DeuteronImpactXY;
          T2DeuteronDCAtoPVZ = (Float_t)DeuteronImpactZ;
          T2DeuteronDCAtoSVXY = (Float_t)DeuteronTrack->GetD(SVposition[0],SVposition[1],lMagneticField);
          T2DeuteronTPCNcls = (Float_t)DeuteronTrack->GetTPCNcls();
          T2DeuteronTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(DeuteronTrack,AliPID::kDeuteron);
          T2PionPx = (Float_t)pPionAt[0];
          T2PionPy = (Float_t)pPionAt[1];
          T2PionPz = (Float_t)pPionAt[2];
          T2PionDCAtoPVXY = (Float_t)PionImpactXY;
          T2PionDCAtoPVZ = (Float_t)PionImpactZ;
          T2PionDCAtoSVXY = (Float_t)PionTrack->GetD(SVposition[0],SVposition[1],lMagneticField);
          T2PionTPCNcls = (Float_t)PionTrack->GetTPCNcls();
          T2PionTPCnSigma = (Float_t)fPIDResponse->NumberOfSigmasTPC(PionTrack,AliPID::kPion);
          PiDDecayTree->Fill();
        }
      }
    }
  }*/
  PostData(1,fOutputList);
  PostData(2,PiKDCascadeDecayTree);
  PostData(3,PiDDecayTree);
}
void AliAnalysisTaskHyperFinder3Body::Terminate(Option_t* option)
{

}
Double_t AliAnalysisTaskHyperFinder3Body::PropagateToDCA(AliESDv0 *v, AliExternalTrackParam *t, Double_t b,Double_t *xb) {
  //--------------------------------------------------------------------
  // This function returns the DCA between the V0 and the track
  //--------------------------------------------------------------------
  Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
  Double_t r[3]; t->GetXYZ(r);
  Double_t x1=r[0], y1=r[1], z1=r[2];
  Double_t p[3]; t->GetPxPyPz(p);
  Double_t px1=p[0], py1=p[1], pz1=p[2];
  
  Double_t x2,y2,z2;     // position and momentum of V0
  Double_t px2,py2,pz2;
  
  v->GetXYZ(x2,y2,z2);
  v->GetPxPyPz(px2,py2,pz2);
  
  // calculation dca
  
  Double_t dd= Det(x2-x1,y2-y1,z2-z1,px1,py1,pz1,px2,py2,pz2);
  Double_t ax= Det(py1,pz1,py2,pz2);
  Double_t ay=-Det(px1,pz1,px2,pz2);
  Double_t az= Det(px1,py1,px2,py2);
  
  Double_t dca=TMath::Abs(dd)/TMath::Sqrt(ax*ax + ay*ay + az*az);
  if (dca > 2.0) return 1.e+33;
  
  //points of the DCA
  Double_t t1 = Det(x2-x1,y2-y1,z2-z1,px2,py2,pz2,ax,ay,az)/
  Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);
  
  x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;

  
  if (x1*x1+y1*y1 > 105*105) return 1.e+33;
  
  //propagate track to the points of DCA
  
  x1=x1*cs1 + y1*sn1;
  *xb = x1;
  
  if (!t->PropagateTo(x1,b)) {
    fhPropagationcheck->Fill(0.5);
    AliError("Propagation failed");
    //    AliErrorF("Propagation failed for X=%f | V0: %f %f %f",x1,x2,y2,z2);
    //    t->Print();
    //
    return 1.e+33;
  }
  fhPropagationcheck->Fill(1.5);
  return dca;
}
Double_t AliAnalysisTaskHyperFinder3Body::Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 2x2 determinant
  //--------------------------------------------------------------------
  return a00*a11 - a01*a10;
}

Double_t AliAnalysisTaskHyperFinder3Body::Det(Double_t a00,Double_t a01,Double_t a02,
                                 Double_t a10,Double_t a11,Double_t a12,
                                 Double_t a20,Double_t a21,Double_t a22) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 3x3 determinant
  //--------------------------------------------------------------------
  return  a00*Det(a11,a12,a21,a22)-a01*Det(a10,a12,a20,a22)+a02*Det(a10,a11,a20,a21);
}