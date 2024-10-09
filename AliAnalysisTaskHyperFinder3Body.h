#ifndef AliAnalysisTaskHyperFinder3Body_H
#define AliAnalysisTaskHyperFinder3Body_H
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "TList.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliEventCuts.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
class AliAnalysisTaskHyperFinder3Body : public AliAnalysisTaskSE
{
  public:
    AliAnalysisTaskHyperFinder3Body();
    AliAnalysisTaskHyperFinder3Body(const char *name);
    virtual ~AliAnalysisTaskHyperFinder3Body();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    Double_t Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const;
    Double_t Det(Double_t a00,Double_t a01,Double_t a02,
                 Double_t a10,Double_t a11,Double_t a12,
                 Double_t a20,Double_t a21,Double_t a22) const;
    Double_t PropagateToDCA(AliESDv0 *vtx,AliExternalTrackParam *trk,Double_t b, Double_t* xb);
    const Double_t PionMass = 0.13957018;
    const Double_t KaonMass = 0.493677;
    const Double_t DeuteronMass = 1.875612;
    const Double_t ProtonMass = 0.938272;
    const Double_t LambdaMass = 1.115683;
    const Double_t XiMass = 1.32171;
    const Double_t OmegaMass = 1.67245;
    const Double_t PhiMass = 1.019461;
  private:
  AliInputEventHandler *fInputHandler;
  AliESDEvent *fESDevent;
  TList *fOutputList;
  AliESDtrackCuts *fESDtrackCuts;
  AliEventCuts feventCut;
  AliPIDResponse *fPIDResponse;
  Int_t nRotation;
  Double_t TrackDCAtoPVcut;
  Bool_t TrackTPCRefit;
  Bool_t TrackAcceptKinkDaughters;
  Double_t TrackTPCNClscut;
  Double_t TrackLengthcut;
  Double_t TrackNCrossedRowscut;
  Double_t TrackEtacut;
  Double_t TrackChi2PerClusterTPCcut;
  Double_t V0Radiusmin;
  Double_t V0Radiusmax;
  Double_t V0DCAdaughtercut;
  Double_t V0CPAcut;
  Double_t V0DCAtoPVcut;
  Double_t CascadeDCAcut;
  Double_t CascadeV0CPAcut;
  Double_t CascadeCPAcut;
  Double_t CascadeRadiusmin;
  Double_t CascadeRadiusmax;
  Double_t Cascadeptmin;
  Double_t Cascadeptmax;
  Double_t CascadeV0Invmassmax;
  Double_t CascadeInvmassmin;
  Double_t CascadeInvmassmax;
  Double_t CascadeInvmassmin2;
  Double_t CascadeInvmassmax2;
  Double_t V0Invmassmin;
  Double_t V0Invmassmax;
  Double_t V0Invmassmin2;
  Double_t V0Invmassmax2;
  Bool_t kOpenTrackImpacCutV0;
  Bool_t kRotatePion;
  Bool_t kRotateKaon;
  Bool_t kRotateDeuteron;

  TH1I* fhEventCount;
  TH1F* fhcentralitybeforecut;
  TH1F* fhcentralityaftercut;
  TH1I* fhCausalitycheck;
  TH1D* fhTrackDCAtoPVxycheck;
  TH2I* fhPiKDChargeCheckC1;
  TH2I* fhPiKDChargeCheckC2;
  TH2I* fhPiKDChargeCheckC3;
  TH2I* fhPiKDChargeCheckC4;
  TH2I* fhPiKDChargeCheckC5;
  TH2I* fhPiKDChargeCheckC6;
  TH2I* fhPiKDChargeCheckC7;
  TH2I* fhPiKDChargeCheckC8;
  TH2I* fhPiDChargeCheckC1;
  TH2I* fhPiDChargeCheckC2;
  TH2I* fhPiDChargeCheckC3;
  TH2I* fhPiDChargeCheckC4;
  TH1I* fhPropagationcheck;
  TH2D* fhCascadePropagationPcheck;
  TH1I* fhDeuteronEventCheck;
  TH1F* fhTrigger;
  TH2D* fhBBPion;
  TH2D* fhBBDeuteron;
  TH2D* fhBBKaon;

  TTree * PiKDCascadeDecayTree;
  Int_t   T1DecayChannel;
  Float_t T1EventPercentile;
  Float_t T1TVDCA;
  Float_t T1TVCPA;
  Float_t T1TVtoPV;
  Float_t T1SVDCA;
  Float_t T1SVCPA;
  Float_t T1SVtoPV;
  Float_t T1SVtoTV;
  Float_t T1KaonPx;
  Float_t T1KaonPy;
  Float_t T1KaonPz;
  Float_t T1KaonDCAtoPVXY;
  Float_t T1KaonTPCNcls;
  Float_t T1KaonTPCnSigma;
  Float_t T1DeuteronPx;
  Float_t T1DeuteronPy;
  Float_t T1DeuteronPz;
  Float_t T1DeuteronDCAtoPVXY;
  Float_t T1DeuteronTPCNcls;
  Float_t T1DeuteronTPCnSigma;
  Float_t T1PionPx;
  Float_t T1PionPy;
  Float_t T1PionPz;
  Float_t T1PionDCAtoPVXY;
  Float_t T1PionTPCNcls;
  Float_t T1PionTPCnSigma;
  
  TTree * PiDDecayTree;
  Int_t   T2DecayChannel;
  Float_t T2EventPercentile;
  Float_t T2SVDCA;
  Float_t T2SVCPA;
  Float_t T2SVtoPV;
  Float_t T2DeuteronPx;
  Float_t T2DeuteronPy;
  Float_t T2DeuteronPz;
  Float_t T2DeuteronDCAtoPVXY;
  Float_t T2DeuteronDCAtoPVZ;
  Float_t T2DeuteronDCAtoSVXY;
  Float_t T2DeuteronTPCNcls;
  Float_t T2DeuteronTPCnSigma;
  Float_t T2PionPx;
  Float_t T2PionPy;
  Float_t T2PionPz;
  Float_t T2PionDCAtoPVXY;
  Float_t T2PionDCAtoPVZ;
  Float_t T2PionDCAtoSVXY;
  Float_t T2PionTPCNcls;
  Float_t T2PionTPCnSigma;


  AliAnalysisTaskHyperFinder3Body(const AliAnalysisTaskHyperFinder3Body&); // not implemented
  AliAnalysisTaskHyperFinder3Body& operator=(const AliAnalysisTaskHyperFinder3Body&); // not implemented
  ClassDef(AliAnalysisTaskHyperFinder3Body, 1); 
};
#endif