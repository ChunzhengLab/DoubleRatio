#ifndef AliAnalysisTaskRdoubleRatio_cxx
#define AliAnalysisTaskRdoubleRatio_cxx

//class TList;
//class TH1F;
//class TH2F;
//class TProfile;
//class AliAnalysisUtils;

//#include "AliAnalysisTaskSE.h"

using std::cout;
using std::endl;

class AliAnalysisTaskRdoubleRatio : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskRdoubleRatio();
  AliAnalysisTaskRdoubleRatio(const char *name);
  virtual ~AliAnalysisTaskRdoubleRatio();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  virtual void SetHarmonic(double harmonic) { mHarmonic = harmonic; }
  virtual void SetFilterbit(unsigned int filterbit) { fFilterBit = filterbit; }
  virtual void SetPtMin(double ptMin) { fPtMin = ptMin; }
  virtual void SetPtMax(double ptMax) { fPtMax = ptMax; }
  virtual void SetEtaMax(double etaMax) { fEtaMax = etaMax; }
  virtual void SetNhitsMin(double nhitsMin) { fNhitsMin = nhitsMin; }
  virtual void SetChi2Max(double chi2Max) { fChi2Max = chi2Max; }
  virtual void SetDeDxMin(double deDxMin) { fDeDxMin = deDxMin; }
  virtual void SetDcaXyMax(double dcaXyMax) { fDcaXyMax = dcaXyMax; }
  virtual void SetDcaZMax(double dcaZMax) { fDcaZMax = dcaZMax; }

private:
  double mHarmonic;
  unsigned int fFilterBit;
  double fPtMin;
  double fPtMax;
  double fEtaMax;
  double fNhitsMin;
  double fChi2Max;
  double fDeDxMin;
  double fDcaXyMax;
  double fDcaZMax;

  int GetRunNumBin(int runNum);
  bool IsGoodV0(AliAODv0 *v0);
  bool IsGoodDaughterTrack(const AliAODTrack *t);

  int runNum;
  int runNumBin;
  double vtx[3];
  int vzBin;
  double cent;
  int centBin;

  TList *mOutputList;
  //Event-wise
  TH1I *hEvtCount;
  TH1I *hRunNumBin;
  TH1D *hCent;
  TH2D *hCentCorr[2];
  TH2D *hVxy[2];
  TH1D *hVz[2];

  //Random Event Plane
  TH3D *hPsi2RDM;
  TH3D *hPsi3RDM;
  //TPC Event Plane
  TH3D *hPsi2TPC;
  TH3D *hPsi3TPC;

  // Track-wise
  TH1D *hPt[2];
  TH1D *hEta[2];
  TH1D *hPhi[2];
  TH1D *hNhits[2];
  TH1D *hDcaXy[2];
  TH1D *hDcaZ[2];
  TH2D *hPDedx;

  //CME
  //N(S) (150,-1.5,1.5)
  TH1D* hNDeltaSinSPsi2[10];
  //N(S_sf)
  TH1D* hNDeltaSinSPsi2_sf[10];
  //N(SVert)
  TH1D* hNDeltaSinSVertPsi2[10];
  //N(SVert_sf)
  TH1D* hNDeltaSinSVertPsi2_sf[10];
  //N(S) (150,-1.5,1.5)
  TH1D* hNDeltaSinSPsi3[10];
  //N(S_sf)
  TH1D* hNDeltaSinSPsi3_sf[10];
  //N(SVert)
  TH1D* hNDeltaSinSVertPsi3[10];
  //N(SVert_sf)
  TH1D* hNDeltaSinSVertPsi3_sf[10];

  //CMW
  //N(S) (150,-1.5,1.5)
  TH1D* hNDeltaCosSPsi2[10];
  //N(S_sf)
  TH1D* hNDeltaCosSPsi2_sf[10];
  //N(SVert)
  TH1D* hNDeltaCosSVertPsi2[10];
  //N(SVert_sf)
  TH1D* hNDeltaCosSVertPsi2_sf[10];
  //N(S) (150,-1.5,1.5)
  TH1D* hNDeltaCosSPsi3[10];
  //N(S_sf)
  TH1D* hNDeltaCosSPsi3_sf[10];
  //N(SVert)
  TH1D* hNDeltaCosSVertPsi3[10];
  //N(SVert_sf)
  TH1D* hNDeltaCosSVertPsi3_sf[10];

  AliAnalysisTaskRdoubleRatio(const AliAnalysisTaskRdoubleRatio &);
  AliAnalysisTaskRdoubleRatio &operator=(const AliAnalysisTaskRdoubleRatio &);

  ClassDef(AliAnalysisTaskRdoubleRatio, 1);
};
#endif