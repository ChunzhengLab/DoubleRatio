#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <vector>
#include <algorithm>
// ROOT classes
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1I.h"
#include "TH3I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TExMap.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TComplex.h"
// Alice analysis base class
#include "AliAnalysisTaskSE.h"
// Alice analysis additional classes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
// Alice AOD classes
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"

#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODHeader.h"

// Alice classes
#include "AliCentrality.h"
#include "AliEventplane.h"
//#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
// Alice MC classes
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
// Alice "V" classes
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVVZERO.h"
// Alice PID classes
#include "AliAODPid.h"
#include "AliAODpidUtil.h"
#include "AliPID.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
//#include "AliMultSelection.h"

#include "AliAnalysisTaskRdoubleRatio.h"

using std::cout;
using std::endl;
using std::vector;

ClassImp(AliAnalysisTaskRdoubleRatio);

//---------------------------------------------------
AliAnalysisTaskRdoubleRatio::AliAnalysisTaskRdoubleRatio() : AliAnalysisTaskSE(),
mHarmonic(2.),
fFilterBit(1),
fPtMin(0.2),
fPtMax(5.),
fEtaMax(0.8),
fNhitsMin(80),
fChi2Max(4.),
fDeDxMin(10),
fDcaXyMax(3.),
fDcaZMax(3.)
{
  runNum       = -999;
  runNumBin    = -999;
  for (int i=0; i<3; ++i) vtx[i] = -999;
  vzBin        = -999;
  cent         = -999;
  centBin      = -999;

  mOutputList = NULL;
  // Event-wise
  hEvtCount  = NULL;
  hRunNumBin = NULL;
  hCent      = NULL;
  for (int i=0; i<2; ++i) hCentCorr[i] = NULL;
  for (int i=0; i<2; ++i) hVxy[i]      = NULL;
  for (int i=0; i<2; ++i) hVz[i]       = NULL;

  //Random Event Plane
  hPsi2RDM = NULL;
  hPsi3RDM = NULL;

  //TPC Event Plane
  hPsi2TPC = NULL;
  hPsi3TPC = NULL;

  // Track-wise
  for (int i=0; i<2; ++i) hPt[i]    = NULL;
  for (int i=0; i<2; ++i) hEta[i]   = NULL;
  for (int i=0; i<2; ++i) hPhi[i]   = NULL;
  for (int i=0; i<2; ++i) hNhits[i] = NULL;
  for (int i=0; i<2; ++i) hDcaXy[i] = NULL;
  for (int i=0; i<2; ++i) hDcaZ[i]  = NULL;
  hPDedx = NULL;

  //CME
  //N(S)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi2[i] = NULL;
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi2_sf[i] = NULL;
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi2[i] = NULL;
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi2_sf[i] = NULL;
  //N(S)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi3[i] = NULL;
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi3_sf[i] = NULL;
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi3[i] = NULL;
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi3_sf[i] = NULL;

  //CMW
  //N(S)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi2[i] = NULL;
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi2_sf[i] = NULL;
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi2[i] = NULL;
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi2_sf[i] = NULL;
  //N(S)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi3[i] = NULL;
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi3_sf[i] = NULL;
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi3[i] = NULL;
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi3_sf[i] = NULL;
}

//---------------------------------------------------
AliAnalysisTaskRdoubleRatio::AliAnalysisTaskRdoubleRatio(const char *name) : AliAnalysisTaskSE(name),
mHarmonic(2.),
fFilterBit(1),
fPtMin(0.2),
fPtMax(5.),
fEtaMax(0.8),
fNhitsMin(80),
fChi2Max(4.),
fDeDxMin(10),
fDcaXyMax(3.),
fDcaZMax(3.)
{
  runNum       = -999;
  runNumBin    = -999;
  for (int i=0; i<3; ++i) vtx[i] = -999;
  vzBin        = -999;
  cent         = -999;
  centBin      = -999;

  mOutputList = NULL;
  // Event-wise
  hEvtCount  = NULL;
  hRunNumBin = NULL;
  hCent      = NULL;
  for (int i=0; i<2; ++i) hCentCorr[i] = NULL;
  for (int i=0; i<2; ++i) hVxy[i]      = NULL;
  for (int i=0; i<2; ++i) hVz[i]       = NULL;

  //Random Event Plane
  hPsi2RDM = NULL;
  hPsi3RDM = NULL;

  //TPC Event Plane
  hPsi2TPC = NULL;
  hPsi3TPC = NULL;

  // Track-wise
  for (int i=0; i<2; ++i) hPt[i]    = NULL;
  for (int i=0; i<2; ++i) hEta[i]   = NULL;
  for (int i=0; i<2; ++i) hPhi[i]   = NULL;
  for (int i=0; i<2; ++i) hNhits[i] = NULL;
  for (int i=0; i<2; ++i) hDcaXy[i] = NULL;
  for (int i=0; i<2; ++i) hDcaZ[i]  = NULL;
  hPDedx = NULL;

  //CME
  //N(S)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi2[i] = NULL;
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi2_sf[i] = NULL;
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi2[i] = NULL;
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi2_sf[i] = NULL;
  //N(S)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi3[i] = NULL;
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi3_sf[i] = NULL;
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi3[i] = NULL;
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi3_sf[i] = NULL;

  //CMW
  //N(S)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi2[i] = NULL;
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi2_sf[i] = NULL;
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi2[i] = NULL;
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi2_sf[i] = NULL;
  //N(S)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi3[i] = NULL;
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi3_sf[i] = NULL;
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi3[i] = NULL;
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi3_sf[i] = NULL;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//---------------------------------------------------
AliAnalysisTaskRdoubleRatio::~AliAnalysisTaskRdoubleRatio()
{
}

//---------------------------------------------------
void AliAnalysisTaskRdoubleRatio::Terminate(Option_t *) 
{
}

//---------------------------------------------------
void AliAnalysisTaskRdoubleRatio::UserCreateOutputObjects()
{
  mOutputList = new TList();
  mOutputList->SetName(GetName());
  mOutputList->SetOwner(kTRUE);

  // event-wise
  hEvtCount = new TH1I("evtCount","",20, 1, 21);
  hEvtCount->GetXaxis()->SetBinLabel(1,"All");
  hEvtCount->GetXaxis()->SetBinLabel(2,"Info");
  hEvtCount->GetXaxis()->SetBinLabel(3,"Evt");
  hEvtCount->GetXaxis()->SetBinLabel(4,"Cent");

  hEvtCount->GetXaxis()->SetBinLabel(10,"Manager");
  hEvtCount->GetXaxis()->SetBinLabel(11,"Handler");
  hEvtCount->GetXaxis()->SetBinLabel(12,"AOD");
  hEvtCount->GetXaxis()->SetBinLabel(13,"PID");
  hEvtCount->GetXaxis()->SetBinLabel(14,"Utils");
  // hEvtCount->GetXaxis()->SetBinLabel(15,"MultSel");

  hEvtCount->GetXaxis()->SetBinLabel(16,"Random plane");
  hEvtCount->GetXaxis()->SetBinLabel(17,"TPC plane");
  hEvtCount->GetXaxis()->SetBinLabel(18,"VZERO plane");
  hEvtCount->GetXaxis()->SetBinLabel(19,"ZDC plane");
  hEvtCount->GetXaxis()->SetBinLabel(20,"Loops end");
  mOutputList->Add(hEvtCount);

  // 10h
  TString runNumList[91]={
    "139510","139507","139505","139503","139465","139438","139437","139360","139329","139328","139314","139310",
    "139309","139173","139107","139105","139038","139037","139036","139029","139028","138872","138871","138870",
    "138837","138732","138730","138666","138662","138653","138652","138638","138624","138621","138583","138582",
    "138579","138578","138534","138469","138442","138439","138438","138396","138364","138275","138225","138201",
    "138197","138192","138190","137848","137844","137752","137751","137724","137722","137718","137704","137693",
    "137692","137691","137686","137685","137639","137638","137608","137595","137549","137546","137544","137541",
    "137539","137531","137530","137443","137441","137440","137439","137434","137432","137431","137430","137243",
    "137236","137235","137232","137231","137230","137162","137161"
  };
  // 11h
  // TString runNum[39]={"170387","170040","170268","170228","170207","169838","170159","170204","170311","170084",
  //                     "169835","170088","170593","170203","170270","169846","170163","170388","170155","170083",
  //                     "170572","169837","169855","170306","170269","170089","170309","170091","170081","170230",
  //                     "170085","170315","170027","170193","170312","170313","170308","169858","169859"};
  hRunNumBin = new TH1I("runNumBin","",100,0,100);
  for (int i=0; i<91; ++i) {    
    hRunNumBin->GetXaxis()->SetBinLabel(i+1,runNumList[i].Data());
  }
  mOutputList->Add(hRunNumBin);

  hCent = new TH1D("centrality","",100,0,100);
  mOutputList->Add(hCent);
  hCentCorr[0] = new TH2D("centcorr0","",100,0,100,100,0,100);
  hCentCorr[1] = new TH2D("centcorr1","",100,0,100,100,0,100);
  for (int i=0; i<2; ++i) mOutputList->Add(hCentCorr[i]);

  hVxy[0] = new TH2D("vxy0","",100,-0.5,0.5,100,-0.5,0.5);
  hVxy[1] = new TH2D("vxy1","",100,-0.5,0.5,100,-0.5,0.5);
  hVz[0]  = new TH1D("vz0","",200,-50,50);
  hVz[1]  = new TH1D("vz1","",200,-50,50);
  for (int i=0; i<2; ++i) mOutputList->Add(hVxy[i]);
  for (int i=0; i<2; ++i) mOutputList->Add(hVz[i]);


  //Random Event Plane
  hPsi2RDM = new TH3D("hPsiRDM","",100, 0, 100, 10, 0, 10, 180, 0, TMath::Pi());
  hPsi3RDM = new TH3D("hPsiRDM","",100, 0, 100, 10, 0, 10, 180, 0, 2./3.*TMath::Pi());
  mOutputList->Add(hPsi2RDM);
  mOutputList->Add(hPsi3RDM);

  //TPC Event Plane
  hPsi2TPC = new TH3D("hPsiTPC","",100, 0, 100, 10, 0, 10, 180, 0, TMath::Pi());
  hPsi3TPC = new TH3D("hPsiTPC","",100, 0, 100, 10, 0, 10, 180, 0, 2./3.*TMath::Pi());
  mOutputList->Add(hPsi2TPC);
  mOutputList->Add(hPsi3TPC);

  // track-wise
  hPt[0] = new TH1D("hPtBeforeCut", "", 200, 0., 20.);
  hPt[1] = new TH1D("hPtAfterCut", "", 200, 0., 20.);
  for (int i=0; i<2; ++i) mOutputList->Add(hPt[i]);
  hEta[0] = new TH1D("hEtaBeforeCut", "", 200, -10., 10.);
  hEta[1] = new TH1D("hEtaAfterCut",  "", 200, -10., 10.);
  for (int i=0; i<2; ++i) mOutputList->Add(hEta[i]);
  hPhi[0] = new TH1D("hPhiBeforeCut", "", 400, 0, 2*TMath::Pi());
  hPhi[1] = new TH1D("hPhiAfterCut", "", 400, 0, 2*TMath::Pi());
  for (int i=0; i<2; ++i) mOutputList->Add(hPhi[i]);
  hNhits[0] = new TH1D("hNhitsBeforeCut", "", 200, 0., 200.);
  hNhits[1] = new TH1D("hNhitsAfterCut",  "", 200, 0., 200.);
  for (int i=0; i<2; ++i) mOutputList->Add(hNhits[i]);
  hDcaXy[0] = new TH1D("hDcaXyBeforeCut", "", 100, 0., 10.);
  hDcaXy[1] = new TH1D("hDcaXyAfterCut",  "", 100, 0., 10.);
  for (int i=0; i<2; ++i) mOutputList->Add(hDcaXy[i]);
  hDcaZ[0] = new TH1D("hDcaZBeforeCut", "", 100, 0., 10.);
  hDcaZ[1] = new TH1D("hDcaZAfterCut",  "", 100, 0., 10.);
  for (int i=0; i<2; ++i) mOutputList->Add(hDcaZ[i]);
  hPDedx = new TH2D("hPDedx", "", 400, -10., 10., 400, 0, 1000);
  mOutputList->Add(hPDedx);

  //CME
  //N(S)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi2[i] = new TH1D(Form("hNDeltaSinSPsi2_cent%i",i),"",150,-1.5,1.5);
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi2_sf[i] = new TH1D(Form("hNDeltaSinSPsi2_sf_cent%i",i),"",150,-1.5,1.5);
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi2[i] = new TH1D(Form("hNDeltaSinSVertPsi2_cent%i",i),"",150,-1.5,1.5);
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi2_sf[i] = new TH1D(Form("hNDeltaSinSVertPsi2_sf_cent%i",i),"",150,-1.5,1.5);

  //N(S)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi3[i] = new TH1D(Form("hNDeltaSinSPsi3_cent%i",i),"",150,-1.5,1.5);
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSPsi3_sf[i] = new TH1D(Form("hNDeltaSinSPsi3_sf_cent%i",i),"",150,-1.5,1.5);
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi3[i] = new TH1D(Form("hNDeltaSinSVertPsi3_cent%i",i),"",150,-1.5,1.5);
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaSinSVertPsi3_sf[i] = new TH1D(Form("hNDeltaSinSvertPsi3_sf_cent%i",i),"",150,-1.5,1.5);

  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSPsi2[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSPsi2_sf[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSVertPsi2[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSVertPsi2_sf[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSPsi3[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSPsi3_sf[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSVertPsi3[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSVertPsi3_sf[i]);

  //CMW
  //N(S)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi2[i] = new TH1D(Form("hNDeltaCosSPsi2_cent%i",i),"",150,-1.5,1.5);
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi2_sf[i] = new TH1D(Form("hNDeltaCosSPsi2_sf_cent%i",i),"",150,-1.5,1.5);
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi2[i] = new TH1D(Form("hNDeltaCosSVertPsi2_cent%i",i),"",150,-1.5,1.5);
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi2_sf[i] = new TH1D(Form("hNDeltaCosSVertPsi2_sf_cent%i",i),"",150,-1.5,1.5);

  //N(S)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi3[i] = new TH1D(Form("hNDeltaCosSPsi3_cent%i",i),"",150,-1.5,1.5);
  //N(S_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSPsi3_sf[i] = new TH1D(Form("hNDeltaCosSPsi3_sf_cent%i",i),"",150,-1.5,1.5);
  //N(SVert)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi3[i] = new TH1D(Form("hNDeltaCosSVertPsi3_cent%i",i),"",150,-1.5,1.5);
  //N(SVert_sf)
  for (int i=0; i<10; ++i) hNDeltaCosSVertPsi3_sf[i] = new TH1D(Form("hNDeltaCosSvertPsi3_sf_cent%i",i),"",150,-1.5,1.5);

  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSPsi2[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSPsi2_sf[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSVertPsi2[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSVertPsi2_sf[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSPsi3[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSPsi3_sf[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSVertPsi3[i]);
  for (int i=0; i<10; ++i) mOutputList->Add(hNDeltaCosSVertPsi3_sf[i]);

  PostData(1,mOutputList);
}

//---------------------------------------------------
void AliAnalysisTaskRdoubleRatio::UserExec(Option_t *)
{
  hEvtCount->Fill(1);

  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else hEvtCount->Fill(10);
  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else hEvtCount->Fill(11);
  AliAODEvent* fAOD = (AliAODEvent*)InputEvent();
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  } else hEvtCount->Fill(12);
  AliPIDResponse* fPID = handler->GetPIDResponse();
  if (!fPID) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else hEvtCount->Fill(13);
  AliAnalysisUtils* fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else hEvtCount->Fill(14);
  // AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
  // if (!fMultSel) {
  //   AliError(Form("%s: Could not get AliMultSelection", GetName()));
  // } else hEvtCount->Fill(15);
  if (!manager || !handler || !fAOD || !fUtils) return;
  hEvtCount->Fill(2);

  //------------------
  // event-wise
  //------------------

  // runNumber
  runNum = fAOD->GetRunNumber();
  runNumBin = GetRunNumBin(runNum);
  if (runNumBin<0) return;
  hRunNumBin->Fill(runNumBin);

  // vertex
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  vtx[0] = (double)fVtx->GetX();
  vtx[1] = (double)fVtx->GetY();
  vtx[2] = (double)fVtx->GetZ();
  double vzSPD  = fAOD->GetPrimaryVertexSPD()->GetZ();
  hVxy[0]->Fill(vtx[0], vtx[1]);
  hVz[0]->Fill(vtx[2]);
  if (fabs(vtx[0])<1e-6 || fabs(vtx[1])<1e-6 || fabs(vtx[2])<1e-6) return;


  if (fabs(vtx[2])>10) return;
  if (fabs(vtx[2]-vzSPD)>0.5) return;
  hVxy[1]->Fill(vtx[0], vtx[1]);
  hVz[1]->Fill(vtx[2]);
  for (int i = 0; i < 20; ++i) {
    if (vtx[2] > -10+i*1 && vtx[2] < -10+(i+1)*1) {vzBin = i; break;}
  }
  if (vzBin<0) return;

  // pileup
  fUtils->SetUseOutOfBunchPileUp(true);
  fUtils->SetUseMVPlpSelection(true);
  // fUtils->SetMinPlpContribMV(5);
  bool isPileup = fUtils->IsPileUpEvent(fAOD);
  // bool isPileup = fUtils->IsPileUpMV(fAOD); // pp, p-Pb
  if (isPileup) return;
  hEvtCount->Fill(3); 

  // centrality
  // run1
  cent = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
  double cent2 = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
  double cent3 = fAOD->GetCentrality()->GetCentralityPercentile("TRK");
  hCentCorr[0]->Fill(cent,cent2);
  if (fabs(cent-cent2)>7.5) return;
  if (fabs(cent-cent3)>5) return; // ANA-280
  if (cent<0 || cent>=100) return;
  hCentCorr[1]->Fill(cent,cent2);
  centBin = (int)cent/10;
  hCent->Fill(cent);
  hEvtCount->Fill(4);
  
  //------------------
  //* loop trk
  //------------------
  vector<double> vecPhi;
  vector<int> vecCharge;

  //TPC QxQy
  double sumCos2Phi = 0.;
  double sumSin2Phi = 0.;
  double sumCos3Phi = 0.;
  double sumSin3Phi = 0.;

  int nTrk = fAOD->GetNumberOfTracks();
  for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrk);
    if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }
    if (!track->TestFilterBit(fFilterBit)) continue;
    double pt     = track->Pt();
    double eta    = track->Eta();
    double phi    = track->Phi();
    int    charge = track->Charge();
    int    nhits  = track->GetTPCNcls();
    double dedx   = track->GetTPCsignal();
    double chi2   = track->Chi2perNDF();
    hPt[0]->Fill(pt);
    hEta[0]->Fill(eta);
    hPhi[0]->Fill(phi);
    hNhits[0]->Fill(nhits);

    if (pt<fPtMin || pt>fPtMax) continue;
    if (fabs(eta)>fEtaMax) continue;
    if (fabs(nhits)<fNhitsMin) continue;
    if (chi2>fChi2Max) continue;
    if (dedx<fDeDxMin) continue;

    hPt[1]->Fill(pt);
    hEta[1]->Fill(eta);
    hPhi[1]->Fill(phi);
    hNhits[1]->Fill(nhits);
    hPDedx->Fill(track->P()*charge, dedx);

    //------------------
    // dca
    //------------------
    double mag = fAOD->GetMagneticField();
    double dcaxy  = 999.;
    double dcaz   = 999.;
    double r[3];
    double dca[2];
    double cov[3];
    bool proptodca = track->PropagateToDCA(fVtx, mag, 100., dca, cov);
    if (track->GetXYZ(r)) {
      dcaxy = r[0];
      dcaz  = r[1];
    } else {
      double dcax = r[0] - vtx[0];
      double dcay = r[1] - vtx[1];
      dcaz  = r[2] - vtx[2];
      dcaxy = sqrt(dcax*dcax + dcay*dcay);
      // dcaxy = dca[0];
    }
    hDcaXy[0]->Fill(dcaxy);
    if (fabs(dcaxy)>fDcaXyMax) continue;//DCAxy cut
    hDcaXy[1]->Fill(dcaxy);
    hDcaZ[0]->Fill(dcaz);
    if (fabs(dcaz)>fDcaZMax) continue;//DCAz cut
    hDcaZ[1]->Fill(dcaz);

    vecPhi.push_back(phi);
    vecCharge.push_back(charge);

    //Qn
    double cos2phi = TMath::Cos(2 * phi);
    double sin2phi = TMath::Sin(2 * phi);
    double cos3phi = TMath::Cos(3 * phi);
    double sin3phi = TMath::Sin(3 * phi);
    sumCos2Phi += cos2phi;
    sumSin2Phi += sin2phi;
    sumCos3Phi += cos3phi;
    sumSin3Phi += sin3phi;
  }

  //RDM Plane
  double psi2RDM = gRandom->Uniform(0,TMath::Pi());
  hPsi2RDM -> Fill(runNumBin, centBin, psi2RDM);
  double psi3RDM = gRandom->Uniform(0,2/3*TMath::Pi());
  hPsi2RDM -> Fill(runNumBin, centBin, psi3RDM);
  hEvtCount->Fill(16);
  
  //TPC Plane
  TVector2 Q2TPC;
  Q2TPC.Set(sumCos2Phi,sumSin2Phi);
  double psi2TPC = Q2TPC.Phi()/2.;
  hPsi2TPC -> Fill(runNumBin, centBin, psi2TPC);
  TVector2 Q3TPC;
  Q3TPC.Set(sumCos3Phi,sumSin3Phi);
  double psi3TPC = Q3TPC.Phi()/3.;
  hPsi3TPC -> Fill(runNumBin, centBin, psi3TPC);
  hEvtCount->Fill(17);

  vector<int> vecCharge_sf;
  vecCharge_sf.assign(vecCharge.begin(),vecCharge.end());
  random_shuffle(vecCharge_sf.begin(), vecCharge_sf.end());

  //Sin
  double sumSin2DeltaPhiPos = 0.;
  double sumSin2DeltaPhiPos_sf = 0.;
  double sumSin2DeltaPhiNeg = 0.;
  double sumSin2DeltaPhiNeg_sf = 0.;

  double sumSin3DeltaPhiPos = 0.;
  double sumSin3DeltaPhiPos_sf = 0.;
  double sumSin3DeltaPhiNeg = 0.;
  double sumSin3DeltaPhiNeg_sf = 0.;

  double sumSin2DeltaPhiVertPos = 0.;
  double sumSin2DeltaPhiVertPos_sf = 0.;
  double sumSin2DeltaPhiVertNeg = 0.;
  double sumSin2DeltaPhiVertNeg_sf = 0.;

  double sumSin3DeltaPhiVertPos = 0.;
  double sumSin3DeltaPhiVertPos_sf = 0.;
  double sumSin3DeltaPhiVertNeg = 0.;
  double sumSin3DeltaPhiVertNeg_sf = 0.;

  //Cos
  double sumCos2DeltaPhiPos = 0.;
  double sumCos2DeltaPhiPos_sf = 0.;
  double sumCos2DeltaPhiNeg = 0.;
  double sumCos2DeltaPhiNeg_sf = 0.;

  double sumCos3DeltaPhiPos = 0.;
  double sumCos3DeltaPhiPos_sf = 0.;
  double sumCos3DeltaPhiNeg = 0.;
  double sumCos3DeltaPhiNeg_sf = 0.;

  double sumCos2DeltaPhiVertPos = 0.;
  double sumCos2DeltaPhiVertPos_sf = 0.;
  double sumCos2DeltaPhiVertNeg = 0.;
  double sumCos2DeltaPhiVertNeg_sf = 0.;

  double sumCos3DeltaPhiVertPos = 0.;
  double sumCos3DeltaPhiVertPos_sf = 0.;
  double sumCos3DeltaPhiVertNeg = 0.;
  double sumCos3DeltaPhiVertNeg_sf = 0.;

  int nPos = 0;
  int nNeg = 0;

  for (vector<double>::size_type iTrk = 0; iTrk < vecPhi.size(); iTrk++) {
    double phi = vecPhi[iTrk];
    int  charge  = vecCharge[iTrk];
    int  charge_sf = vecCharge_sf[iTrk];

    TVector2 Q2TPC_tmp;
    TVector2 Q3TPC_tmp;
    Q2TPC_tmp.Set(sumCos2Phi - TMath::Cos(2 * phi), sumSin2Phi - TMath::Sin(2 * phi));
    Q3TPC_tmp.Set(sumCos3Phi - TMath::Cos(3 * phi), sumSin3Phi - TMath::Sin(3 * phi));
    double psi2TPC_tmp = Q2TPC_tmp.Phi()/2;
    double psi3TPC_tmp = Q3TPC_tmp.Phi()/3;

    if(charge > 0) 
    {
      sumSin2DeltaPhiPos += TMath::Sin(phi-psi2TPC_tmp);
      sumSin3DeltaPhiPos += TMath::Sin(3/2*(phi-psi3TPC_tmp));
      sumSin2DeltaPhiVertPos += TMath::Sin(phi-psi2TPC_tmp+TMath::Pi()/2);
      sumSin3DeltaPhiVertPos += TMath::Sin(3/2*(phi-psi3TPC_tmp+TMath::Pi()/3));

      sumCos2DeltaPhiPos += TMath::Cos(2*(phi-psi2TPC_tmp));
      sumCos3DeltaPhiPos += TMath::Cos(3*(phi-psi3TPC_tmp));
      sumCos2DeltaPhiVertPos += TMath::Cos(2*(phi-psi2TPC_tmp+TMath::Pi()/2));
      sumCos3DeltaPhiVertPos += TMath::Cos(3*(phi-psi3TPC_tmp+TMath::Pi()/3));
      nPos++;
    }
    if(charge_sf > 0) 
    {
      sumSin2DeltaPhiPos_sf += TMath::Sin(phi-psi2TPC_tmp);
      sumSin3DeltaPhiPos_sf += TMath::Sin(3/2*(phi-psi3TPC_tmp));
      sumSin2DeltaPhiVertPos_sf += TMath::Sin(phi-psi2TPC_tmp+TMath::Pi()/2);
      sumSin3DeltaPhiVertPos_sf += TMath::Sin(3/2*(phi-psi3TPC_tmp+TMath::Pi()/3));

      sumCos2DeltaPhiPos_sf += TMath::Cos(2*(phi-psi2TPC_tmp));
      sumCos3DeltaPhiPos_sf += TMath::Cos(3*(phi-psi3TPC_tmp));
      sumCos2DeltaPhiVertPos_sf += TMath::Cos(2*(phi-psi2TPC_tmp+TMath::Pi()/2));
      sumCos3DeltaPhiVertPos_sf += TMath::Cos(3*(phi-psi3TPC_tmp+TMath::Pi()/3));
    }
    if(charge < 0) 
    {
      sumSin2DeltaPhiNeg += TMath::Sin(phi-psi2TPC_tmp);
      sumSin3DeltaPhiNeg += TMath::Sin(3/2*(phi-psi3TPC_tmp));
      sumSin2DeltaPhiVertNeg += TMath::Sin(phi-psi2TPC_tmp+TMath::Pi()/2);
      sumSin3DeltaPhiVertNeg += TMath::Sin(3/2*(phi-psi3TPC_tmp+TMath::Pi()/3));

      sumCos2DeltaPhiNeg += TMath::Cos(2*(phi-psi2TPC_tmp));
      sumCos3DeltaPhiNeg += TMath::Cos(3*(phi-psi3TPC_tmp));
      sumCos2DeltaPhiVertNeg += TMath::Cos(2*(phi-psi2TPC_tmp+TMath::Pi()/2));
      sumCos3DeltaPhiVertNeg += TMath::Cos(3*(phi-psi3TPC_tmp+TMath::Pi()/3));
      nNeg++;

    }
    if(charge_sf < 0) 
    {
      sumSin2DeltaPhiNeg_sf += TMath::Sin(phi-psi2TPC_tmp);
      sumSin3DeltaPhiNeg_sf += TMath::Sin(3/2*(phi-psi3TPC_tmp));
      sumSin2DeltaPhiVertNeg_sf += TMath::Sin(phi-psi2TPC_tmp+TMath::Pi()/2);
      sumSin3DeltaPhiVertNeg_sf += TMath::Sin(3/2*(phi-psi3TPC_tmp+TMath::Pi()/3));

      sumCos2DeltaPhiNeg_sf += TMath::Cos(2*(phi-psi2TPC_tmp));
      sumCos3DeltaPhiNeg_sf += TMath::Cos(3*(phi-psi3TPC_tmp));
      sumCos2DeltaPhiVertNeg_sf += TMath::Cos(2*(phi-psi2TPC_tmp+TMath::Pi()/2));
      sumCos3DeltaPhiVertNeg_sf += TMath::Cos(3*(phi-psi3TPC_tmp+TMath::Pi()/3));
    }
  }
    if(nPos * nNeg == 0) return;
    //Sin
    sumSin2DeltaPhiPos /= (double)nPos;
    sumSin2DeltaPhiPos_sf /= (double)nPos;
    sumSin2DeltaPhiVertPos /= (double)nPos;
    sumSin2DeltaPhiVertPos_sf /= (double)nPos;

    sumSin2DeltaPhiNeg /= (double)nNeg;
    sumSin2DeltaPhiNeg_sf /= (double)nNeg;
    sumSin2DeltaPhiVertNeg /= (double)nNeg;
    sumSin2DeltaPhiVertNeg_sf /= (double)nNeg;

    sumSin3DeltaPhiPos /= (double)nPos;
    sumSin3DeltaPhiPos_sf /= (double)nPos;
    sumSin3DeltaPhiVertPos /= (double)nPos;
    sumSin3DeltaPhiVertPos_sf /= (double)nPos;

    sumSin3DeltaPhiNeg /= (double)nNeg;
    sumSin3DeltaPhiNeg_sf /= (double)nNeg;
    sumSin3DeltaPhiVertNeg /= (double)nNeg;
    sumSin3DeltaPhiVertNeg_sf /= (double)nNeg;

    //Cos
    sumCos2DeltaPhiPos /= (double)nPos;
    sumCos2DeltaPhiPos_sf /= (double)nPos;
    sumCos2DeltaPhiVertPos /= (double)nPos;
    sumCos2DeltaPhiVertPos_sf /= (double)nPos;

    sumCos2DeltaPhiNeg /= (double)nNeg;
    sumCos2DeltaPhiNeg_sf /= (double)nNeg;
    sumCos2DeltaPhiVertNeg /= (double)nNeg;
    sumCos2DeltaPhiVertNeg_sf /= (double)nNeg;

    sumCos3DeltaPhiPos /= (double)nPos;
    sumCos3DeltaPhiPos_sf /= (double)nPos;
    sumCos3DeltaPhiVertPos /= (double)nPos;
    sumCos3DeltaPhiVertPos_sf /= (double)nPos;

    sumCos3DeltaPhiNeg /= (double)nNeg;
    sumCos3DeltaPhiNeg_sf /= (double)nNeg;
    sumCos3DeltaPhiVertNeg /= (double)nNeg;
    sumCos3DeltaPhiVertNeg_sf /= (double)nNeg;

    //Sin
    double deltaSinS2 = sumSin2DeltaPhiPos - sumSin2DeltaPhiNeg;
    double deltaSinS2_sf = sumSin2DeltaPhiPos_sf - sumSin2DeltaPhiNeg_sf;
    double deltaSinS2Vert = sumSin2DeltaPhiVertPos - sumSin2DeltaPhiVertNeg;
    double deltaSinS2Vert_sf = sumSin2DeltaPhiVertPos_sf - sumSin2DeltaPhiVertNeg_sf;

    double deltaSinS3 = sumSin3DeltaPhiPos - sumSin3DeltaPhiNeg;
    double deltaSinS3_sf = sumSin3DeltaPhiPos_sf - sumSin3DeltaPhiNeg_sf;
    double deltaSinS3Vert = sumSin3DeltaPhiVertPos - sumSin3DeltaPhiVertNeg;
    double deltaSinS3Vert_sf = sumSin3DeltaPhiVertPos_sf - sumSin3DeltaPhiVertNeg_sf;

    //Cos
    double deltaCosS2 = sumCos2DeltaPhiPos - sumCos2DeltaPhiNeg;
    double deltaCosS2_sf = sumCos2DeltaPhiPos_sf - sumCos2DeltaPhiNeg_sf;
    double deltaCosS2Vert = sumCos2DeltaPhiVertPos - sumCos2DeltaPhiVertNeg;
    double deltaCosS2Vert_sf = sumCos2DeltaPhiVertPos_sf - sumCos2DeltaPhiVertNeg_sf;

    double deltaCosS3 = sumCos3DeltaPhiPos - sumCos3DeltaPhiNeg;
    double deltaCosS3_sf = sumCos3DeltaPhiPos_sf - sumCos3DeltaPhiNeg_sf;
    double deltaCosS3Vert = sumCos3DeltaPhiVertPos - sumCos3DeltaPhiVertNeg;
    double deltaCosS3Vert_sf = sumCos3DeltaPhiVertPos_sf - sumCos3DeltaPhiVertNeg_sf;

    //N(deltaSinS)
    hNDeltaSinSPsi2[centBin]->Fill(deltaSinS2);
    //N(deltaSinS_sf)
    hNDeltaSinSPsi2_sf[centBin]->Fill(deltaSinS2_sf);
    //N(deltaSinSVert)
    hNDeltaSinSVertPsi2[centBin]->Fill(deltaSinS2Vert);
    //N(deltaSinSVert_sf)
    hNDeltaSinSVertPsi2_sf[centBin]->Fill(deltaSinS2Vert_sf);
    //N(deltaSinS)
    hNDeltaSinSPsi3[centBin]->Fill(deltaSinS3);
    //N(deltaSinS_sf)
    hNDeltaSinSPsi3_sf[centBin]->Fill(deltaSinS3_sf);
    //N(deltaSinSVert)
    hNDeltaSinSVertPsi3[centBin]->Fill(deltaSinS3Vert);
    //N(deltaSinSVert_sf)
    hNDeltaSinSVertPsi3_sf[centBin]->Fill(deltaSinS3Vert_sf);

    //N(deltaSinS)
    hNDeltaCosSPsi2[centBin]->Fill(deltaCosS2);
    //N(deltaSinS_sf)
    hNDeltaCosSPsi2_sf[centBin]->Fill(deltaCosS2_sf);
    //N(deltaSinSVert)
    hNDeltaCosSVertPsi2[centBin]->Fill(deltaCosS2Vert);
    //N(deltaSinSVert_sf)
    hNDeltaCosSVertPsi2_sf[centBin]->Fill(deltaCosS2Vert_sf);
    //N(deltaSinS)
    hNDeltaCosSPsi3[centBin]->Fill(deltaCosS3);
    //N(deltaSinS_sf)
    hNDeltaCosSPsi3_sf[centBin]->Fill(deltaCosS3_sf);
    //N(deltaSinSVert)
    hNDeltaCosSVertPsi3[centBin]->Fill(deltaCosS3Vert);
    //N(deltaSinSVert_sf)
    hNDeltaCosSVertPsi3_sf[centBin]->Fill(deltaCosS3Vert_sf);

  hEvtCount->Fill(20);
  PostData(1,mOutputList);
}

//---------------------------------------------------
int AliAnalysisTaskRdoubleRatio::GetRunNumBin(int runNum)
{
  int runNumBin=-1;
  // 10h
  int runNumList[91]={
    139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310,
    139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870,
    138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582,
    138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201,
    138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693,
    137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541,
    137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243,
    137236, 137235, 137232, 137231, 137230, 137162, 137161
  };
  // 11h
  // int runNumList[39]={170387,170040,170268,170228,170207,169838,170159,170204,170311,170084,
  //                     169835,170088,170593,170203,170270,169846,170163,170388,170155,170083,
  //                     170572,169837,169855,170306,170269,170089,170309,170091,170081,170230,
  //                     170085,170315,170027,170193,170312,170313,170308,169858,169859};
  for (int i = 0; i < 91; ++i) {
    if (runNum==runNumList[i]) {runNumBin=i; break;}
    else continue;
  }
  return runNumBin;
}