#if !defined(__CLING__) || defined(__ROOTCLING__)
//#if !defined(__CINT__) || defined(__MAKECINT__)
#include "DataFormatsHMP/Cluster.h"
#include "HMPIDReconstruction/Clusterer.h"
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TApplication.h>
#include <TF1.h>
#include <TH2F.h>
#include <TLine.h>
#include <TList.h>
#include <TROOT.h> // gRoot
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>

#include <TRandom.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <vector>
#include <fairlogger/Logger.h>

// C++ header files and libraries
#include <chrono>
#include <thread>
#include <ctime>
#include <fstream>
#include <iostream>
#include <gsl/gsl>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#endif


using std::this_thread::sleep_for;
using std::vector, std::cout, std::cin, std::endl;
using o2::hmpid::Cluster, o2::hmpid::Digit, o2::hmpid::Trigger, o2::hmpid::Clusterer;

std::vector<o2::hmpid::Digit> mDigitsFromFile,
*mDigitsFromFilePtr = &mDigitsFromFile;
std::vector<o2::hmpid::Trigger> mTriggersFromFile,
*mTriggersFromFilePtr = &mTriggersFromFile;
// void SaveFolder(string inpFile);
void changeFont();

vector<string> dig2Clus(const std::string &fileName, vector<Cluster>& clusters, vector<Trigger>& clusterTriggers, vector<Digit>& digits);
bool mReadFile = false;
std::string mSigmaCutPar;
float mSigmaCut[7] = {4, 4, 4, 4, 4, 4, 4};

std::unique_ptr<TFile> mFile; ///< input file containin the tree
std::unique_ptr<TTree> mTree; ///< input tree



std::unique_ptr<Clusterer> mRec; // ef: changed to smart-pointer


void initFileIn(const std::string &fileName);

void strToFloatsSplit(std::string s, std::string delimiter, float *res,
                      int maxElem = 7);

#include <filesystem>
#include <iostream>
#include <string>
namespace fs = std::filesystem;
void readDigits();
void fillDigMap(vector<Digit>& digits);
void readClusters(int nEvents) {
  changeFont();
  auto folderName = (gSystem->GetWorkingDirectory());

  auto gwd = gSystem->GetWorkingDirectory();
  std::string path = static_cast<string> (gwd);

  auto runNumber = (gwd.substr(gwd.length() - 9, gwd.length()));

  std::array<std::unique_ptr<TH1F>, 7> digCharge, hMipCharge;
  std::array<std::unique_ptr<TH2F>, 7> digMap;

  // changeFont();			     // specify folder to save files in
  // SaveFolder("clusterChambers");   // apply custom canvas figure options


  vector<Cluster> clusters; vector<Trigger> clusterTriggers;
  vector<Digit> digits;
  int fname = 000000;
  vector<string> fileInfo;


  bool fileFound = false;
  for (const auto &entry : fs::directory_iterator(path)) {
    const auto& pathName = static_cast<string>(entry.path());
    //std::cout << pathName << std::endl;
    if (pathName.size() < 5) continue;
    const auto &fileEnd =
        (pathName.substr(pathName.length() - 5, pathName.length()));
    if (fileEnd != ".root") continue;
    std::cout << " dig2clus " << pathName << std::endl;
    

    auto tmp = pathName.substr(0,  pathName.find_last_of("\\/"));

    auto folderName = tmp.substr(tmp.length()-6);
    fname = stoi(folderName);
    if(folderName.length() == 6){
      if(std::all_of(folderName.begin(), folderName.end(), ::isdigit)){
        fileFound = true;
        std::cout << " folderName " << folderName << std::endl;
        fileInfo = dig2Clus(pathName, clusters, clusterTriggers, digits);
	char* fn = strdup(folderName.c_str());
	std::cout << " fname " << fn << std::endl;
        //readDigits();
      }
    }
    //cout << " Digits file << with {} Clusters, {} Triggers", , clusters.size(), clusterTriggers.size());	
  }
  if(!fileFound){
    cout << "No fitting file found!";
  }
  

  for (int i = 0; i < 7; i++) {

    const char* canStringMip = Form("MIP-Charge %i", i);
    //hMipCharge[i].reset(new TH1F(canStringMip, canStringMip, 50, 200., 2200.));
    hMipCharge[i].reset(new TH1F(canStringMip, canStringMip, 50, 200., 2200.));
    hMipCharge[i]->SetXTitle("Charge (ADC channel)");
    hMipCharge[i]->SetYTitle("Entries/40 ADC");
    hMipCharge[i]->SetStats(kTRUE);

    const char* canStringSize = Form("Digit Charge %i", i);
    digCharge[i].reset(new TH1F(canStringSize, canStringSize, 2000, 100., 2100.));
    digCharge[i]->SetXTitle("Charge (ADC channel)");
    digCharge[i]->SetYTitle("Entries/40 ADC");



    const char* canDigMap = Form("Digit Map %i", i);
    digMap[i].reset(new TH2F(canDigMap, canDigMap, 160, 0, 159, 144, 0, 143));
    //digMap[i].reset(new TH2F(canDigMap, canDigMap, 160*0.8, 0, 159*0.8, 144*0.8, 0, 160*0.8));
    digMap[i]->SetXTitle("x [cm]");
    digMap[i]->SetYTitle("y [cm]");
  }
 

  std::array<std::unique_ptr<TCanvas>, 3> canvas;  

  (canvas[0]).reset(new TCanvas(Form("MIP Cluster-Charge %i",fname), Form("MIP Cluster-Charge %i",fname),1200, 1200));  



  (canvas[1]).reset(new TCanvas(Form("Digit Charge %i",fname), Form("Digit Charge %i",fname), 1200, 1200));



  (canvas[2]).reset(new TCanvas(Form("Digit-Map %i",fname), Form("Digit-Map %i",fname), 1200, 1200));


  canvas[0]->SetLeftMargin(.1);
  canvas[1]->SetLeftMargin(.15);
  canvas[2]->SetLeftMargin(.125);


  Int_t nTotTriggers = 0;

    std::unique_ptr<Trigger> pTgr;
    std::unique_ptr<Cluster> pClu;
    std::unique_ptr<Cluster> pCluEvt;
    std::unique_ptr<Cluster> cluster;

    //vector<Cluster> *pClusters = &clusters;
    vector<Cluster> oneEventClusters;
    //vector<Trigger> *pTrigger = &clusterTriggers;

    const int digSize = digits.size();
    Printf("digit size = %i",digSize);
    int padChX = 0, padChY = 0, module = 0;



    for(const auto& dig : digits){
      Digit::pad2Absolute(dig.getPadID(), &module, &padChX, &padChY);
      digCharge[module]->Fill(dig.getQ());
      digMap[module]->Fill(padChX, padChY, dig.getQ());// ef removed Charge , dig.getQ());      // fill y-coordinate
    }

    const int clusterSize = clusters.size();
    Printf("clusters size = %i", clusterSize);


    float minCharge = 9999.0f; float maxCharge = 0.0f;


    for (auto& clus : clusters) {
      const auto& charge = clus.q();
      const auto& clusSize = clus.size();
      const auto& x = clus.x();
      const auto& y = clus.y();
      const auto& module = clus.ch();   

      if(minCharge < charge) {minCharge = charge;}
      if(minCharge < charge) {minCharge = charge;}


      if (clusSize >= 3 && clusSize <= 7){
        hMipCharge[module]->Fill(charge);
      }
    }

    nTotTriggers += clusterTriggers.size();
    const int triggerSize = clusterTriggers.size();
    Printf("trigger size from clusters = %i", triggerSize);
  
  //auto folderName = fname.c_str();
  std::array<std::unique_ptr<TPaveText>, 3> tpvs;
  

  //fileInfo
  const auto f1 = (fileInfo[0]).c_str();
  const auto f2 = (fileInfo[1]).c_str();
  const auto f3 = (fileInfo[2]).c_str();

  //const char* runLabel = Form("%i  Duration = %s", fname, f1);
  const char* runLabel = Form("%i", fname);
  for(auto& tpv: tpvs){
    tpv.reset(new TPaveText(0.05, .05, .9, .9));
    tpv->AddText(runLabel);
  }

  tpvs[0]->AddText("MIP Clusters Charge");
  tpvs[1]->AddText("Digits-Charge");
  tpvs[2]->AddText("Digits-Map");


  for(auto& tpv: tpvs){
    tpv->AddText(f1);
    tpv->AddText(f2);
    tpv->AddText(f3);
  }

  const int tpvSize = static_cast<int>(tpvs.size());
  for(int i = 0; i < tpvSize; i++){
    canvas[i]->Divide(3, 3); 
    canvas[i]->cd(3);
    tpvs[i]->Draw();
  }
  changeFont();    

  gStyle->SetStatX(0.975);
  gStyle->SetStatY(0.975);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3); 
    

  const int posArr[] = {9, 8, 6, 5, 4, 2, 1};
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== Digit MAP =========================
    auto pad5 = static_cast<TPad*>(canvas[2]->cd(pos));
    const auto& pTotalDigs = static_cast<float>(100.0f*digMap[iCh]->GetEntries()/digSize);
    digMap[iCh]->SetTitle(Form("Chamber %i  of total = %02.1f", iCh, pTotalDigs));
    digMap[iCh]->Draw();
  }
  gStyle->SetOptStat("e");





  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== MIP Charge =========================
    auto pad0 = static_cast<TPad*>(canvas[0]->cd(pos)); // Constant*TMath::Landau(1, [MPV], sigma, 0)
    hMipCharge[iCh]->Fit("landau", "I"); // I = fit by integral
    //hMipCharge[iCh]->SetTitle(Form("Constant %03.1f \n MPV %03.1f Sigma %03.1f", Constant, MPV, Sigma));
    hMipCharge[iCh]->Draw();
    hMipCharge[iCh]->FitPanel();
  }

  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3); 
  gStyle->SetOptStat("emr");
  gStyle->SetLabelOffset(0.00525, "y");

  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];

    // ========== Digit Charge =========================
    auto pad3 = static_cast<TPad*>(canvas[1]->cd(pos));
    digCharge[iCh]->Draw();
  }


  gStyle->SetOptStat("eim");
  gStyle->SetLabelOffset(0.008, "y");

  canvas[0]->SetLeftMargin(.1);
  canvas[1]->SetLeftMargin(.15);
  canvas[2]->SetLeftMargin(.125);

  canvas[0]->SaveAs(Form("MIP_Cluster_Charge_%i_.png",fname));
  canvas[1]->SaveAs(Form("Digit_Charge_%i_.png",fname));
  canvas[2]->SaveAs(Form("Digit_Map_%i_.png",fname));

  canvas[0]->Show();
  canvas[1]->Show();
  canvas[2]->Show();

  sleep_for(5000ms);

  bool userInput = false;
  while(!userInput){
    sleep_for(5000ms);
    string uInputString;
    sleep_for(50ms);
    cin >> uInputString;
    if(uInputString == "C"){

      canvas[0]->Show();
      canvas[1]->Show();
      canvas[2]->Show();
      sleep_for(10000ms);

      while(uInputString != "Q"){
        cin >> uInputString;
        sleep_for(10000ms);
      }
    }
    sleep_for(1000ms);
    if(userInput){
      canvas[0]->Close();
      canvas[1]->Close();
      canvas[2]->Close();
      cout << "Got End fro User.. Exiting!";
      canvas[0]->SaveAs(Form("MIP_Cluster_Charge_%i_.png",fname));
      canvas[1]->SaveAs(Form("Digit_Charge_%i_.png",fname));
      canvas[2]->SaveAs(Form("Digit_Map_%i_.png",fname));
    }
  }
}

/*
void fillDigMap(vector<Digit>& digits);
{
  const int digSize = pDigits->size();
  Printf("digit size = %i",digSize);
  int padChX = 0, padChY = 0, module = 0;

  // Loop through digits in file

  for(const auto& dig : digits){
    Digit::pad2Absolute(dig.getPadID(), &module, &padChX, &padChY);
    Digit::pad2Absolute(dig.getPadID(), &module, &padChX, &padChY);
    hCharge[module]->Fill(dig.getQ()); // getQ gets charge
    h_xCoord[module]->Fill(padChX);      // fill x-coordinate
    h_yCoord[module]->Fill(padChY);      // fill y-coordinate
  }
}*/



void strToFloatsSplit(std::string s, std::string delimiter, float *res,
                      int maxElem) {
  int index = 0;
  size_t pos_start = 0;
  size_t pos_end;
  size_t delim_len = delimiter.length();
  std::string token;
  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res[index++] = std::stof(token);
    if (index == maxElem) {
      return;
    }
  }
  res[index++] = (std::stof(s.substr(pos_start)));
  return;
}


vector<string> dig2Clus(const std::string &fileName, vector<Cluster>& clusters, vector<Trigger>& clusterTriggers, vector<Digit>& digits) {
  long mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
  uint32_t firstTrigger, lastTrigger = 0;
  
  mRec.reset(new o2::hmpid::Clusterer()); // ef: changed to smart-pointer

  // mExTimer.start();
  initFileIn(fileName);

  std::cout << "[HMPID DClusterization - run() ] Enter ...";

  
  clusters.clear();
  clusterTriggers.clear();
  cout << "[HMPID DClusterization - run() ] Entries  = " << mTree->GetEntries() << endl; 

  // check if more entries in tree
  if (mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {
    //mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
    //firstTrigger, lastTrigger = 0;
  } else {
    auto entry = mTree->GetReadEntry() + 1;
    assert(entry < mTree->GetEntries());
    mTree->GetEntry(entry);
        
    firstTrigger = (mTriggersFromFilePtr->at(0)).getOrbit();
    lastTrigger = (mTriggersFromFilePtr->back()).getOrbit();

    for (const auto &trig : *mTriggersFromFilePtr) {
      if (trig.getNumberOfObjects()) {
        gsl::span<const Digit> trigDigits{
            mDigitsFromFilePtr->data() + trig.getFirstEntry(),
            size_t(trig.getNumberOfObjects())};
        const size_t clStart = clusters.size();
        mRec->Dig2Clu(trigDigits, clusters, mSigmaCut, true);
        clusterTriggers.emplace_back(trig.getIr(), clStart,
                                     clusters.size() - clStart);
      }
    }

    cout << "Received" << mTriggersFromFilePtr->size() << "triggers with " << mDigitsFromFilePtr->size() << "digits -> clusters = " << clusters.size();

    digits = *mDigitsFromFilePtr;
    if(digits.size() == 0){
      digits = mDigitsFromFile;
    }
    mDigitsReceived = mDigitsFromFilePtr->size();
    mClustersReceived = clusters.size();
    mTriggersReceived = mTriggersFromFilePtr->size();
  }


  const int numTriggers = static_cast<int>(mTriggersReceived);
  const int numDigits = static_cast<int>(mDigitsReceived);
  const int numClusters = static_cast<int>(mClustersReceived);
  

  const auto& trigInfo = Form("Triggers %i" , numTriggers); 
  const auto& digClusInfo = Form("Digits %i Clusters %i",
              numDigits, numClusters);
  //const auto trigger = Form("First Entry %i Last %i", firstTrigger, lastTrigger);
  
  int durSec = static_cast<int>((lastTrigger - firstTrigger)/1000);
  auto durMin = durSec%60;
  
  const auto durInfo = Form("Duration of triggers = %i min", durMin);
  return  {trigInfo, digClusInfo, durInfo};
}
void initFileIn(const std::string &filename) {
 
  long mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
  // Create the TFIle
  mFile = std::make_unique<TFile>(filename.c_str(), "OLD");
  assert(mFile && !mFile->IsZombie());

  mTree.reset((TTree *)mFile->Get("o2sim"));
  if (!mTree) {
    mTree.reset((TTree *)mFile->Get("o2hmp"));
  }

  if (!mTree) {
    LOG(error)
        << "HMPID DigitToClusterSpec::init() : Did not find o2sim tree in "
        << filename.c_str();
    throw std::runtime_error("HMPID DigitToClusterSpec::init() : Did not find "
                             "o2sim file in digits tree");
  }

  if ((mTree->GetBranchStatus("HMPDigit")) == 1) {
    mTree->SetBranchAddress("HMPDigit", &mDigitsFromFilePtr);
  } else if ((mTree->GetBranchStatus("HMPIDDigits")) == 1) {
    mTree->SetBranchAddress("HMPIDDigits", &mDigitsFromFilePtr);
  } else {
    throw std::runtime_error(
        "HMPID DigitToClusterSpec::init() : Error in branches!");
  }

  mTree->SetBranchAddress("InteractionRecords", &mTriggersFromFilePtr);
  mTree->Print("toponly");
}









void changeFont()
{

  /*
  std::unique_ptr<TStyle> mStyle; 
  mStyle.reset(new TStyle("canvasStyle", "Canvas Root Styles"));
  */ 

  gStyle->SetStatX(0.925);
  gStyle->SetStatY(0.925);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.25);
  gStyle->SetStatFontSize(0.065);
  gStyle->SetLegendTextSize(0.065);//

  gStyle->SetTitleSize(.055, "xzy");
  gStyle->SetTitleOffset(.925, "xz");//.95
  gStyle->SetTitleOffset(1, "y");//1.1
  gStyle->SetTitleFontSize(.05);
  //gStyle->SetTitleFont(16, "xz");
  
  gStyle->SetLabelOffset(0.0065, "y");
  gStyle->SetLabelFont(22, "xyz");
  gStyle->SetLabelSize(.055, "xyz"); //.0525 // verdi av akser


  //mStyle->SetStyle("canvasStyle");
}




