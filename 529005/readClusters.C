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

vector<string> dig2Clus(const std::string &fileName, vector<Cluster>& clusters, vector<Trigger> clusterTriggers);
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

void readClusters(int nEvents) {
  changeFont();
  auto folderName = (gSystem->GetWorkingDirectory());

  auto gwd = gSystem->GetWorkingDirectory();
  std::string path = static_cast<string> (gwd);

  auto runNumber = (gwd.substr(gwd.length() - 9, gwd.length()));



  std::array<std::unique_ptr<TH1F>, 7> hCharge, hMipCharge, hSize;
  std::array<std::unique_ptr<TH2F>, 7> hMap, hMipMap;

  TH1 *hMipCharge2[7];

  // changeFont();			     // specify folder to save files in
  // SaveFolder("clusterChambers");   // apply custom canvas figure options


  vector<Cluster> clusters; vector<Trigger> clusterTriggers;
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
        fileInfo = dig2Clus(pathName, clusters, clusterTriggers);
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

    const char* canStringhMap = Form("Map %i", i);
    hMap[i].reset(new TH2F(canStringhMap, canStringhMap, 160, 0, 159, 144, 0, 143));
    hMap[i]->SetXTitle("X (cm)");
    hMap[i]->SetYTitle("Y (cm)");


    const char* canStringCharge = Form("Charge %i", i);
    hCharge[i].reset(new TH1F(canStringCharge, canStringCharge, 2000, 100., 2100.));
    hCharge[i]->SetXTitle("Charge (ADC channel)");
    hCharge[i]->SetYTitle("Entries");

    const char* canStringMip = Form("MIP-Charge %i", i);
    hMipCharge[i].reset(new TH1F(canStringMip, canStringMip, 50, 200., 2200.));
    hMipCharge[i]->SetXTitle("Charge (ADC channel)");
    hMipCharge[i]->SetYTitle("Entries/40 ADC");
    //hMipCharge[i]->SetLineColor(kBlack);
    hMipCharge[i]->SetStats(kTRUE);
    //hMipCharge[i]->StatOverflows(kTRUE);


    const char* canStringMip2 = Form("MIP-Charge2 %i", i);
    hMipCharge2[i]= new TH1F(canStringMip2, canStringMip2, 50, 200., 2200.);
    hMipCharge2[i]->SetXTitle("Charge (ADC channel)");
    hMipCharge2[i]->SetYTitle("Entries/40 ADC");
    hMipCharge2[i]->SetStats(kTRUE);


    const char* canStringMipMap = Form("MIP-Charge Map %i", i);
    hMipMap[i].reset(new TH2F(canStringMipMap, canStringMipMap, 160, 0, 159, 144, 0, 143));
    hMipMap[i]->SetXTitle("X (cm)");
    hMipMap[i]->SetYTitle("Y (cm)");

    const char* canStringSize = Form("Size %i", i);
    hSize[i].reset(new TH1F(canStringSize, canStringSize, 20, 0., 20.));
    hSize[i]->SetXTitle("Cluster size");
    hSize[i]->SetYTitle("Entries");
  }
 

  std::array<std::unique_ptr<TCanvas>, 5> canvas;  
  
  (canvas[0]).reset(new TCanvas(Form("Cluster-Map %i",fname), Form("Cluster-Map %i",fname), 1200, 1200));
  (canvas[1]).reset(new TCanvas(Form("Cluster-Charge %i",fname), Form("Cluster-Charge %i",fname), 1200, 1200));
  (canvas[2]).reset(new TCanvas(Form("MIP Cluster-Charge %i",fname), Form("MIP Cluster-Charge %i",fname),1200, 1200));  
  (canvas[3]).reset(new TCanvas(Form("Cluster-Size %i",fname), Form("Cluster-Size %i",fname), 1200, 1200));
  (canvas[4]).reset(new TCanvas(Form("MIP Cluster-Map %i",fname), Form("MIP Cluster-Map %i",fname), 1200, 1200));

  TCanvas *mipCanvasRaw = new TCanvas(Form("MipLandau %i",fname), Form("MipLandau %i",fname), 1200, 1200);



  Int_t nTotTriggers = 0;


  //for (int k = 0; k < nEvents; k++) {

    std::unique_ptr<Trigger> pTgr;
    std::unique_ptr<Cluster> pClu;
    std::unique_ptr<Cluster> pCluEvt;
    std::unique_ptr<Cluster> cluster;

    //vector<Cluster> *pClusters = &clusters;
    vector<Cluster> oneEventClusters;
    //vector<Trigger> *pTrigger = &clusterTriggers;



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


      hCharge[module]->Fill(charge);
      if (clusSize >= 3 && clusSize <= 7){
        hMipCharge[module]->Fill(charge);
        hMipCharge2[module]->Fill(charge);
        hMipMap[module]->Fill(x, y);
      }
      hMap[module]->Fill(x, y);

      hSize[module]->Fill(clusSize);
    }

    nTotTriggers += clusterTriggers.size();
    const int triggerSize = clusterTriggers.size();
    Printf("trigger size from clusters = %i", triggerSize);
  
  //auto folderName = fname.c_str();
  std::array<std::unique_ptr<TPaveText>, 5> tpvs;
  

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

  tpvs[0]->AddText("Cluster-Map");
  tpvs[1]->AddText("Clusters Charge");
  tpvs[2]->AddText("MIP Clusters Charge");
  tpvs[3]->AddText("Clusters Size");
  tpvs[4]->AddText("MIP Cluster-Map");

  for(auto& tpv: tpvs){
    tpv->AddText(f1);
    tpv->AddText(f2);
    tpv->AddText(f3);
  }

  //efRaw
  TPaveText *tpvMipRaw = new TPaveText(0.05, .05, .9, .9);
  tpvMipRaw->AddText(runLabel);
  tpvMipRaw->AddText("MIP Landau Raw");
  tpvMipRaw->AddText(f1); tpvMipRaw->AddText(f2); tpvMipRaw->AddText(f3);
  mipCanvasRaw->Divide(3, 3); mipCanvasRaw->cd(3);
  tpvMipRaw->Draw();//*/

  const int tpvSize = static_cast<int>(tpvs.size());
  for(int i = 0; i < tpvSize; i++){
    canvas[i]->Divide(3, 3); 
    canvas[i]->cd(3);
    tpvs[i]->Draw();
  }
  


  std::array<unique_ptr<TF1>, 7> tf1Fit;
  TF1::InitStandardFunctions();
  


  const int posArr[] = {9, 8, 6, 5, 4, 2, 1};
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    
    auto pad0 = static_cast<TPad*>(canvas[0]->cd(pos));


    hMap[iCh]->SetMarkerStyle(3);
    hMap[iCh]->Draw();
    gStyle->SetOptStat("mi");
    gStyle->SetStatX(0.975);
    gStyle->SetStatY(0.975);
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.1);
    
    auto pad1 = static_cast<TPad*>(canvas[1]->cd(pos));
    hCharge[iCh]->Draw();
    
    changeFont();


    auto pad2 = static_cast<TPad*>(canvas[2]->cd(pos));
    hMipCharge[iCh]->SetStats(kTRUE);
    gStyle->SetOptStat("rmi");
    //gStyle->SetOptStat();
    unique_ptr<TF1> fitFunc;
    fitFunc.reset(new TF1(Form("FitFunc%i", iCh), "landau", 0., 4000));
    hMipCharge[iCh]->Fit(fitFunc.get(), Form("FitTH1F%i", iCh), "", 20., 2220.);

    cout << "Values : " << endl;
    const auto& Constant = (fitFunc.get())->GetParameter(0);
    const auto& MPV = (fitFunc.get())->GetParameter(1);
    const auto& Sigma = (fitFunc.get())->GetParameter(2);
    cout << Constant << endl;   
    cout << MPV << endl;   
    cout << Sigma << endl;
  
    hMipCharge[iCh]->SetTitle(Form("Constant %03.1f \n MPV %03.1f Sigma %03.1f", Constant, MPV, Sigma));
    hMipCharge[iCh]->Draw();
    /*std::unique_ptr<TPaveText> tPaveStat;
    tPaveStat.reset(new TPaveText(0.05, .05, .9, .9));
    tPaveStat->AddText(Form("Constant %f", 1.0));
    
    pad2->Update(); 
    tPaveStat->Draw();
    pad2->Modified();*/
    //const auto& fitString = Form("Constant %f MPV %f Sigma %f", Constant, MPV, Sigma);


    //pad2->Add();

    /*
       
     

    TPaveStats *fitStats = static_cast<TPaveStats*>(pad2->GetPrimitive("stats"));
    fitStats->SetName("temp");
    TList *listOfStats = fitStats->GetListOfLines(); 

    TLatex* constStat = new TLatex(0, 0, Form("MPV = %f ", MPV));
    constStat->SetTextFont(42);
    constStat->SetTextSize(0.04);

    // BRYTER HER:::
    //listOfStats->Add(constStat);*/ 

    //unique_ptr<TLatex> constStat; 
    //constStat.reset(new TLatex(0, 0, Form("MPV = %f ", MPV)));



    //hMipCharge[iCh]->SetStats(0);  
    //pad2->Modified();

    //auto a = (hMipCharge[iCh])->GetFunction(Form("FitFunc%i", iCh));

    //cout << *(a->GetFormula()) << endl; 

    auto pad3 = static_cast<TPad*>(canvas[3]->cd(pos));
    hCharge[iCh]->Draw();

    auto pad4 = static_cast<TPad*>(canvas[4]->cd(pos));
    hMipMap[iCh]->Draw();
   
    //*efRaw
    auto pad5 = static_cast<TPad*>(mipCanvasRaw->cd(pos));
    hMipCharge2[iCh]->Fit("landau");
    hMipCharge2[iCh]->Draw();//*/


  }
  
  canvas[0]->SaveAs(Form("clusterMap_%i_.eps",fname));
  canvas[1]->SaveAs(Form("clusterCharge_%i_.eps",fname));
  canvas[2]->SaveAs(Form("mipClustesCharge_%i_.eps",fname));
  canvas[3]->SaveAs(Form("ClusterSize%i_.eps",fname));
  canvas[4]->SaveAs(Form("MipMap%i_.eps",fname));
  mipCanvasRaw->SaveAs(Form("MipLandau%i_.eps",fname));

  canvas[0]->SaveAs(Form("clusterMap_%i_.png",fname));
  canvas[1]->SaveAs(Form("clusterCharge_%i_.png",fname));
  canvas[2]->SaveAs(Form("mipClustesCharge_%i_.png",fname));
  canvas[3]->SaveAs(Form("mipClustesCharge_2%i_.png",fname));
  canvas[4]->SaveAs(Form("MipMap%i_.png",fname));


  sleep_for(5000ms);

  bool userInput = false;
  while(!userInput){
    sleep_for(5000ms);
    string uInputString;
    sleep_for(50ms);
    cin >> uInputString;
    if(uInputString == "C"){
      userInput = true;
    }
    sleep_for(1000ms);
    if(userInput){
      cout << "Got End fro User.. Exiting!";
      mipCanvasRaw->SaveAs(Form("MipLandau%i_.png",fname));
      canvas[0]->SaveAs(Form("clusterMap_%i_.eps",fname));
      canvas[1]->SaveAs(Form("clusterCharge_%i_.eps",fname));
      canvas[2]->SaveAs(Form("mipClustesCharge_%i_.eps",fname));
      canvas[3]->SaveAs(Form("ClusterSize%i_.eps",fname));
      canvas[4]->SaveAs(Form("MipMap%i_.eps",fname));
      mipCanvasRaw->SaveAs(Form("MipLandau%i_.eps",fname));

      canvas[0]->SaveAs(Form("clusterMap_%i_.png",fname));
      canvas[1]->SaveAs(Form("clusterCharge_%i_.png",fname));
      canvas[2]->SaveAs(Form("mipClustesCharge_%i_.png",fname));
      canvas[3]->SaveAs(Form("mipClustesCharge_2%i_.png",fname));
      canvas[4]->SaveAs(Form("MipMap%i_.png",fname));
    }
  }
}

//********************************************************************************************************************

void readDigits() 
{

  const auto& filename = "hmpDig70.root";
  int nEvent = 5;
  std::unique_ptr<TFile> fileDigits;
  fileDigits.reset(static_cast<TFile*>(TFile::Open(filename)));
  // Cast to TTree*, get tree by key "o2sim"

  if(!fileDigits){
    Printf("Problem fetching file %s !", filename);
    Printf("Could not get Digits-file, will return!");
    return;
  }

  std::unique_ptr<TTree> treeDigits;
  treeDigits.reset((TTree*)fileDigits->Get("o2sim"));
  if(!treeDigits){
    treeDigits.reset((TTree *)fileDigits->Get("o2hmp"));
  }

  if(!treeDigits){
    Printf("Problem fetching Tree in file %s !", filename);
    Printf("Could not get Digits-Tree, will return!");
    return;
  }


  std::array<std::unique_ptr<TH1F>, 7> h_xCoord, h_yCoord, hCharge;
  std::array<std::unique_ptr<TH2F>, 7> hMap;

  // Label the histograms
  for (int i = 0; i < 7; i++) {
    // define element number i in the pointer-array

    hMap[i].reset(new TH2F(Form("Digits Map %i", i), Form("Digits Map %i", i), 160,
                       0, 159, 144, 0, 143));
    hMap[i]->SetXTitle("pad X [cm]");
    hMap[i]->SetYTitle("pad Y [cm]");


    hCharge[i].reset(new TH1F(Form("Digits Charge %i", i),
                          Form("Digits Charge %i", i), 2000, 100., 2100.));
    hCharge[i]->SetXTitle("Charge (ADC channel)");
    hCharge[i]->SetYTitle("Entries");

    h_xCoord[i].reset(new TH1F(Form("Digits X-location Histogram %i", i),
                      Form("Digits X-location Histogram %i", i), 2000, 10., 159.));
    h_xCoord[i]->SetXTitle("X [cm]");
    h_xCoord[i]->SetYTitle("Entries");

    h_yCoord[i].reset(new TH1F(Form("Digits Y-location Histogram %i", i),
                 Form("Digits Y-location Histogram %i", i), 2000, 10., 144.));
    h_yCoord[i]->SetXTitle("Y [cm]");
    h_yCoord[i]->SetYTitle("Entries");
  }

  // 		         // apply custom canvas figure options
  // SaveFolder("digitChambers"); // specify folder to save files in

  // Define canvases for plotting the figures
  std::unique_ptr<TCanvas> c1, c2, c3, c4;

  c1.reset( new TCanvas("c1", "c1", 2000, 1200)); c1->Divide(3, 3);
  c2.reset( new TCanvas("c2", "c2", 2000, 1200)); c2->Divide(3, 3);
  c3.reset( new TCanvas("c3", "c3", 2000, 1200)); c3->Divide(3, 3);
  c4.reset( new TCanvas("c4", "c4", 2000, 1200)); c4->Divide(3, 3);


  // Define positions for the plots for the chambers in the canvases
  const int posArr[] = {9, 8, 6, 5, 4, 2, 1};

  std::unique_ptr<Trigger> pTgr;  // pointer to Trigger-object
  std::unique_ptr<Digit> pDig;    // pointer to Digit-object
  std::unique_ptr<Digit> pDigEvt; // pointer to Digit-object
  Digit digit;    // declaration of digit-object


  vector<Digit> *pDigits = nullptr;
  vector<Digit> oneEventDigits;
  vector<Trigger> *pTriggers = nullptr;
  vector<Digit> digits = *pDigits;
  vector<Trigger> triggers = *pTriggers;

  if ((mTree->GetBranchStatus("HMPDigit")) == 1) {
    mTree->SetBranchAddress("HMPDigit", &pDigits);
  } else if ((mTree->GetBranchStatus("HMPIDDigits")) == 1) {
    mTree->SetBranchAddress("HMPIDDigits", &pDigits);
  } else {
    throw std::runtime_error(
        "HMPID DigitToClusterSpec::init() : Error in branches!");
  }
  treeDigits->SetBranchAddress("InteractionRecords", &pTriggers);


  // Number of entries in tree
  const int treeDigSize = treeDigits->GetEntries();
  Printf("tree entries = %i", treeDigSize);

  treeDigits->GetEntry(0);


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
  
  const int triggerSize = static_cast<int>(pTriggers->size());
  Printf("trigger size from digits = %i", triggerSize);

  int i = 0;
  for (const auto& trg : triggers) {
    oneEventDigits.clear(); // empty vector

    for (int j = trg.getFirstEntry(); j <= trg.getLastEntry(); j++) {
      digit = static_cast<Digit>(digits[j]);
      oneEventDigits.push_back(digit);
    }

    if (i == nEvent) {
      for (const auto& digEvt : oneEventDigits) {

        Digit::pad2Absolute(digEvt.getPadID(), &module, &padChX,
                                       &padChY);
        hMap[module]->Fill(padChX, padChY, digEvt.getQ());
      }
    } // end if

    i++;
  }

  int pos;
  for (int iCh = 0; iCh < 7; iCh++) {
    pos = posArr[iCh];
    c2->cd(pos);
    hCharge[iCh]->SetStats(false);

    hCharge[iCh]->Draw();

    c1->cd(pos);
    hMap[iCh]->SetStats(false);
    hMap[iCh]->Draw("Colz");

    c3->cd(pos);
    h_xCoord[iCh]->SetStats(false);
    h_xCoord[iCh]->Draw();

    c4->cd(pos);
    h_yCoord[iCh]->SetStats(false);
    h_yCoord[iCh]->Draw();
  }

  c1->SaveAs("digitmap.eps");
  c2->SaveAs("digitcharge.eps");
  c3->SaveAs("digitX.eps");
  c4->SaveAs("digitY.eps");
}

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


vector<string> dig2Clus(const std::string &fileName, vector<Cluster>& clusters, vector<Trigger> clusterTriggers) {
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

  gStyle->SetTitleSize(.0505, "xzy");
  gStyle->SetTitleOffset(.925, "xz");//.95
  gStyle->SetTitleOffset(1, "y");//1.1
  gStyle->SetTitleFontSize(.05);
  //gStyle->SetTitleFont(16, "xz");
  
  gStyle->SetLabelOffset(0.00625, "y");
  gStyle->SetLabelFont(22, "xyz");
  gStyle->SetLabelSize(.05, "xyz"); //.0525 // verdi av akser


  //mStyle->SetStyle("canvasStyle");
}




