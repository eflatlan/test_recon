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


#include "CommonDataFormat/InteractionRecord.h"

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

using o2::InteractionRecord;

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

bool padDig[7][160][144] = {{{true}}};

void setPadChannel(bool (&padDigOff)[7][160][144], int chamber, int xLow, int xHigh, int yLow, int yHigh);

void readDigits();
void fillDigMap(vector<Digit>& digits);


void readClusters(int nEvents) 
{
  changeFont();
  auto folderName = (gSystem->GetWorkingDirectory());

  auto gwd = gSystem->GetWorkingDirectory();
  std::string path = static_cast<string> (gwd);

  auto runNumber = (gwd.substr(gwd.length() - 9, gwd.length()));

  std::array<std::unique_ptr<TGraph>, 7> trigGraph;
  
  std::array<std::unique_ptr<TH1F>, 7> digCharge, hMipCharge, digPerEvent;
  std::array<std::unique_ptr<TH2F>, 7> digMap, digMapAvg, digMapCount;



  std::unique_ptr<TH1F> triggerTimeFreqHist, digPerEventFreq;
  std::unique_ptr<TGraph> trigTime;
  //double padDigits[7][160][144];
  // changeFont();			     // specify folder to save files in
  // SaveFolder("clusterChambers");   // apply custom canvas figure options


  // do not fill high values for selected values
  bool padDigOff[7][160][144] = {{{true}}};
  for(int chamber = 0; chamber < 7; chamber++){
    for(int x = 0; x < 160; x++){
      for(int y = 0; y < 144; y++){
        padDig[chamber][x][y] = true;
      }
    }
  }
  


  setPadChannel(padDigOff, 6, 141, 150, 105, 110); // chamber, xLow, xHigh, yLow,  yHigh
  
  for(int chamber = 0; chamber < 7; chamber++){
  for(int x = 0; x < 160; x++){
    for(int y = 0; y < 144; y++){
      if(!padDig[chamber][x][y])
      {cout << "False padDigOff " << chamber << " x " << x << " y " << y << endl;}
    }
  }}

  vector<Cluster> clusters; vector<Trigger> clusterTriggers;
  vector<Digit> digits;
  int fname = 000000;
  vector<string> fileInfo;
  int numTriggers; 

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
        numTriggers = clusterTriggers.size();
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

  float avgDigits = static_cast<float>(1.0f*digits.size()/numTriggers);
  

  for (int i = 0; i < 7; i++) {

    const char* canStringMip = Form("MIP-Charge %i", i);
    hMipCharge[i].reset(new TH1F(canStringMip, canStringMip, 50, 200., 2200.));
    //hMipCharge[i].reset(new TH1F(canStringMip, canStringMip, 500, 200., 2200.));
    hMipCharge[i]->SetXTitle("Charge (ADC channel)");
    hMipCharge[i]->SetYTitle("Entries/40 ADC");
    hMipCharge[i]->SetStats(kTRUE);

    const char* canStringSize = Form("Digit Charge %i", i);
    digCharge[i].reset(new TH1F(canStringSize, canStringSize, 50, 0., 400.));
    digCharge[i]->SetXTitle("Charge (ADC channel)");
    digCharge[i]->SetYTitle("Entries/40 ADC");
    digCharge[i]->SetLabelOffset(0.0065, "y");


    const char* canDigMap = Form("Digit Map %i", i);
    digMap[i].reset(new TH2F(canDigMap, canDigMap, 160, 0, 159, 144, 0, 143));
    digMap[i]->SetXTitle("x [cm]");
    digMap[i]->SetYTitle("y [cm]");

    const char* canDigAvg = Form("Avg Charge Per Pad %i", i);
    digMapAvg[i].reset(new TH2F(canDigAvg, canDigAvg, 160, 0, 159, 144, 0, 143));
    //digMap[i].reset(new TH2F(canDigMap, canDigMap, 160*0.8, 0, 159*0.8, 144*0.8, 0, 160*0.8));
    digMapAvg[i]->SetXTitle("x [cm]");
    digMapAvg[i]->SetYTitle("y [cm]");


    //const char* canDigMap = Form("Digit Map %i", i);
    //strigGraph.reset(new TGraph(, , 160, 0, 159, 144, 0, 143));

  const char* digEvtFreqStr = Form("Digits Per Event Frequency%i",i);
  digPerEvent[i].reset(new TH1F(digEvtFreqStr, digEvtFreqStr, 500, 0., 1.));
  //digPerEvent[i].reset(new TH1F(digEvtFreqStr, digEvtFreqStr, 500, 100*(avgDigits-150.)/(144*160), 100*(avgDigits+150.)/(144*160)));
    digPerEvent[i]->SetXTitle("Number of digits");
    digPerEvent[i]->SetYTitle("Frequencies");
   
  }
  

  const char* trigTimeStr = Form("Trigger Time Freq");
  triggerTimeFreqHist.reset(new TH1F(trigTimeStr, trigTimeStr, numTriggers, 0., 500000.));
  triggerTimeFreqHist->SetXTitle("Trigger Time");
  triggerTimeFreqHist->SetYTitle("Frequency");
 
  const char* trigEvtStr = Form("Trigger Time; Event Number; Delta Time");
  trigTime.reset(new TGraph);
  trigTime->SetTitle(trigEvtStr);

  std::array<std::unique_ptr<TCanvas>, 7> canvas;  

  (canvas[0]).reset(new TCanvas(Form("MIP Cluster-Charge %i",fname), Form("MIP Cluster-Charge %i",fname),1200, 1200));  

  (canvas[1]).reset(new TCanvas(Form("Digit Charge %i",fname), Form("Digit Charge %i",fname), 1200, 1200));

  (canvas[2]).reset(new TCanvas(Form("Digit-Map %i",fname), Form("Digit-Map %i",fname), 1200, 1200));

  (canvas[3]).reset(new TCanvas(Form("Digit-Map Avg %i",fname), Form("Digit-Map Avg %i",fname), 1200, 1200));

  (canvas[4]).reset(new TCanvas(Form("Digits Per Event %i",fname), Form("Digits Per Event %i",fname), 1200, 1200));
  


  (canvas[5]).reset(new TCanvas(Form("Trigger Time %i",fname), Form("Trigger Time %i",fname), 1200, 1200));

  (canvas[6]).reset(new TCanvas(Form("Digits Freq %i",fname), Form("Digits Freq  %i",fname), 1200, 1200));


  canvas[0]->SetLeftMargin(.1+canvas[0]->GetLeftMargin());
  canvas[1]->SetLeftMargin(.175+canvas[1]->GetLeftMargin());
  canvas[2]->SetLeftMargin(.1+canvas[2]->GetLeftMargin());
  canvas[3]->SetLeftMargin(.1+canvas[3]->GetLeftMargin());
  canvas[4]->SetLeftMargin(.1+canvas[4]->GetLeftMargin());
  Int_t nTotTriggers = 0;
     
    vector<Digit> oneEventDigits;
    std::unique_ptr<Trigger> pTgr;
    std::unique_ptr<Cluster> pClu, pCluEvt, cluster;
    //std::unique_ptr<Cluster> pCluEvt;
    //std::unique_ptr<Cluster> cluster;

    //vector<Cluster> *pClusters = &clusters;
    vector<Cluster> oneEventClusters;
    //vector<Trigger> *pTrigger = &clusterTriggers;

    const int digSize = digits.size();
    Printf("digit size = %i",digSize);
    int padChX = 0, padChY = 0, module = 0;

    int trigNum = 0;
    double tprev;
    Trigger trigPrev;

    std::array<float, 7> avgDig;

    const auto rel = 100.0/(144*160);
    for(const auto& trig : mTriggersFromFile){
      oneEventDigits.clear();
      const int numDigPerTrig = trig.getNumberOfObjects();
      const int firstTrig = trig.getFirstEntry();
      const int lastTrig = trig.getLastEntry();


      const auto& orbit = trig.getOrbit();
      const auto& bc = trig.getBc();
      const auto& time = InteractionRecord::bc2ns(bc, orbit);
      


      // må endres til per chamber

      const auto& tDif = (trig.getIr()).differenceInBCNS(trigPrev.getIr());

      
      if(trigNum > 0 && trigNum < 50){

        if(tDif < -1000000 || tDif > 1000000){ 
          //trigTime->SetPoint(trigNum, trigNum, 0.0f);
        } else {trigTime->SetPoint(trigNum, trigNum, tDif);}
        //cout << "Filled " << trigNum << " Time " << tDif;
      }

      triggerTimeFreqHist->Fill(tDif);

      //digPerEventFreq->Fill(numDigPerTrig*rel);

      //cout << trigNum << "Filled histograms w triggerTimeFreqHist " << tDif << endl << " digPerFreq " << numDigPerTrig*rel << endl;
      std::array<int, 7> cntCh = {0,0,0,0,0,0,0};

      for(int j = firstTrig; j < lastTrig; j++){
        const auto& dig = digits[j];
        Digit::pad2Absolute(dig.getPadID(), &module, &padChX, &padChY);
        //cout << "Fill Chamber, trNum " <<  module << " " << trigNum << endl;
        //cout << "Total " << cntCh[module] << endl;
	cntCh[module]++;
	avgDig[module]++;
	
      }
      for(auto c : cntCh){cout << "  " << c << endl;}
 	cout << endl;


      for(int ch = 0; ch < 7; ch++){
        digPerEvent[ch]->Fill(rel*cntCh[ch]);
      }


      tprev = time;
      trigPrev = trig;
      trigNum++; 

    }



    for(const auto& dig : digits){
      Digit::pad2Absolute(dig.getPadID(), &module, &padChX, &padChY);

      if(padDig[module][padChX][padChY] == true){
        digCharge[module]->Fill(dig.getQ());
      }

      digMap[module]->Fill(padChX, padChY, dig.getQ());
      digMapAvg[module]->Fill(padChX, padChY, dig.getQ()/5);//numTriggers
      //digMap[module]->Fill(padChX, padChY, padDigOff[module][padChX][padChY]);
      //padDigits[module][padChX][padChY] += dig.getQ();
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
  std::array<std::unique_ptr<TPaveText>, 7> tpvs;
  

  //fileInfo
  const auto f1 = (fileInfo[0]).c_str();
  const auto f2 = (fileInfo[1]).c_str();
  const auto f3 = (fileInfo[2]).c_str();
  const auto f4 = (fileInfo[3]).c_str();

  //const char* runLabel = Form("%i  Duration = %s", fname, f1);
  const char* runLabel = Form("%i", fname);


  vector<const char*> tpvTexts{"MIP Clusters Charge", "Digits-Charge", "Digits-Map", "Digits-Map Avg", "Digits Per Event", "", ""};

  int j = 0;
  for(auto& tpv: tpvs){
    tpv.reset(new TPaveText(0.05, .05, .9, .9));
    tpv->AddText(Form("%s %s", runLabel, tpvTexts[j++]));
  }

  /*
  tpvs[0]->AddText("MIP Clusters Charge");
  tpvs[1]->AddText("Digits-Charge");
  tpvs[2]->AddText("Digits-Map"); */ 


  for(auto& tpv: tpvs){
    tpv->AddText(f1);
    tpv->AddText(f2);
    tpv->AddText(f3);
    tpv->AddText(f4);
  }



  const int tpvSize = static_cast<int>(tpvs.size());
  for(int i = 0; i < tpvSize; i++){
    canvas[i]->Divide(3, 3); 
    canvas[i]->cd(3);
    tpvs[i]->Draw();
  }



  const int posArr[] = {9, 8, 6, 5, 4, 2, 1};
  changeFont();    

  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3); 
    


  std::unique_ptr<TCanvas> temp1;
  temp1.reset(new TCanvas(Form("temp1%i",fname), Form("temp1%i",fname),1200, 1200));

  temp1->Divide(2,1);
  { 

    temp1->cd(2);
    //tpvs[0]->Draw();
    digPerEvent[0]->Draw();

    auto pad2 = static_cast<TPad*>(temp1->cd(1));
    triggerTimeFreqHist->Draw();
  }
  temp1->Show();
  temp1->SaveAs(Form("TriggerFreq_%i_.png",fname));


  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    auto pad5 = static_cast<TPad*>(canvas[4]->cd(pos));
    digPerEvent[iCh]->Draw();
  }


  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== Digit MAP =========================


    auto pad5 = static_cast<TPad*>(canvas[4]->cd(pos));
    digPerEvent[iCh]->Draw();
    //pad5->SetLeftMargin(+.025+pad5->GetLeftMargin());


    /*const auto& pTotalDigs = static_cast<float>(100.0f*digMap[iCh]->GetEntries()/digSize);
    
    digMap[iCh]->SetLabelOffset(digMap[iCh]->GetLabelOffset("y")-0.0015, "y");
    digMap[iCh]->SetTitleOffset(digMap[iCh]->GetTitleOffset("y")-0.0015, "y");
    digMap[iCh]->SetTitleOffset(digMap[iCh]->GetTitleOffset("x")-0.0005, "x");

    pad5->SetBottomMargin(.0015+pad5->GetBottomMargin());
    pad5->SetRightMargin(-.0025+pad5->GetRightMargin());
    digMap[iCh]->SetTitle(Form("Chamber %i Percentage of total = %02.0f", iCh, pTotalDigs));
    digMap[iCh]->SetMarkerStyle(3);
    digMap[iCh]->Draw("Colz"); */
  }


  // avg digits charge
  for (int iCh = 0; iCh < 7; iCh++) {
    //(*digMapAvg[iCh]) = (*digMap[iCh])/(*digMapCount[iCh]);
    const auto& pos = posArr[iCh];
    // ========== Digit MAP =========================
    auto pad5 = static_cast<TPad*>(canvas[3]->cd(pos));
    //pad5->SetLeftMargin(+.025+pad5->GetLeftMargin());
    const auto& pTotalDigs = static_cast<float>(100.0f*digMapAvg[iCh]->GetEntries()/digSize);

    pad5->SetBottomMargin(.0015+pad5->GetBottomMargin());
    pad5->SetRightMargin(-.0025+pad5->GetRightMargin());
    digMapAvg[iCh]->SetTitle(Form("Chamber Avg %i Percentage of total = %02.0f", iCh, pTotalDigs));
    digMapAvg[iCh]->Draw("Colz");

    const auto maxPos = digPerEvent[iCh]->GetMaximumBin();
    //digPerEvent[iCh]->SetBins(500, 0., 0.5);
    digPerEvent[iCh]->Draw();

    cout << "avg "  << avgDig[iCh] << " Num digs " << digPerEvent[iCh]->GetEntries() << endl;
  }



  gStyle->SetOptStat("e");
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6); 
    
  gStyle->SetOptStat("e");
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6); 
    
  /*
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];

    auto pad5 = static_cast<TPad*>(canvas[4]->cd(pos));
    //pad5->SetLeftMargin(+.025+pad5->GetLeftMargin());
    pad5->SetBottomMargin(.0015+pad5->GetBottomMargin());
    pad5->SetRightMargin(-.0025+pad5->GetRightMargin());
    digPerEvent[iCh]->Draw();
  } */
  
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== MIP Charge =========================
    auto pad0 = static_cast<TPad*>(canvas[0]->cd(pos)); // Constant*TMath::Landau(1, [MPV], sigma, 0)

    pad0->SetLeftMargin(-.0025+pad0->GetLeftMargin());
    pad0->SetRightMargin(-.005+pad0->GetRightMargin());
    pad0->SetBottomMargin(.0025+pad0->GetBottomMargin());
    hMipCharge[iCh]->Fit("landau", "I"); // I = fit by integral
    hMipCharge[iCh]->SetLabelOffset(hMipCharge[iCh]->GetLabelOffset("y")-0.0025, "y");
    hMipCharge[iCh]->SetTitleOffset(0.8, "y");
    //hMipCharge[iCh]->SetTitle(Form("Constant %03.1f \n MPV %03.1f Sigma %03.1f", Constant, MPV, Sigma));
    hMipCharge[iCh]->Draw();

  }

  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6); 
  gStyle->SetOptStat("eimr");
  gStyle->SetLabelOffset(0.00525, "y");

  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];

    // ========== Digit Charge =========================
    auto pad3 = static_cast<TPad*>(canvas[1]->cd(pos));
    pad3->SetBottomMargin(.0025+pad3->GetBottomMargin());
    pad3->SetLeftMargin(.05+pad3->GetLeftMargin());
    digCharge[iCh]->SetLabelOffset(digCharge[iCh]->GetLabelOffset("y")+0.0015, "y");
    digCharge[iCh]->SetTitleOffset(1.2,"y");
    pad3->SetRightMargin(-.0025+pad3->GetRightMargin());
    digCharge[iCh]->Draw();
  }

  gStyle->SetStatH(0.2); 

  gStyle->SetOptStat("eim");
  gStyle->SetLabelOffset(0.008, "y");




  canvas[0]->SaveAs(Form("MIP_Cluster_Charge_%i_.png",fname));
  canvas[1]->SaveAs(Form("Digit_Charge_%i_.png",fname));
  canvas[2]->SaveAs(Form("Digit_Map_%i_.png",fname));
  canvas[3]->SaveAs(Form("Digit_Map_Avg_%i_.png",fname));
  canvas[4]->SaveAs(Form("DigitsPerEvent_%i_.png",fname));
  //canvas[5]->SaveAs(Form("TriggerTime_%i_.png",fname));
  //canvas[6]->SaveAs(Form("TriggerTime_%i_.png",fname));

  canvas[0]->Show();
  canvas[1]->Show();
  canvas[2]->Show();
  canvas[3]->Show();
  canvas[4]->Show();
  //canvas[5]->Show();
  //canvas[6]->Show();

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
      canvas[3]->Show();
      canvas[4]->Show();
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
      canvas[3]->Close();
      canvas[4]->Close();
      cout << "Got End fro User.. Exiting!";
      canvas[0]->SaveAs(Form("MIP_Cluster_Charge_%i_.png",fname));
      canvas[1]->SaveAs(Form("Digit_Charge_%i_.png",fname));
      canvas[2]->SaveAs(Form("Digit_Map_%i_.png",fname));
    }
  }
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
  double durMin, durSec = 0.0;
  // check if more entries in tree
  if (mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {
    //mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
    //firstTrigger, lastTrigger = 0;
  } else {
    auto entry = mTree->GetReadEntry() + 1;
    assert(entry < mTree->GetEntries());
    mTree->GetEntry(entry);
        
    Trigger trigFirst = mTriggersFromFilePtr->at(0);
    Trigger trigLast = mTriggersFromFilePtr->back();   

    const auto& tDif = (trigLast.getIr()).differenceInBCNS(trigFirst.getIr());

    durSec = static_cast<double>((tDif)/1000000000.0);
    durMin = static_cast<double>((durSec)/60.0);

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

    cout << " Received " << mTriggersFromFilePtr->size() << " triggers with " << mDigitsFromFilePtr->size() << " digits -> clusters = " << clusters.size();

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

  const float digClusRatio = static_cast<float>(1.0f*numDigits/numClusters);
  const float digTrigRatio = static_cast<float>(1.0f*numDigits/numTriggers);
  const float triggerFrequency = static_cast<float>(1.0f*numTriggers/durSec);

  cout << "digClusRatio" << digClusRatio << endl;
  cout << "digTrigRatio" << digTrigRatio << endl;

  const auto& ratioInfo = Form("Dig/Clus = %.2f Dig/Triggers= %.0f", digClusRatio, digTrigRatio); 

  const auto& trigInfo = Form("Triggers %i, Frequency [Hz]= %.2f " , numTriggers, triggerFrequency); 
  const auto& digClusInfo = Form("Digits %i Clusters %i",
              numDigits, numClusters);
  //const auto trigger = Form("First Entry %i Last %i", firstTrigger, lastTrigger);
  

  
  const auto durInfo = Form("Duration of triggers = %.2f min", durMin);
  return  {trigInfo, digClusInfo, ratioInfo, durInfo};
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
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.25);
  gStyle->SetStatFontSize(0.065);
  gStyle->SetLegendTextSize(0.08);//

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



void setPadChannel(bool (&padDigOff)[7][160][144], int chamber, int xLow, int xHigh, int yLow, int yHigh)
{
  
  if(xLow > xHigh || yLow > yHigh){ 
    cout << "Wrong Comparison" << endl;
    return;
  }

  if(xLow > 159 || xLow < 0 || xHigh > 159 || xHigh < 0 ){ 
    cout << "Wrong Comparison" << endl;
    return;
  } 

  if(yLow > 144 || yHigh < 0 || yLow > 144 || yHigh < 0){ 
    cout << "Wrong Comparison" << endl;
    return;
  } 

  for(int x = xLow; x <= xHigh; x++){
    for(int y = yLow; y <= yHigh; y++){
      padDigOff[chamber][x][y] = false;
      padDig[chamber][x][y] = false;
    }
  }
}






