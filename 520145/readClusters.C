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
#include <math.h>
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


vector<double> lastTimes, firstTimes, timeOfEvents; 
void sortTimes(vector<double>& timeOfEvents, int nEvents); //


double firstTrg, lastTrg;

TH1F* trigSort = new TH1F("Trigger Time Histogram", "Trigger Time Histogram", 50, 0., 1000000000.);
TH1F* trigSort2 = new TH1F("Trigger Frequency Histogram", "Trigger Frequency Histogram", 50, 0., 20000.);

void sortTriggers(vector<Trigger>& sortedTriggers);
//void sortTriggers(vector<Trigger>& sortedTriggers, TGraph& trigTimeSortStd);
double largestDiff = std::numeric_limits<double>::min();
double largestNegDiff = std::numeric_limits<double>::max();

vector<string> dig2Clus(const std::string &fileName, vector<Cluster>& clusters, vector<Trigger>& clusterTriggers, vector<Digit>& digits);

bool mReadFile = false;
std::string mSigmaCutPar;
float mSigmaCut[7] = {4, 4, 4, 4, 4, 4, 4};

std::unique_ptr<TFile> mFile; ///< input file containin the tree
std::unique_ptr<TTree> mTree; ///< input tree

int chargeBelow4[7][5];


std::unique_ptr<Clusterer> mRec; // ef: changed to smart-pointer
int trigTimeCount,trigTimeCount2, trigTimeCount3 = 0;

void initFileIn(const std::string &fileName);

void strToFloatsSplit(std::string s, std::string delimiter, float *res,
                      int maxElem = 7);


struct TriggerTimeInf
{
  double  timeInNs;
  int triggerIndex;

  TriggerTimeInf(double _timeInNs, int _triggerIndex)
    : timeInNs(_timeInNs)
    , triggerIndex(_triggerIndex)
  {}

};

vector<TriggerTimeInf> triggerInfoVec;

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
  
  std::array<std::unique_ptr<TH1F>, 7> digCharge, hMipCharge, digPerEvent, digCharges;
  std::array<std::unique_ptr<TH2F>, 7> digMap, digMapAvg, digMapSel, test, mapCharge4;


  std::array<std::unique_ptr<TH1F>, 3> triggerTimeFreqHist;
  TH1F trigSortHist;

  std::unique_ptr<TH1F> digPerEventFreq;
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

        // can not sort if only 1 trigger
        if(clusterTriggers.size() > 1) {
          sortTimes(timeOfEvents, nEvents);
        }

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
    return;
    std::exit(0);
  }

  float avgDigits = static_cast<float>(1.0f*digits.size()/numTriggers);
  
  for (int i = 0; i < 7; i++) {

    const char* canStringMip = Form("MIP-Charge %i", i);
    hMipCharge[i].reset(new TH1F(canStringMip, canStringMip, 50, 200., 2200.));
    //hMipCharge[i].reset(new TH1F(canStringMip, canStringMip, 500, 200., 2200.));
    hMipCharge[i]->SetXTitle("Charge (ADC channel)");
    hMipCharge[i]->SetYTitle("Entries/40 ADC");
    hMipCharge[i]->SetStats(kTRUE);

    const char* canStringSize = Form("Logaritmic Digit Charge %i", i);
    digCharge[i].reset(new TH1F(canStringSize, canStringSize, 50, 0., 400.));
    digCharge[i]->SetXTitle("Charge (ADC channel)");
    digCharge[i]->SetYTitle("Entries/40 ADC");
    digCharge[i]->SetLabelOffset(0.0065, "y");

    const char* canStringSizes = Form("Digit Charge %i", i);
    digCharges[i].reset(new TH1F(canStringSizes, canStringSizes, 10, 0., 10.));
    digCharges[i]->SetXTitle("Charge (ADC channel)");
    digCharges[i]->SetYTitle("Entries/40 ADC");
    digCharges[i]->SetLabelOffset(0.0065, "y");


    const char* canDigMap = Form("Digit Map %i", i);
    digMap[i].reset(new TH2F(canDigMap, canDigMap, 160, 0, 159, 144, 0, 143));
    digMap[i]->SetXTitle("x [cm]");
    digMap[i]->SetYTitle("y [cm]");

    const char* mapCharge4Str = Form("Chamber %i Digits with Charge < 4 ", i);
    mapCharge4[i].reset(new TH2F(mapCharge4Str, mapCharge4Str, 160, 0, 159, 144, 0, 143));
    mapCharge4[i]->SetXTitle("x [cm]");
    mapCharge4[i]->SetYTitle("y [cm]");
    

    const char* canDigSel = Form("Selected Pads %i", i);
    digMapSel[i].reset(new TH2F(canDigSel, canDigSel, 160, 0, 159, 144, 0, 143));
    //digMap[i].reset(new TH2F(canDigMap, canDigMap, 160*0.8, 0, 159*0.8, 144*0.8, 0, 160*0.8));
    digMapSel[i]->SetXTitle("x [cm]");
    digMapSel[i]->SetYTitle("y [cm]");

    const char* canDigAvg = Form("Avg Charge Per Pad %i", i);
    digMapAvg[i].reset(new TH2F(canDigAvg, canDigAvg, 160, 0, 159, 144, 0, 143));
    //digMap[i].reset(new TH2F(canDigMap, canDigMap, 160*0.8, 0, 159*0.8, 144*0.8, 0, 160*0.8));
    digMapAvg[i]->SetXTitle("x [cm]");
    digMapAvg[i]->SetYTitle("y [cm]");

    const char* digEvtFreqStr = Form("Digits Per Event Frequency%i",i);
    digPerEvent[i].reset(new TH1F(digEvtFreqStr, digEvtFreqStr, 500, 0., .5));
    //digPerEvent[i].reset(new TH1F(digEvtFreqStr, digEvtFreqStr, 500, 100*(avgDigits-150.)/(144*160), 100*(avgDigits+150.)/(144*160)));
    digPerEvent[i]->SetXTitle("Occupancy [%]");
    digPerEvent[i]->SetYTitle("Frequencies");

  }


  for(int chamber = 0; chamber < 7; chamber++){
    for(int x = 0; x < 160; x++){
      for(int y = 0; y < 144; y++){
        if(!padDig[chamber][x][y]){
          digMapSel[chamber]->Fill(x, y, 500000.);
          //cout << "False padDigOff " << chamber << " x " << x << " y " << y << endl;
        } else {
          digMapSel[chamber]->Fill(x, y, 10.);
        }
      }
    }
  }

  for(int chamber = 0; chamber < 7; chamber++){
    for(int x = 0; x < 160; x++){
      for(int y = 0; y < 144; y++){

        // verify that pads are selected off:
        if(digMapSel[chamber]->GetBinContent(x,y)==50.){
          // cout << "False padDigOff " << chamber << " x " << x << " y " << y << endl;
        }
      }
    }
  }


  for(int i = 0; i<3; i++){

    if(i==0){
     const char* trigTimeStr = Form("Trigger Time Freq%i",i);
     triggerTimeFreqHist[i].reset(new TH1F(trigTimeStr, trigTimeStr, 50, 0, largestDiff));
    } else if (i==1) {
     const char* trigTimeStr = Form("Event Delta-Time [nS]");
     triggerTimeFreqHist[i].reset(new TH1F(trigTimeStr, trigTimeStr, 10, 0, largestDiff));
     triggerTimeFreqHist[i]->SetXTitle("Event Delta-time [nS]");
    } else if (i==2) {
     const char* trigTimeStr = Form("Event Frequency");
     triggerTimeFreqHist[i].reset(new TH1F(trigTimeStr, trigTimeStr, 50, 0, 30000));
     triggerTimeFreqHist[i]->SetXTitle("Event Frequncy [Hz]");
    } 
    triggerTimeFreqHist[i]->SetYTitle("Number Of Entries");
  }
 
  const int nBins = static_cast<int>((lastTrg-firstTrg)*pow(10,-5));
  cout << "Number of Bins " << nBins << endl;
  cout << "1 " << firstTrg << " 2 " << lastTrg << endl;
  trigSort->SetBins(nBins, firstTrg, lastTrg);

  const char* trigEvtStr = Form("Graph of Event Delta-time; Time Event Occured [nS]; Event Delta Time [nS]");
  trigTime.reset(new TGraph);
  trigTime->SetTitle(trigEvtStr);  
  
  std::array<std::unique_ptr<TCanvas>, 8> canvas;  

  (canvas[0]).reset(new TCanvas(Form("MIP Cluster-Charge %i",fname), Form("MIP Cluster-Charge %i",fname),1200, 1200));  

  (canvas[1]).reset(new TCanvas(Form("Digit Charge Log%i",fname), Form("Digit Charge Log %i",fname), 1200, 1200));
  (canvas[1])->SetLogy();
  (canvas[1])->SetLogx();

  (canvas[2]).reset(new TCanvas(Form("Digit-Map %i",fname), Form("Digit-Map %i",fname), 1200, 1200));

  (canvas[3]).reset(new TCanvas(Form("Digit-Map Avg %i",fname), Form("Digit-Map Avg %i",fname), 1200, 1200));

  (canvas[4]).reset(new TCanvas(Form("Digits Per Event %i",fname), Form("Digits Per Event %i",fname), 1200, 1200));
  (canvas[4])->SetLogy();

  (canvas[5]).reset(new TCanvas(Form("PlaceHolder %i",fname), Form("PlaceHolder %i",fname),1200, 1200));
  (canvas[6]).reset(new TCanvas(Form("PlaceHolder2 %i",fname), Form("PlaceHolder2 %i",fname),1200, 1200)); 
  (canvas[7]).reset(new TCanvas(Form("PlaceHolder3 %i",fname), Form("PlaceHolder3 %i",fname),1200, 1200));  


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

      const auto& tDif = (trig.getIr()).differenceInBCNS(trigPrev.getIr());

      if(trigNum > 0 && time > pow(10,6)){
        const auto& freq = (pow(10, 9))/tDif;         
        triggerTimeFreqHist[2]->Fill(freq);
        triggerTimeFreqHist[1]->Fill(tDif);
        triggerTimeFreqHist[0]->Fill(tDif);
        trigTime->SetPoint(trigNum-1, static_cast<double>(time), tDif);
        trigSort->Fill(time); 
      }

      if(tDif > pow(10, 8)){
        cout << " Exceeded limit Trigger Number " << trigNum << " Tdiff " << tDif << endl ;  cout << "Count " << trigTimeCount3 << endl;
        trigTimeCount3++;

      }

      if(tDif > pow(10, 7)){
        cout << " Exceeded limit Trigger Number " << trigNum << " Tdiff " << tDif << endl ;  cout << "Count " << trigTimeCount2 << endl;
        trigTimeCount2++;

      }

      if(tDif > pow(10, 6)){
        cout << " Exceeded limit Trigger Number " << trigNum << " Tdiff " << tDif << endl ;  cout << "Count " << trigTimeCount << endl;
        trigTimeCount++;

      }

      if(time < pow(10, 6)){
        cout << "Time "<< time << " Trigger Number " << trigNum << " Tdiff " << tDif << endl ;
      }

      std::array<int, 7> cntCh = {0,0,0,0,0,0,0};

      for(int j = firstTrig; j < lastTrig; j++){
        const auto& dig = digits[j];
        Digit::pad2Absolute(dig.getPadID(), &module, &padChX, &padChY);
        //cout << "Fill Chamber, trNum " <<  module << " " << trigNum << endl;
        //cout << "Total " << cntCh[module] << endl;
	cntCh[module]++;
	avgDig[module]++;	
      }

      for(int ch = 0; ch < 7; ch++){
        digPerEvent[ch]->Fill(rel*cntCh[ch]);
      }

      tprev = time;
      trigPrev = trig;
      trigNum++; 

    }



    for(const auto& dig : digits){
      Digit::pad2Absolute(dig.getPadID(), &module, &padChX, &padChY);
      const auto& charge = dig.getQ();
      if(padDig[module][padChX][padChY] == true){
        digCharge[module]->Fill(charge);
        if(charge < 4){
          chargeBelow4[module][static_cast<int>(charge)] += 1; 
          mapCharge4[module]->Fill(padChX, padChY, charge);
        }
        if(charge <= 10){
          digCharges[module]->Fill(charge);
        }
      }

      digMap[module]->Fill(padChX, padChY, charge);
      digMapAvg[module]->Fill(padChX, padChY, charge/numTriggers);
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
  std::array<std::unique_ptr<TPaveText>, 4> tpvs2;

  //fileInfo
  const auto f1 = (fileInfo[0]).c_str();
  const auto f2 = (fileInfo[1]).c_str();
  const auto f3 = (fileInfo[2]).c_str();
  const auto f4 = (fileInfo[3]).c_str();

  //const char* runLabel = Form("%i  Duration = %s", fname, f1);
  const char* runLabel = Form("%i", fname);


  vector<const char*> tpvTexts{"MIP Clusters Charge", "Digits-Charge, logx logy", "Digits-Map", "Digits-Map Avg", "Digits Per Event", "Map of Evaluated Areas", "Digits Charge Small Scale"};

  vector<const char*> tpvTexts2{"Event Info", "Trigger Info",  "Map of Evaluated Areas","Charge Below 4"};


  int j = 0;
  for(auto& tpv: tpvs2){
    tpv.reset(new TPaveText(0.05, .05, .9, .9));
    tpv->AddText(Form("%s %s", runLabel, tpvTexts2[j++]));
  }

  for(auto& tpv: tpvs2){
    tpv->AddText(f1);
    tpv->AddText(f2);
    tpv->AddText(f3);
    tpv->AddText(f4);
  }

  j = 0;
  for(auto& tpv: tpvs){
    tpv.reset(new TPaveText(0.05, .05, .9, .9));
    tpv->AddText(Form("%s %s", runLabel, tpvTexts[j++]));
  }

  for(auto& tpv: tpvs){
    tpv->AddText(f1);
    tpv->AddText(f2);
    tpv->AddText(f3);
    tpv->AddText(f4);
  }

  const int tpvSize = static_cast<int>(tpvs.size());
  for(int i = 0; i < tpvSize -1; i++){
    canvas[i]->Divide(3, 3); 
    canvas[i]->cd(3);
    tpvs[i]->Draw();
  }

  const int posArr[] = {9, 8, 6, 5, 4, 2, 1};
  changeFont();    

  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.925);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3); 
    
  std::unique_ptr<TCanvas> temp1;
  temp1.reset(new TCanvas(Form("Event Information %i",fname), Form("Event Information %i",fname),1200, 2000));
  temp1->Divide(2,2);
  temp1->cd(2);
  tpvs2[0]->Draw();
    


  for (int i = 0; i < 3; i++)
  { 
    TPad* pad2;
    if(i < 1){ 
      pad2 = static_cast<TPad*>(temp1->cd(i+1));
      pad2->SetLeftMargin(.01+pad2->GetLeftMargin());
      //trigTime->SetTitleOffset(trigTime->GetTitleOffset("y")*0.8, "y");
      trigTime->Draw("A*");
    } else {
      pad2 = static_cast<TPad*>(temp1->cd(i+2));
      pad2->SetLeftMargin(.0275+pad2->GetLeftMargin());
      pad2->SetBottomMargin(.0375+pad2->GetBottomMargin());
      pad2->SetRightMargin(.0375+pad2->GetRightMargin());
      triggerTimeFreqHist[i]->SetTitleOffset(triggerTimeFreqHist[i]->GetTitleOffset("y")*1.2, "xy");
       
      triggerTimeFreqHist[i]->SetTitleSize(triggerTimeFreqHist[i]->GetTitleSize("x")*0.95, "xy");
      triggerTimeFreqHist[i]->SetLabelSize(triggerTimeFreqHist[i]->GetLabelSize("x")*0.925, "x");
      triggerTimeFreqHist[i]->SetLabelSize(triggerTimeFreqHist[i]->GetLabelSize("y")*0.925, "y");
      triggerTimeFreqHist[i]->Draw();
    }

  }

  gStyle->SetOptStat("eim");
  gStyle->SetStatX(0.95);
  temp1->Show();
  temp1->SaveAs(Form("Event Information_%i_.png",fname));



  /*
   **********************************
   Sorted Triggers   
   **********************************
  */

  std::unique_ptr<TCanvas> temp2;
  temp2.reset(new TCanvas(Form("Trigger Frequency%i",fname), Form("Trigger Frequency%i",fname),1200, 2000));
  temp2->Divide(2,2);
  temp2->cd(2);
  tpvs2[1]->Draw();
  
  auto pad5 = static_cast<TPad*>(temp2->cd(1));

  //trigTime->SetMinimum(pow(10,12));
  //trigTime->Draw("AC*");
  trigTime->Draw("A*");




  /* Change size of trigSort ?
  triggerTimeFreqHist[i]->SetTitleOffset(triggerTimeFreqHist[i]->GetTitleOffset("y")*1.2, "xy");
  triggerTimeFreqHist[i]->SetTitleSize(triggerTimeFreqHist[i]->GetTitleSize("x")*0.95, "xy");
  triggerTimeFreqHist[i]->SetLabelSize(triggerTimeFreqHist[i]->GetLabelSize("x")*0.925, "x");
  triggerTimeFreqHist[i]->SetLabelSize(triggerTimeFreqHist[i]->GetLabelSize("y")*0.925, "y");*/
 
  for(int i = 0; i < 2; i++){
    TPad* pad = static_cast<TPad*>(temp2->cd(3+i));
    pad->SetLeftMargin(.0275+pad->GetLeftMargin());
    pad->SetBottomMargin(.0375+pad->GetBottomMargin());
    pad->SetRightMargin(.0375+pad->GetRightMargin());

    if(i==0){
      trigSort->SetXTitle("Event Time in LHC nS");
      trigSort->SetYTitle("Numer of Entries"); 
      trigSort->Draw();
    } else {
      trigSort2->SetXTitle("Event Frequency");
      trigSort2->SetYTitle("Numer of Entries");
      trigSort2->Draw();
    }
  }
  
  temp2->SaveAs(Form("Trigger Frequency Hist and Graph %i.png",fname));


  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== Digit MAP =========================


    auto pad5 = static_cast<TPad*>(canvas[2]->cd(pos));
    //digPerEvent[iCh]->Draw();
    pad5->SetLeftMargin(+.025+pad5->GetLeftMargin());


    const auto& pTotalDigs = static_cast<float>(100.0f*digMap[iCh]->GetEntries()/digSize);
    
    digMap[iCh]->SetLabelOffset(digMap[iCh]->GetLabelOffset("y")-0.0015, "y");
    digMap[iCh]->SetTitleOffset(digMap[iCh]->GetTitleOffset("y")-0.0015, "y");
    digMap[iCh]->SetTitleOffset(digMap[iCh]->GetTitleOffset("x")-0.0005, "x");

    pad5->SetBottomMargin(.0015+pad5->GetBottomMargin());
    pad5->SetRightMargin(.125+pad5->GetRightMargin());
    digMap[iCh]->SetTitle(Form("Chamber %i Percentage of total = %02.0f", iCh, pTotalDigs));
    digMap[iCh]->SetMarkerStyle(3);
    digMap[iCh]->Draw("Colz");

    digMap[iCh]->SetStats(kFALSE);
  }

  gStyle->SetStatX(0.85);
  gStyle->SetOptStat("e");
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6); 

  // avg digits charge
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== Digit MAP =========================
    auto pad5 = static_cast<TPad*>(canvas[3]->cd(pos));
    //pad5->SetLeftMargin(+.025+pad5->GetLeftMargin());
    const auto& pTotalDigs = static_cast<float>(100.0f*digMapAvg[iCh]->GetEntries()/digSize);

    pad5->SetBottomMargin(.0015+pad5->GetBottomMargin());
    pad5->SetRightMargin(.125+pad5->GetRightMargin());
    digMapAvg[iCh]->SetTitle(Form("Chamber Avg %i Percentage of total = %02.0f", iCh, pTotalDigs));
    digMapAvg[iCh]->SetStats(kFALSE);
    digMapAvg[iCh]->Draw("Colz");
  }
    
  gStyle->SetOptStat("e");
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6);     
  gStyle->SetStatX(0.85);

  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== MIP Charge =========================
    auto pad0 = static_cast<TPad*>(canvas[0]->cd(pos));
    pad0->SetLeftMargin(-.0025+pad0->GetLeftMargin());
    pad0->SetRightMargin(-.005+pad0->GetRightMargin());
    pad0->SetBottomMargin(.0025+pad0->GetBottomMargin());
    hMipCharge[iCh]->Fit("landau", "I"); // I = fit by integral
    hMipCharge[iCh]->SetLabelOffset(hMipCharge[iCh]->GetLabelOffset("y")-0.0025, "y");
    hMipCharge[iCh]->SetTitleOffset(0.8, "y");
    //hMipCharge[iCh]->SetTitle(Form("Constant %03.1f \n MPV %03.1f Sigma %03.1f", Constant, MPV, Sigma));
    hMipCharge[iCh]->Draw();
  }

  gStyle->SetStatX(0.95);
  //drawMipCharge(hMipCharge)

  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6); 
  gStyle->SetOptStat("eimr");
  gStyle->SetLabelOffset(0.00525, "y");

  // avg digits charge
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== Digit MAP =========================
    auto pad5 = static_cast<TPad*>(canvas[4]->cd(pos));
    pad5->SetLeftMargin(+.025+pad5->GetLeftMargin());
    pad5->SetLogy(1);
    pad5->SetBottomMargin(.0015+pad5->GetBottomMargin());
    pad5->SetRightMargin(-.0025+pad5->GetRightMargin());
    digPerEvent[iCh]->SetTitleOffset(digPerEvent[iCh]->GetTitleOffset("y")+0.025, "y");
    digPerEvent[iCh]->Draw();
  }
  gStyle->SetStatX(0.95);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.6); 
  gStyle->SetOptStat("eimr");
  gStyle->SetLabelOffset(0.00525, "y");


  std::unique_ptr<TCanvas> digMapSelCanv;
  digMapSelCanv.reset(new TCanvas(Form("Pads turned off by user%i",fname), Form("Pads turned off by user%i",fname), 1200, 1200));
  digMapSelCanv->Divide(3,3);

  digMapSelCanv->cd(3);
  tpvs2[2]->Draw();
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== Digit Charge =========================
    auto pad3 = static_cast<TPad*>(digMapSelCanv->cd(pos));
    pad3->SetBottomMargin(.0025+pad3->GetBottomMargin());
    pad3->SetLeftMargin(.065+pad3->GetLeftMargin());
    digMapSel[iCh]->SetLabelOffset(digMapSel[iCh]->GetLabelOffset("y")+0.0015, "y");
    digMapSel[iCh]->SetTitleOffset(1.3,"y");
    pad3->SetRightMargin(-.0025+pad3->GetRightMargin());
    digMapSel[iCh]->SetMarkerStyle(3);
    digMapSel[iCh]->SetStats(kFALSE);
    digMapSel[iCh]->Draw("Colz");
  }

  digMapSelCanv->Show();


  std::unique_ptr<TCanvas> digMapLowCan;
  digMapLowCan.reset(new TCanvas(Form("Charge below 4 %i",fname), Form("Charge Below 4 %i",fname), 1200, 1200));
  digMapLowCan->Divide(3,3);
  digMapLowCan->cd(3);
  tpvs2[3]->Draw();
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== Digit Charge =========================
    auto pad3 = static_cast<TPad*>(digMapLowCan->cd(pos));
   
    pad3->SetBottomMargin(.0025+pad3->GetBottomMargin());
    pad3->SetLeftMargin(.065+pad3->GetLeftMargin());
    mapCharge4[iCh]->SetLabelOffset(mapCharge4[iCh]->GetLabelOffset("y")+0.0015, "y");
    mapCharge4[iCh]->SetTitleOffset(1.3,"y");
    pad3->SetRightMargin(-.0025+pad3->GetRightMargin());
    mapCharge4[iCh]->SetStats(kFALSE);
    mapCharge4[iCh]->SetMarkerStyle(3);
    mapCharge4[iCh]->Draw("Colz");
  }

  digMapLowCan->Show();
  digMapLowCan->SaveAs(Form("Digit Map of Low Charge%i.png",fname));

  std::unique_ptr<TCanvas> t;
  t.reset(new TCanvas(Form("Digits SmallRange%i",fname), Form("Digits SmallRange %i",fname), 1200, 1200));
  t->Divide(3,3);
  t->cd(3);
  (tpvs.back())->Draw();

  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];
    // ========== Digit Charge =========================
    auto pad3 = static_cast<TPad*>(t->cd(pos));
    pad3->SetBottomMargin(.0025+pad3->GetBottomMargin());
    pad3->SetLeftMargin(.065+pad3->GetLeftMargin());
    digCharges[iCh]->SetLabelOffset(digCharges[iCh]->GetLabelOffset("y")+0.0015, "y");
    digCharges[iCh]->SetTitleOffset(1.3,"y");
    pad3->SetRightMargin(-.0025+pad3->GetRightMargin());
    digCharges[iCh]->Draw();
  }  

  for (int iCh = 0; iCh < 7; iCh++) {
    const auto& pos = posArr[iCh];

    // ========== Digit Charge =========================
    auto pad3 = static_cast<TPad*>(canvas[1]->cd(pos));
    (canvas[1])->SetLogy();
    (canvas[1])->SetLogx();
    pad3->SetLogy(1);
    pad3->SetLogx(1);
    pad3->SetBottomMargin(.025+pad3->GetBottomMargin());
    pad3->SetLeftMargin(.02+pad3->GetLeftMargin());
    digCharge[iCh]->SetLabelOffset(digCharge[iCh]->GetLabelOffset("y")+0.0015, "y");
    digCharge[iCh]->SetTitleOffset(1.05,"y");
    digCharge[iCh]->SetTitleOffset(1.1,"x");
    digCharge[iCh]->SetTitleSize(digCharge[iCh]->GetTitleSize("x")*1.05, "xy");
    digCharge[iCh]->SetLabelSize(digCharge[iCh]->GetLabelSize("x")*1.05, "xy");

    pad3->SetRightMargin(-.0025+pad3->GetRightMargin());
   
    digCharge[iCh]->Draw();

   for(int charge = 0; charge <= 4; charge++){
      //cout << " Ch, Charge, num " << iCh << " " << charge << " " << chargeBelow4[iCh][charge] << endl; 
   }
  }
  cout << "trigTimeCount Ended at " << trigTimeCount<< endl;
  cout << "trigTimeCount2 Ended at " << trigTimeCount2<< endl;
  cout << "trigTimeCount3 Ended at " << trigTimeCount3<< endl;
  gStyle->SetStatH(0.2); 
  gStyle->SetStatX(0.95);
  gStyle->SetOptStat("eim");
  gStyle->SetLabelOffset(0.008, "y");


  canvas[0]->SaveAs(Form("MIP_Cluster_Charge_%i_.png",fname));

  canvas[1]->SaveAs(Form("Digit_Charge_Log%i_.png",fname));
  canvas[2]->SaveAs(Form("Digit_Map_%i_.png",fname));
  canvas[3]->SaveAs(Form("Digit_Map_Avg_%i_.png",fname));
  canvas[4]->SaveAs(Form("DigitsPerEvent_%i_.png",fname));

  t->SaveAs(Form("Digit_Charge_SmallRange%i_.png",fname));
  digMapLowCan->SaveAs(Form("Digits below 4 %i_.png",fname));
  canvas[0]->Show();
  canvas[1]->Show();
  t->Show();
  canvas[2]->Show();
  canvas[3]->Show();
  canvas[4]->Show();


  sleep_for(5000ms);


  //return;
  //std::exit(0);


  bool userInput = false;
  while(!userInput){
    sleep_for(5000ms);
   

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
      canvas[0]->Close();canvas[1]->Close();canvas[2]->Close();
      canvas[3]->Close();canvas[4]->Close();cout << "Got End fro User.. Exiting!";
      canvas[0]->SaveAs(Form("MIP_Cluster_Charge_%i_.png",fname));
      canvas[1]->SaveAs(Form("Digit_Charge_%i_.png",fname));
      canvas[2]->SaveAs(Form("Digit_Map_%i_.png",fname));
    }
  }
}

vector<string> dig2Clus(const std::string &fileName, vector<Cluster>& clusters, vector<Trigger>& clusterTriggers, vector<Digit>& digits)
{
  long mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
  uint32_t firstTrigger, lastTrigger = 0;
  
  mRec.reset(new o2::hmpid::Clusterer()); // ef: changed to smart-pointer


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
        
    // isClustersAlready();
    
    for (const auto &trig : *mTriggersFromFilePtr) {
      if (trig.getNumberOfObjects()) {
        gsl::span<const Digit> trigDigits{
            mDigitsFromFilePtr->data() + trig.getFirstEntry(),
            size_t(trig.getNumberOfObjects())};
        const size_t clStart = clusters.size();
        mRec->Dig2Clu(trigDigits, clusters, mSigmaCut, true); //ef:uncomment
        clusterTriggers.emplace_back(trig.getIr(), clStart,
                                     clusters.size() - clStart);
      }
    }           

    cout << " Received " << mTriggersFromFilePtr->size() << " triggers with " << mDigitsFromFilePtr->size() << " digits -> clusters = " << clusters.size() << endl;

    digits = *mDigitsFromFilePtr;
    if(digits.size() == 0){
      digits = mDigitsFromFile;
    }
    mDigitsReceived = mDigitsFromFilePtr->size();
    mClustersReceived = clusters.size();
    mTriggersReceived = mTriggersFromFilePtr->size();
  }

  cout << "ef marker remove" << endl ;
  const int numTriggers = static_cast<int>(mTriggersReceived);
  const int numDigits = static_cast<int>(mDigitsReceived);
  const int numClusters = static_cast<int>(mClustersReceived);

  const float digClusRatio = static_cast<float>(1.0f*numDigits/numClusters);
  const float digTrigRatio = static_cast<float>(1.0f*numDigits/numTriggers);
  const auto& ratioInfo = Form("Dig/Clus = %.2f Dig/Events= %.0f", digClusRatio, digTrigRatio); 

   
  const auto& digClusInfo = Form("Digits %i Clusters %i",
              numDigits, numClusters);

  if(mDigitsFromFilePtr->size() < 2) {
    const float triggerFrequency = static_cast<float>(1.0f*numTriggers/durSec);
    const auto& trigInfo = Form("Events %i, Average Frequency [Hz] = %.2f " , numTriggers, triggerFrequency);
    return  {trigInfo, digClusInfo, ratioInfo, Form("Not enough Events (%i) for Frequency", numTriggers)};
  }

  // sort triggers by time
  sortTriggers(*mTriggersFromFilePtr);

  Trigger trigFirst = mTriggersFromFilePtr->at(0);
  Trigger trigLast = mTriggersFromFilePtr->back();   

  const auto& tDif = (trigLast.getIr()).differenceInBCNS(trigFirst.getIr());

  durSec = static_cast<double>((tDif)/1000000000.0);
  durMin = static_cast<double>((durSec)/60.0);
  const auto durInfo = Form("Duration of Events = %.2f min", durMin);
  cout << "digClusRatio " << digClusRatio << endl;
  cout << "digTrigRatio " << digTrigRatio << endl;
  const float triggerFrequency = static_cast<float>(1.0f*numTriggers/durSec);
  const auto& trigInfo = Form("Events %i, Frequency [Hz] = %.2f " , numTriggers, triggerFrequency); 
  
  int trigNum = 0;
  Trigger trigPrev;

  const auto& fTrig = mTriggersFromFile[0];
  const auto& lTrig = mTriggersFromFilePtr->back();
  
  const auto& irFirst = fTrig.getIr();
  const auto& irLast = lTrig.getIr();

  const auto& orbitFirst = fTrig.getOrbit();
  const auto& bcFirst = fTrig.getBc();

  const auto& orbitLast = lTrig.getOrbit();
  const auto& bcLast = lTrig.getBc();

  const auto& nsFirst = InteractionRecord::bc2ns(bcFirst, orbitFirst);
  const auto& nsLast = InteractionRecord::bc2ns(bcLast, orbitLast);

  const auto prevTimeNsTrig = nsFirst;
  int bcPrev; unsigned int orbitPrev;

  cout << "First, Ns " << nsFirst << " BC " << bcFirst  << " Orbit " << orbitFirst << endl;
  cout << "Last , Ns " << nsLast <<  " BC " << bcLast  << " Orbit " << orbitLast << endl;

  for(const auto &trig : *mTriggersFromFilePtr){
    const int numDigPerTrig = trig.getNumberOfObjects();
    const int firstTrig = trig.getFirstEntry();
    const int lastTrig = trig.getLastEntry();

    auto tDelta = (trig.getIr()).differenceInBCNS(trigPrev.getIr());
   
    const auto& orbit = trig.getOrbit();
    const auto& bc = trig.getBc();
    const auto timeTrigger = InteractionRecord::bc2ns(bc, orbit);
    
    triggerInfoVec.emplace_back(timeTrigger, trigNum);
    trigPrev = trig;
    if(trigNum > 0){
      tDelta = InteractionRecord::bc2ns(bc, orbit) - InteractionRecord::bc2ns(bcPrev, orbitPrev);

      // just to check theyre equal:
      //cout << "diff trig ir "<< (trig.getIr()).bc2ns() - timeTrigger << endl;
    }
    trigNum++;
    bcPrev = bc; orbitPrev = orbit;
  }
  
  cout << " largest difference " << largestDiff << endl;
  cout << " largest negative   " << largestNegDiff << endl;
  

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

    LOG(warn)
        << "HMPID DigitToClusterSpec::init() : Did not find o2sim tree in "
        << filename.c_str() << endl;
    return;
    std::exit(0);
  }

  if ((mTree->GetBranchStatus("HMPDigit")) == 1) {
    mTree->SetBranchAddress("HMPDigit", &mDigitsFromFilePtr);
  } else if ((mTree->GetBranchStatus("HMPIDDigits")) == 1) {
    mTree->SetBranchAddress("HMPIDDigits", &mDigitsFromFilePtr);
  } else {
   LOG(warn)
        << "HMPID DigitToClusterSpec::init() : Error in branches!" << endl;
    return;
    std::exit(0);
  }

  mTree->SetBranchAddress("InteractionRecords", &mTriggersFromFilePtr);
  mTree->Print("toponly");
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

void sortTimes(vector<double>& timeOfEvents, int nEvents) // lastTimes firstTimes
{ 

  std::sort(timeOfEvents.begin(),  timeOfEvents.end(), [](const auto& a, const auto& b)
  {
    return (a < b);
  }); 

  const int timeLimit = timeOfEvents.size() - nEvents;
  int i = 0;
  for(const auto& time : timeOfEvents){
    if(i < timeLimit) {firstTimes.emplace_back(time); }
    else {lastTimes.emplace_back(time); } 

    i++;
  }

  double avgTimeF, avgTimeL;
  for(const auto& f: firstTimes){
    //cout << "time : " << f << endl;
    avgTimeF += f;
  } avgTimeF = avgTimeF/timeLimit;

  for(const auto& f: lastTimes){
    cout << "time : " << f << endl;
    avgTimeF += f;
  } avgTimeL = avgTimeL/nEvents;
  if(firstTimes.size() > 2) {
    cout << "firstTimes[timeLimit-2]" << firstTimes[timeLimit-2] << endl;
    cout << " min min2 nax times First : " << firstTimes[0]<< " " << firstTimes[1] << " " << firstTimes.back() <<endl;
  }

  if(lastTimes.size() > 2) {
    cout << " min min2 nax times Last : " << lastTimes[0] << " " << lastTimes[1] << " " << lastTimes.back() <<endl;
  }
  cout << " avg times First : " << avgTimeF << " Last : " << avgTimeL << endl;
}

void sortTriggers(vector<Trigger>& sortedTriggers)
{
  cout << "Sorting triggers length =  " << sortedTriggers.size() << endl;

  std::sort(sortedTriggers.begin(),  sortedTriggers.end(), [](const auto& a, const auto& b)
  {
    return (a.getIr()).bc2ns() < (b.getIr()).bc2ns();
  });

  cout << "Iterating through Sorted triggers " << endl;  

  int trigNum = 0;
  Trigger trigPrev;
   
  const auto firstTrig = sortedTriggers[0];
  const auto lastTrig = sortedTriggers.back();
  const auto& tDifTotal = (lastTrig.getIr()).differenceInBCNS(firstTrig.getIr());
  
  cout << "Diff Between first and last Event " << tDifTotal << endl;  
  cout << " Avg Frequency = " << sortedTriggers.size()/(tDifTotal*pow(10,9)) << endl;

  firstTrg = (firstTrig.getIr()).bc2ns();
  lastTrg = (lastTrig.getIr()).bc2ns();
  for(const auto& trig : sortedTriggers){
    
    const auto& tS = (trig.getIr()).bc2ns();
    const auto& tE = (trigPrev.getIr()).bc2ns();

    const auto& tDif2 = tE-tS;
    const auto& tDif = (trig.getIr()).differenceInBCNS(trigPrev.getIr());


    //cout << "  tDif " << tDif << endl;
    //cout << "  tDif2 " << tDif2 << endl;

    if(trigNum > 0){
      timeOfEvents.emplace_back(tDif);
      const auto& freq = (pow(10, 9))/tDif;     
      trigSort2->Fill(freq);

      if(tDif>largestDiff){
        largestDiff = tDif;
      }
      if(tDif<largestNegDiff){
        largestNegDiff = tDif;
      }
    }
    trigPrev = trig;
    trigNum++;
  }

  cout << "Triggers Sorted" << endl;  
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


void changeFont()
{
  /*
  std::unique_ptr<TStyle> mStyle; 
  mStyle.reset(new TStyle("canvasStyle", "Canvas Root Styles"));
  */ 
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.925);
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



