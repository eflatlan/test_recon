#if !defined(__CLING__) || defined(__ROOTCLING__)
//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TList.h>
#include <TLine.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TGraph.h>
#include <fstream>
#include <vector>
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include "DataFormatsHMP/Cluster.h"
#include <TROOT.h> // gRoot

// C++ header files and libraries
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <chrono>
#include <ctime>    
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#endif

  void SaveFolder(string inpFile);
  void changeFont();
  bool mReadFile = false;
  std::string mSigmaCutPar;
  float mSigmaCut[7] = {4, 4, 4, 4, 4, 4, 4};

  std::unique_ptr<TFile> mFile; ///< input file containin the tree
  std::unique_ptr<TTree> mTree; ///< input tree

  std::vector<o2::hmpid::Digit> mDigitsFromFile, *mDigitsFromFilePtr = &mDigitsFromFile;
  std::vector<o2::hmpid::Trigger> mTriggersFromFile, *mTriggersFromFilePtr = &mTriggersFromFile;

  std::unique_ptr<o2::hmpid::Clusterer> mRec; // ef: changed to smart-pointer
  long mDigitsReceived;
  long mClustersReceived;

  void initFileIn(const std::string& fileName);


  void strToFloatsSplit(std::string s, std::string delimiter, float* res,
                        int maxElem = 7);

void readClusters(int nEvents)
{ 


  TH1F *hCharge[7], *hMipCharge[7], *hSize[7];
  TH2F *hMap[7];
      
  for(int i=0; i<7; i++) {
     hMap[i] = new TH2F(Form("Clusters Map chamber%i",i),Form("Cluster Map chamber%i",i), 160, 0, 159, 144, 0, 143);
     hMap[i]->SetXTitle("X (cm)");
     hMap[i]->SetYTitle("Y (cm)"); 
     
     hCharge[i] = new TH1F(Form("Clusters Charge chamber%i",i),Form("Cluster Charge chamber%i",i),2000, 100., 2100.);
     hCharge[i]->SetXTitle("Charge (ADC channel)");               
     hCharge[i]->SetYTitle("Entries");
     
     hMipCharge[i] = new TH1F(Form("Mip Clusters Charge chamber%i",i),Form("Mip Cluster Charge chamber%i",i),50, 200., 2200.);
     hMipCharge[i]->SetXTitle("Charge (ADC channel)");               
     hMipCharge[i]->SetYTitle("Entries/40 ADC");
     hMipCharge[i]->SetLineColor(kBlack);
     
     hSize[i] = new TH1F(Form("Cluster Size chmaber%i",i),Form("Cluster Size chamber%i",i),20, 0., 20.);
     hSize[i]->SetXTitle("Cluster size");               
     hSize[i]->SetYTitle("Entries");               
  }

  //changeFont();			     // specify folder to save files in
  //SaveFolder("clusterChambers");   // apply custom canvas figure options

  auto folderName = (gSystem->GetWorkingDirectory());
  auto gwd = gSystem->GetWorkingDirectory();

  Printf("RunNumber %s", gSystem->GetDirName());

  auto runNumber = (gwd.substr(gwd.length()-9, gwd.length()));
  /*
  Printf("RunNumber %s", runNumber);
  Printf("RunNumber %s", runNumber.c_str());

  Printf("gwd %s",gwd.c_str());
  Printf("pwd %s",gSystem->pwd());
  Printf("homedir %s",gSystem->HomeDirectory()); */
  Printf("Empty clusters");
  
  TCanvas *c1 = new TCanvas("c1",("Cluster-Map " + runNumber).c_str(),2000,1200); 	
  TCanvas *c2 = new TCanvas("c2", ("Cluster-Charge " + runNumber).c_str(),2000,1200); 	
  TCanvas *c3 = new TCanvas("c3", ("MIP Cluster-Charge " + runNumber).c_str(),2000,1200); 	
  TCanvas *c4 = new TCanvas("c4", ("Cluster-Size " + runNumber).c_str(),2000,1200);  
  /*
  TCanvas *c1 = new TCanvas("c1","c1",1000,800); 
  TCanvas *c2 = new TCanvas("c2","c2",1000,800); 
  TCanvas *c3 = new TCanvas("c3","c3",1000,800); 
  TCanvas *c4 = new TCanvas("c4","c4",1000,800); */ 
   
  Int_t pos[] = {9,8,6,5,4,2,1};
  
  Int_t nTotTriggers = 0;
  


  for(int k = 0; k<nEvents; k++) {
    
    Printf("hmpclus%02i.root", k+1);
    std::unique_ptr<TFile> fileClusters{ TFile::Open(Form("hmpclus%02i.root", k+1))};
    std::unique_ptr<TTree> treeClusters;
    treeClusters.reset((TTree*)fileClusters->Get("o2sim"));
    o2::hmpid::Trigger *pTgr;
    o2::hmpid::Cluster *pClu;
    o2::hmpid::Cluster *pCluEvt;
    o2::hmpid::Cluster cluster;
  
    std::vector<o2::hmpid::Cluster> *clusters = nullptr;
    std::vector<o2::hmpid::Cluster> oneEventClusters;
    std::vector<o2::hmpid::Trigger> *trigger = nullptr;
  
    if(!treeClusters) {
      Printf("Empty clusters");
      continue;
    } else {
       treeClusters->SetBranchAddress("HMPIDClusters",&clusters);
       treeClusters->SetBranchAddress("InteractionRecords",&trigger);
       Printf("Got clusters");
    }
    

    const int treeClusSize = treeClusters->GetEntries();
    Printf("tree entries = %i", treeClusSize);
    
    treeClusters->GetEntry(0);


    const int clusterSize = clusters->size();
    Printf("clusters size = %i", clusterSize);
  
    int module = 0;
    
    for(int j = 0; j < clusterSize; j++) {
      
      pClu = (o2::hmpid::Cluster*)&clusters[j];

      module = pClu->ch();
                             
      hCharge[module]->Fill(pClu->q());
    
      if(pClu->size() >=3 && pClu->size()<=7) hMipCharge[module]->Fill(pClu->q());
    
      hSize[module]->Fill(pClu->size());
   }
       
   nTotTriggers+=trigger->size();
   const int triggerSize = trigger->size();
   Printf("trigger size from clusters = %i", triggerSize); 
   
   for(int i = 0; i < triggerSize; i++) {
     oneEventClusters.clear();
     pTgr = static_cast<o2::hmpid::Trigger*> (&trigger[i]);

     const int pTrgrFirst = pTgr->getFirstEntry();
     const int pTrgrLast = pTgr->getLastEntry();
     // for(int j = pTrgrFirst; j <= pTrgrLast; j++) {
     for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) { 
       //cluster = static_cast<o2::hmpid::Cluster>(clusters->at(j)); // ->at(j) to [j]
       cluster = static_cast<o2::hmpid::Cluster> (clusters->at(j)); // ->at(j) to [j]
       oneEventClusters.push_back(cluster);
     }              
   } // exit for2s

  } // exit for3

  c1->Divide(3,3); c2->Divide(3,3);
  c3->Divide(3,3); c4->Divide(3,3);

  for(int iCh = 0; iCh<7; iCh++){
     c1->cd(pos[iCh]); 
     hMap[iCh]->SetMarkerStyle(3);
     hMap[iCh]->Draw();
     
     c2->cd(pos[iCh]); 
     hCharge[iCh]->Draw();
     
     c3->cd(pos[iCh]); 
     hMipCharge[iCh]->Draw();
     
     c4->cd(pos[iCh]); 
     hSize[iCh]->Draw();
  }
  //Printf("Number of triggers = %i", nTotTriggers);     
  //Printf(

  c1->SaveAs(("clusterMap_" +  runNumber + "_.eps").c_str());
  c2->SaveAs(("clusterCharge_"+ runNumber+ "_.eps").c_str());
  c3->SaveAs(("mipClustesCharge_"+  runNumber+ "_.eps").c_str());
  c4->SaveAs(("clusterSize_" + runNumber+ "_.eps").c_str());

  c1->SaveAs(("clusterMap_" +  runNumber + "_.png").c_str());
  c2->SaveAs(("clusterCharge_"+ runNumber+ "_.png").c_str());
  c3->SaveAs(("mipClustesCharge_"+  runNumber+ "_.png").c_str());
  c4->SaveAs(("clusterSize_" + runNumber+ "_.png").c_str());
}    
   

//********************************************************************************************************************

void readDigits(char* filename, int nEvent)
{
  TFile *fileDigits = TFile::Open(filename);
  // Cast to TTree*, get tree by key "o2sim"
  TTree *treeDigits = (TTree*)fileDigits->Get("o2sim"); 
  
  // Initialize an array of TH1F-pointers 
  // (TH1=1D histogram, F specifies Float-values )
  TH1F *hCharge[7], *h_xCoord[7], *h_yCoord[7];    //  number of digits per y-coordinate
  TH2F *hMap[7];	// 2d map of digits
      
  // Label the histograms
  for(int i=0; i<7; i++) {
    // define element number i in the pointer-array

    hMap[i] = new TH2F(Form("Digits Map %i",i),\
    Form("Digits Map %i",i), 160, 0, 159, 144, 0, 143);  
    hMap[i]->SetXTitle("pad X [cm]");
    hMap[i]->SetYTitle("pad Y [cm]"); 
   
    hCharge[i] = new TH1F(Form("Digits Charge %i",i),\
    Form("Digits Charge %i",i),2000, 100., 2100.);
    hCharge[i]->SetXTitle("Charge (ADC channel)");               
    hCharge[i]->SetYTitle("Entries");         

    h_xCoord[i] = new TH1F(Form("Digits X-location Histogram %i",i),\
    Form("Digits X-location Histogram %i",i),2000, 10., 159.);
    h_xCoord[i]->SetXTitle("X [cm]");               
    h_xCoord[i]->SetYTitle("Entries");  


    h_yCoord[i] = new TH1F(Form("Digits Y-location Histogram %i",i),\
    Form("Digits Y-location Histogram %i",i),2000, 10., 144.);
    h_yCoord[i]->SetXTitle("Y [cm]");               
    h_yCoord[i]->SetYTitle("Entries");    
     
  }
  
  changeFont();		         // apply custom canvas figure options
  //SaveFolder("digitChambers"); // specify folder to save files in

  // Define canvases for plotting the figures  
  TCanvas *c1 = new TCanvas("c1","c1",2000,1200); 	
  TCanvas *c2 = new TCanvas("c2","c2",2000,1200); 	
  TCanvas *c3 = new TCanvas("c3","c3",2000,1200); 	
  TCanvas *c4 = new TCanvas("c4","c4",2000,1200); 	

  // Define positions for the plots for the chambers in the canvases 
  Int_t pos[]= {9,8,6,5,4,2,1};
  
  o2::hmpid::Trigger *pTgr;  // pointer to Trigger-object
  o2::hmpid::Digit *pDig;    // pointer to Digit-object
  o2::hmpid::Digit *pDigEvt; // pointer to Digit-object
  o2::hmpid::Digit digit;    // declaration of digit-object
  
  std::vector<o2::hmpid::Digit> *digits = nullptr;    //  vector of digit-pointers
  std::vector<o2::hmpid::Digit> oneEventDigits;       //  vector of Digit-objects
  std::vector<o2::hmpid::Trigger> *trigger = nullptr; //  vector of trigger-pointers
  
  // specify branches of tree to read from 
  // (const char* branchname, void* adress)
  treeDigits->SetBranchAddress("HMPDigit",&digits); 
  treeDigits->SetBranchAddress("InteractionRecords",&trigger);

  // Number of entries in tree
  Printf("tree entries = %i", treeDigits->GetEntries()); 
  
  treeDigits->GetEntry(0);

  Printf("digit size = %i", digits->size());
  
  // initialize padhit x and y coordinates, 
  // and module (i.e.) chamber-number
  int padChX = 0, padChY = 0, module = 0;

  // Loop through digits in file
  for(unsigned int j = 0; j < digits->size(); j++) {
      
    pDig = (o2::hmpid::Digit*)&digits->at(j); 
    // digits->at(j) = use arrow-operator since it is an array of pointers
    //&digits = create pointer pDig to adress of vector-element j in digits 
    // (o2::hmpid::Digit*) = cast to Digit*
    
    // Call member-function pad2Absolute in Digit-class
    // static void pad2Absolute(uint32_t pad, int* Module, int* x, int* y);
    // No need to be called from an object since the member function is static
    o2::hmpid::Digit::pad2Absolute(pDig->getPadID(), &module, &padChX, &padChY);    
    
    // Fill histograms, with the module specifying which chamber, 
    //(= pad in the canvas), should be filled 
    hCharge[module]->Fill(pDig->getQ()); // getQ gets charge
    h_xCoord[module]->Fill(padChX);	     // fill x-coordinate 
    h_yCoord[module]->Fill(padChY);	     // fill y-coordinate
  }
       
  Printf("trigger size from digits = %i", trigger->size()); 
   
	
  // Loop through triggers 
  for(int i = 0; i < trigger->size(); i++) {
    
    oneEventDigits.clear();	// empty vector 
    pTgr = (o2::hmpid::Trigger*)&trigger->at(i);
    
    for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
      // get all digits corresponding to one trigger, and assemble them 
      // in a vector
      digit = (o2::hmpid::Digit)digits->at(j);
      oneEventDigits.push_back(digit);
            
    }
   
    // get 2D-map for the event-number specified as input   
    if(i==nEvent) {
       
      for(unsigned int k = 0; k < oneEventDigits.size(); k++) {
      
        pDigEvt = (o2::hmpid::Digit*)&oneEventDigits.at(k);
                     
        o2::hmpid::Digit::pad2Absolute(pDigEvt->getPadID(), &module, &padChX, &padChY);
         
        hMap[module]->Fill(padChX, padChY, pDigEvt->getQ());                   
      }              
    }// end if 
  }
  // Divide canvases into 3x3     
  c1->Divide(3,3);
  c2->Divide(3,3); 
  c3->Divide(3,3);  	
  c4->Divide(3,3);  
	


 
 for(int iCh = 0; iCh<7; iCh++){
  
    c2->cd(pos[iCh]); 
    hCharge[iCh]->SetStats(false);
    hCharge[iCh]->Draw();
    
    c1->cd(pos[iCh]); 
    hMap[iCh]->SetStats(false);
    hMap[iCh]->Draw("Colz");
     
    c3->cd(pos[iCh]); 
    h_xCoord[iCh]->SetStats(false);
    h_xCoord[iCh]->Draw();
   
    c4->cd(pos[iCh]); 
    h_yCoord[iCh]->SetStats(false);
    h_yCoord[iCh]->Draw();

  }      
  c1->SaveAs("digitmap.eps");
  c2->SaveAs("digitcharge.eps");
  c3->SaveAs("digitX.eps");
  c4->SaveAs("digitY.eps");      
}
    
void strToFloatsSplit(std::string s,
                                            std::string delimiter, float* res,
                                            int maxElem)
{
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

//=======================
//
/*
void init()
{
  LOG(info) << "[HMPID Clusterization - init() ; mReadFile = ] "
            << mReadFile;
  //mSigmaCutPar = ic.options().get<std::string>("sigma-cut");

  if (mSigmaCutPar != "") {
    strToFloatsSplit(mSigmaCutPar, ",", mSigmaCut, 7);
  }

  mDigitsReceived, mClustersReceived = 0;

  mRec.reset(new o2::hmpid::Clusterer()); // ef: changed to smart-pointer

  mExTimer.start();

  // specify location and filename for output in case of writing to file
  if (mReadFile) {
    // Build the file name
    const auto filename = o2::utils::Str::concat_string(
      o2::utils::Str::rectifyDirectory(
        ic.options().get<std::string>("input-dir")),
      ic.options().get<std::string>("hmpid-digit-infile"));
    initFileIn(filename);
  }
	  // outputs
  std::vector<o2::hmpid::Cluster> clusters;
  std::vector<o2::hmpid::Trigger> clusterTriggers;
  LOG(info) << "[HMPID DClusterization - run() ] Enter ...";
  clusters.clear();
  clusterTriggers.clear();
	
	    LOG(info) << "[HMPID DClusterization - run() ] Entries  = " << mTree->GetEntries();

// check if more entries in tree
if (mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {


} else {
auto entry = mTree->GetReadEntry() + 1;
assert(entry < mTree->GetEntries());

mTree->GetEntry(entry);

// =============== create clusters =====================
for (const auto& trig : *mTriggersFromFilePtr) {
if (trig.getNumberOfObjects()) {
  gsl::span<const o2::hmpid::Digit> trigDigits{
    mDigitsFromFilePtr->data() + trig.getFirstEntry(),
    size_t(trig.getNumberOfObjects())};
  size_t clStart = clusters.size();
  mRec->Dig2Clu(trigDigits, clusters, mSigmaCut, true);
  clusterTriggers.emplace_back(trig.getIr(), clStart,
			       clusters.size() - clStart);
}
}

LOGP(info, "Received {} triggers with {} digits -> {} triggers with {} clusters",
   mTriggersFromFilePtr->size(), mDigitsFromFilePtr->size(), clusterTriggers.size(),
   clusters.size());
mDigitsReceived += mDigitsFromFilePtr->size();
} // <end else of num entries>
}*/

/*
void initFileIn(const std::string& filename)
{
  // Create the TFIle
  mFile = std::make_unique<TFile>(filename.c_str(), "OLD");
  assert(mFile && !mFile->IsZombie());

  mTree.reset((TTree*)mFile->Get("o2sim"));

  if (!mTree){
    mTree.reset((TTree*)mFile->Get("o2hmp"));
  }


  if (!mTree) {
    LOG(error)
      << "HMPID DigitToClusterSpec::init() : Did not find o2sim tree in "
      << filename.c_str();
    throw std::runtime_error(
      "HMPID DigitToClusterSpec::init() : Did not find "
      "o2sim file in digits tree");
  }

  if((mTree->GetBranchStatus("HMPDigit")) == 1) {
    mTree->SetBranchAddress("HMPDigit", &mDigitsFromFilePtr);
  } else if((mTree->GetBranchStatus("HMPIDDigits")) == 1) {
    mTree->SetBranchAddress("HMPIDDigits", &mDigitsFromFilePtr);
  } else {
    throw std::runtime_error(
      "HMPID DigitToClusterSpec::init() : Error in branches!");
  }


  mTree->SetBranchAddress("InteractionRecords", &mTriggersFromFilePtr);
  mTree->Print("toponly");
} */
