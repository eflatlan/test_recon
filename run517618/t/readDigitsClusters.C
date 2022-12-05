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
#endif

void readDigits(char* filename, int nEvent)
{
  TFile *fileDigits = TFile::Open(filename);
  
  TTree *treeDigits = (TTree*)fileDigits->Get("o2sim");
  
  TH1F *hCharge[7];
  TH2F *hMap[7];
      
  for(int i=0; i<7; i++) {
     hMap[i] = new TH2F(Form("Digits Map chmaber%i",i),Form("Digits Map chmaber%i",i), 160, 0, 159, 144, 0, 143);
     hMap[i]->SetXTitle("pad X");
     hMap[i]->SetYTitle("pad Y"); 
     
     hCharge[i] = new TH1F(Form("Digits Charge chamber%i",i),Form("Digits Charge chamber%i",i),2000, 100., 2100.);
     hCharge[i]->SetXTitle("Charge (ADC channel)");               
     hCharge[i]->SetYTitle("Entries");               
   }
   
  TCanvas *c1 = new TCanvas("c1","c1",1000,800); 
  TCanvas *c2 = new TCanvas("c2","c2",1000,800); 
   
  Int_t pos[]= {9,8,6,5,4,2,1};
  
  o2::hmpid::Trigger *pTgr;
  o2::hmpid::Digit *pDig;
  o2::hmpid::Digit *pDigEvt;
  o2::hmpid::Digit digit;
  
  std::vector<o2::hmpid::Digit> *digits = nullptr;
  std::vector<o2::hmpid::Digit> oneEventDigits;
  std::vector<o2::hmpid::Trigger> *trigger = nullptr;
  
  treeDigits->SetBranchAddress("HMPDigit",&digits);
  treeDigits->SetBranchAddress("InteractionRecords",&trigger);
  

  const int treeDigitsSize = treeDigits->GetEntries();
  Printf("tree entries = %i", treeDigitsSize);
  
  treeDigits->GetEntry(0);


  const int digitsSize = digits->size();
  Printf("digit size = %i", digitsSize);
  
  int padChX = 0, padChY = 0, module = 0;
    
  for(int j = 0; j < digitsSize; j++) {
      
    pDig = (o2::hmpid::Digit*)&digits->at(j);

    o2::hmpid::Digit::pad2Absolute(pDig->getPadID(), &module, &padChX, &padChY);    
                         
    hCharge[module]->Fill(pDig->getQ());
   }
       

  const int triggerSize = trigger->size();
  Printf("trigger size from digits = %i", triggerSize); 
   
  for(int i = 0; i < triggerSize; i++) {
    
    oneEventDigits.clear();
    //pTgr = (o2::hmpid::Trigger*)&trigger->at(i);
    pTgr = static_cast<o2::hmpid::Trigger*>(&trigger->at(i));
    for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
      
      digit = (o2::hmpid::Digit)digits->at(j);
      oneEventDigits.push_back(digit);
      
     }
      
    if(i==nEvent) {
       const int oneEventDigSize = oneEventDigits.size();
       for(int k = 0; k < oneEventDigSize; k++) {
      
       //pDigEvt = (o2::hmpid::Digit*)&oneEventDigits.at(k);
       pDigEvt = static_cast<o2::hmpid::Digit*>(&(oneEventDigits[k]));
       o2::hmpid::Digit::pad2Absolute(pDigEvt->getPadID(), &module, &padChX, &padChY);
         
       hMap[module]->Fill(padChX, padChY, pDigEvt->getQ());
                    
      }              
     }
   }
       
   c1->Divide(3,3);
   c2->Divide(3,3);

   for(int iCh = 0; iCh<7; iCh++){
     c1->cd(pos[iCh]); 
     hMap[iCh]->Draw("Colz");
     
     c2->cd(pos[iCh]); 
     hCharge[iCh]->Draw();
   }      
}
//************************************************************************************************************************************************************************************************************************************************************************************
void readClusters(char* filename, int nEvent)
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
   
  TCanvas *c1 = new TCanvas("c1","c1",1000,800); 
  TCanvas *c2 = new TCanvas("c2","c2",1000,800); 
  TCanvas *c3 = new TCanvas("c3","c3",1000,800); 
  TCanvas *c4 = new TCanvas("c4","c4",1000,800);
   
  Int_t pos[] = {9,8,6,5,4,2,1};
  
  Int_t nTotTriggers = 0;
  

  std::unique_ptr<TTree> treeClusters;
  for(int k = 0; k<50; k++) {
    
    Printf("hmpclus%02i.root", k+1);
      
    //TFile *fileClusters = TFile::Open(Form("%shmpclus%03i.root", k+1));
    std::unique_ptr<TFile> fileClusters{ TFile::Open(Form("hmpclus%02i.root", k+1))};


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
     }
     
    else       
      {
       treeClusters->SetBranchAddress("HMPIDClusters",&clusters);
       treeClusters->SetBranchAddress("InteractionRecords",&trigger);
      }
    

    const int treeClusSize = treeClusters->GetEntries();
    Printf("tree entries = %i", treeClusSize);
    
    treeClusters->GetEntry(0);


    const int clusterSize = clusters->size();
    Printf("clusters size = %i", clusterSize);
  
    int module = 0;
    
    for(int j = 0; j < clusterSize; j++) {
      
      pClu = (o2::hmpid::Cluster*)&clusters->at(j);

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
     pTgr = (o2::hmpid::Trigger*)&trigger->at(i);
     const int pTrgrFirst = pTgr->getFirstEntry();
     const int pTrgrLast = pTgr->getLastEntry();
     // for(int j = pTrgrFirst; j <= pTrgrLast; j++) {
     for(int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) { 
      //cluster = static_cast<o2::hmpid::Cluster>(clusters->at(j)); // ->at(j) to [j]
      cluster = (o2::hmpid::Cluster) clusters->at(j); // ->at(j) to [j]
      oneEventClusters.push_back(cluster);

     }

       /*if(i==nEvent) {       
       for(unsigned int k = 0; k < oneEventClusters.size(); k++) {
       pCluEvt = (o2::hmpid::Cluster*)&oneEventClusters.at(k);        
       module = pCluEvt->ch();
       hMap[module]->Fill(pCluEvt->x(), pCluEvt->y());}*/              
     } Printf("Exit forLoop ");
   }

  Printf("Exit nestForLoop");

 
  Printf("Number of triggers222 = %i", nTotTriggers);       
   c1->Divide(3,3);
   c2->Divide(3,3);
   c3->Divide(3,3);
   c4->Divide(3,3);

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
   
  Printf("Number of triggers = %i", nTotTriggers);       
}    
   
    
