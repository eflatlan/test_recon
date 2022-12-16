#if !defined(__CLING__) || defined(__ROOTCLING__)
//#if !defined(__CINT__) || defined(__MAKECINT__)
#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TLine.h>
#include <TList.h>
#include <TROOT.h> // gRoot
#include <TRandom.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <vector>

// C++ header files and libraries
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#endif

void SaveFolder(string inpFile);
void changeFont();

void readClusters(int nEvents) {
  TH1F *hCharge[7], *hMipCharge[7], *hSize[7];
  TH2F *hMap[7];

  for (int i = 0; i < 7; i++) {
    hMap[i] =
        new TH2F(Form("Clusters Map chamber%i", i),
                 Form("Cluster Map chamber%i", i), 160, 0, 159, 144, 0, 143);
    hMap[i]->SetXTitle("X (cm)");
    hMap[i]->SetYTitle("Y (cm)");

    hCharge[i] =
        new TH1F(Form("Clusters Charge chamber%i", i),
                 Form("Cluster Charge chamber%i", i), 2000, 100., 2100.);
    hCharge[i]->SetXTitle("Charge (ADC channel)");
    hCharge[i]->SetYTitle("Entries");

    hMipCharge[i] =
        new TH1F(Form("Mip Clusters Charge chamber%i", i),
                 Form("Mip Cluster Charge chamber%i", i), 50, 200., 2200.);
    hMipCharge[i]->SetXTitle("Charge (ADC channel)");
    hMipCharge[i]->SetYTitle("Entries/40 ADC");
    hMipCharge[i]->SetLineColor(kBlack);

    hSize[i] = new TH1F(Form("Cluster Size chmaber%i", i),
                        Form("Cluster Size chamber%i", i), 20, 0., 20.);
    hSize[i]->SetXTitle("Cluster size");
    hSize[i]->SetYTitle("Entries");
  }

  // changeFont();			     // specify folder to save files in
  // SaveFolder("clusterChambers");   // apply custom canvas figure options

  auto folderName = (gSystem->GetWorkingDirectory());
  auto gwd = gSystem->GetWorkingDirectory();

  auto runNumber = (gwd.substr(gwd.length() - 9, gwd.length()));
  /*
  Printf("RunNumber %s", runNumber);
  Printf("RunNumber %s", runNumber.c_str());

  Printf("gwd %s",gwd.c_str());
  Printf("pwd %s",gSystem->pwd());
  Printf("homedir %s",gSystem->HomeDirectory()); */
  Printf("Empty clusters");

  TCanvas *c1 =
      new TCanvas("c1", ("Cluster-Map " + runNumber).c_str(), 2000, 1200);
  TCanvas *c2 =
      new TCanvas("c2", ("Cluster-Charge " + runNumber).c_str(), 2000, 1200);
  TCanvas *c3 = new TCanvas("c3", ("MIP Cluster-Charge " + runNumber).c_str(),
                            2000, 1200);
  TCanvas *c4 =
      new TCanvas("c4", ("Cluster-Size " + runNumber).c_str(), 2000, 1200);
  /*
  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  TCanvas *c3 = new TCanvas("c3","c3",1000,800);
  TCanvas *c4 = new TCanvas("c4","c4",1000,800); */

  Int_t pos[] = {9, 8, 6, 5, 4, 2, 1};

  Int_t nTotTriggers = 0;

  for (int k = 0; k < nEvents; k++) {

    Printf("hmpclus%02i.root", k + 1);
    std::unique_ptr<TFile> fileClusters{
        TFile::Open(Form("hmpclus%02i.root", k + 1))};
    std::unique_ptr<TTree> treeClusters;
    treeClusters.reset((TTree *)fileClusters->Get("o2sim"));
    o2::hmpid::Trigger *pTgr;
    o2::hmpid::Cluster *pClu;
    o2::hmpid::Cluster *pCluEvt;
    o2::hmpid::Cluster cluster;

    std::vector<o2::hmpid::Cluster> *clusters = nullptr;
    std::vector<o2::hmpid::Cluster> oneEventClusters;
    std::vector<o2::hmpid::Trigger> *trigger = nullptr;

    if (!treeClusters) {
      Printf("Empty clusters");
      continue;
    } else {
      treeClusters->SetBranchAddress("HMPIDClusters", &clusters);
      treeClusters->SetBranchAddress("InteractionRecords", &trigger);
      Printf("Got clusters");
    }

    const int treeClusSize = treeClusters->GetEntries();
    Printf("tree entries = %i", treeClusSize);

    treeClusters->GetEntry(0);

    const int clusterSize = clusters->size();
    Printf("clusters size = %i", clusterSize);

    int module = 0;

    for (int j = 0; j < clusterSize; j++) {

      pClu = (o2::hmpid::Cluster *)&clusters->at(j);

      module = pClu->ch();

      hCharge[module]->Fill(pClu->q());

      if (pClu->size() >= 3 && pClu->size() <= 7)
        hMipCharge[module]->Fill(pClu->q());

      hSize[module]->Fill(pClu->size());
    }

    nTotTriggers += trigger->size();
    const int triggerSize = trigger->size();
    Printf("trigger size from clusters = %i", triggerSize);

    for (int i = 0; i < triggerSize; i++) {
      oneEventClusters.clear();
      pTgr = static_cast<o2::hmpid::Trigger *>(&trigger->at(i));

      const int pTrgrFirst = pTgr->getFirstEntry();
      const int pTrgrLast = pTgr->getLastEntry();
      // for(int j = pTrgrFirst; j <= pTrgrLast; j++) {
      for (int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
        // cluster = static_cast<o2::hmpid::Cluster>(clusters->at(j)); //
        // ->at(j) to [j]
        cluster =
            static_cast<o2::hmpid::Cluster>(clusters->at(j)); // ->at(j) to [j]
        oneEventClusters.push_back(cluster);
      }
    } // exit for2s

  } // exit for3

  c1->Divide(3, 3);
  c2->Divide(3, 3);
  c3->Divide(3, 3);
  c4->Divide(3, 3);

  for (int iCh = 0; iCh < 7; iCh++) {
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
  // Printf("Number of triggers = %i", nTotTriggers);
  // Printf(

  c1->SaveAs(("clusterMap_" + runNumber + "_.eps").c_str());
  c2->SaveAs(("clusterCharge_" + runNumber + "_.eps").c_str());
  c3->SaveAs(("mipClustesCharge_" + runNumber + "_.eps").c_str());
  c4->SaveAs(("clusterSize_" + runNumber + "_.eps").c_str());

  c1->SaveAs(("clusterMap_" + runNumber + "_.png").c_str());
  c2->SaveAs(("clusterCharge_" + runNumber + "_.png").c_str());
  c3->SaveAs(("mipClustesCharge_" + runNumber + "_.png").c_str());
  c4->SaveAs(("clusterSize_" + runNumber + "_.png").c_str());
}

// apply custom-made canvas options
void changeFont() {
  TStyle *canvasStyle = new TStyle("canvasStyle", "Canvas Root Styles");
  canvasStyle->SetPalette(1, 0);
  canvasStyle->SetTitleSize(0.085, "xy"); // size of axis title font
  canvasStyle->SetTitleFont(22, "xz");    // font option
  canvasStyle->SetTitleFontSize(0.1);     // size of canvas-title
  canvasStyle->SetTitleOffset(.825, "y"); //  y-axis title-offset from axis
  canvasStyle->SetTitleOffset(1, "z");
  canvasStyle->SetTitleOffset(.95, "x");
  canvasStyle->SetTitleX(.25); // set canvas-title position
  // labels
  canvasStyle->SetLabelOffset(0.005, "y");
  canvasStyle->SetLabelFont(22, "xyz");
  canvasStyle->SetLabelSize(0.085, "xyz"); // size of axis value font
  // canvas
  canvasStyle->SetCanvasColor(0);
  canvasStyle->SetCanvasBorderMode(0);
  canvasStyle->SetCanvasBorderSize(0);
  // set margins
  canvasStyle->SetPadBottomMargin(0.18);
  canvasStyle->SetPadTopMargin(0.05);
  canvasStyle->SetPadLeftMargin(0.13);
  canvasStyle->SetPadRightMargin(0.02);
  gROOT->SetStyle("canvasStyle");
}

void SaveFolder(string inpFile) {
  char *createdFolder;
  auto time = std::chrono::system_clock::now();
  std::time_t time_t = std::chrono::system_clock::to_time_t(time);
  auto c_time = std::ctime(&time_t);

  // Allocate memory for char* createdFolder
  int pathLen = strlen(c_time);
  pathLen = strlen(gSystem->pwd()) + strlen(inpFile.c_str());
  int numOfSigns = 2; // allocate for - and /
  createdFolder = (char *)calloc(pathLen + numOfSigns + 1, sizeof(char));

  // Copy Base-directory into empty createdFolder
  strcpy(createdFolder, gSystem->pwd());
  // add hyphen to the Folder-name
  strcat(createdFolder, "/");
  // Add name of read file to Folder-name
  strcat(createdFolder, inpFile.c_str());
  // add - to the Folder-name
  strcat(createdFolder, "-");
  // append current time to the Folder-name
  strcat(createdFolder, c_time);

  // Make new Directory from the newly
  // created char* createdFolder
  gSystem->MakeDirectory(createdFolder);

  // Move to directory such that new
  // canvas-files will be saved here
  gSystem->ChangeDirectory(createdFolder);
}
