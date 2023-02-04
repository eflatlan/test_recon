#ifndef HMPIDTOOLS_H
#define HMPIDTOOLS_H

#include <vector>
#include <TCanvas.h>
#include <TH2F.h>
#include <TROOT.h> // gRoot


#include <string.h>

//namespace o2::hmpidTools
//{

using o2::hmpid::Cluster, o2::hmpid::Digit, o2::hmpid::Trigger,
    o2::hmpid::Clusterer;
using std::vector, std::cout, std::cin, std::endl;
using std::this_thread::sleep_for;


class HMPIDTools
{

  public:

    HMPIDTools() = default;
    ~HMPIDTools() = default;

    void setLinkStatus(int chamber, int link, bool status);            // chamber [0..6], link 0 = Left; 1 = Right
    void setSectorStatus(int chamber, int sector, bool status);        // chamber [0..6], sector [0..5]

    void setSectorStatus(TH2F* tHist, int chamber, int sector, bool status);        // chamber [0..6], sector [0..5]


    void setRadiatorStatus(int chamber, int radiator, bool status);    // chamber [0..6], radiator [0..2]

    void setLinkListStatus(vector<std::array<int, 2>>, bool status);      // Set A list of Links off
    void setSectorListStatus(vector<std::array<int, 2>>, bool status);    // Set A List of Sectors off
    void setRadiatorListStatus(vector<std::array<int, 2>>, bool status);  // Set A List of Radiators off
 
    void setPadChannel(int chamber, int xLow, int xHigh, int yLow, int yHigh, bool status);

    void setPadChannel(TH2F* tHist, int chamber, int xLow, int xHigh, int yLow, int yHigh, bool status);


    int getLinkStatus(int chamber, int link)
    {
      return sectorStatus[chamber][link];
    }

    int getSectorStatus(int chamber, int sector)
    {
      return sectorStatus[chamber][sector];
    }

    int getRadiatorStatus(int chamber, int radiator)
    {
      return radiatorStatus[chamber][radiator];
    }

    const std::size_t sectorHeight = 144/6;
    const std::size_t sectorWhidth = 160;

    const std::size_t linkWhidth = 160/2;

    void drawSectorStatus()
    {
      std::unique_ptr<TCanvas>  sectorStatusCanvas = make_unique<TCanvas>("Sector Status", "Sector Status", 1200, 1200);

      if(sectorStatusCanvas == nullptr){
        LOG(warn) << "Sector Status Canvas Was Nullptr" << endl;
        return;
        std::exit(0);
      }

      sectorStatusCanvas->Divide(3, 3);
      for(int iCh = 0; iCh < 7; iCh++){

        const char *sectorStatusString = Form("Map of Sector Status%i;Link; Sector", iCh);
        sectorStatusMap[iCh].reset(new TH2F(sectorStatusString, sectorStatusString, 160, 0, 159, 144, 0, 143));
        sectorStatusMap[iCh]->SetTitleSize(sectorStatusMap[iCh]->GetTitleSize("x") * 1.3, "xyz");
        //sectorStatusMap[iCh]->SetBins(160, 0, 159, 144, 0, 143);
        const auto& pos = posArr[iCh];
        const auto& pad = static_cast<TPad*>(sectorStatusCanvas->cd(pos));
        if(sectorStatusMap[iCh]!=nullptr){
          cout << "======" << endl;
          fillSectorStatus(sectorStatusMap[iCh].get(), iCh);
          sectorStatusMap[iCh]->SetStats(kFALSE);
          sectorStatusMap[iCh]->Draw("Colz");
        } else {
          cout << "sectorStatusMap nullptr" << endl;
        }
      }
      sectorStatusCanvas->Show();
      sectorStatusCanvas->SaveAs("SectorStatus.png");
    }

    void fillSectorStatus(TH2F* tHist, int chamber)
    { 
      
      for(int iSec = 0; iSec < 6; iSec++){
        const auto& sectorStatus = static_cast<double>(getSectorStatus(chamber, iSec));
        //sectorStatusMap[chamber]->Fill(0.0, static_cast<double>(iSec), sectorStatus); // x, y, value
        setSectorStatus(tHist, chamber, iSec, getSectorStatus(chamber, iSec));

        LOGP(info, "Chamber {} Sector {} Status {}", chamber, iSec, sectorStatus);
        if(getSectorStatus(chamber, iSec)) { 
          std::cout << "X" ;
        } else { std::cout << "O";}

      } 
    }


    void drawLinkStatus()
    {  
      std::unique_ptr<TCanvas>  linkStatusCanvas = make_unique<TCanvas>("Link Status", "Link Status", 1200, 1200);

      for(int iCh = 0; iCh < 7; iCh++)
      {
        const char *linkStatusString = Form("Map of Link Status%i;Link; Sector", iCh);
        linkStatusMap[iCh].reset(new TH2F(linkStatusString, linkStatusString, 1, 0, 1, 6, 0, 5));
        linkStatusMap[iCh]->SetTitleSize(linkStatusMap[iCh]->GetTitleSize("x") * 1.3, "xyz");
      }
    }

    void drawChannelStatus()
    {
      std::unique_ptr<TCanvas>  channelStatusCanvas = make_unique<TCanvas>("Pad Channel Status", "Pad Channel Status", 1200, 1200);

      for(int iCh = 0; iCh < 7; iCh++)
      {
        const char *channelStatusString = Form("Map of Pad-Channel Status%i; x [cm];y [cm]", iCh);
        channelStatusMap[iCh].reset(new TH2F(channelStatusString, channelStatusString, 160, 0, 159, 144, 0, 143));
        channelStatusMap[iCh]->SetTitleSize(channelStatusMap[iCh]->GetTitleSize("x") * 1.3, "xyz");
      }
    }

  private:
   bool channelStatus[7][160][144] = {{{true}}};
    //bool sectorStatus[7][160][144] = {{{true}}};
    //bool linkStatus[7][160][144] = {{{true}}};

    int sectorStatus[7][6] = {{1}};
    int linkStatus[7][2] = {{1}};
    int radiatorStatus[7][3] = {{1}};
    
    // position of pads in canvas according to HMPID-modules
    static constexpr int posArr[] = {9, 8, 6, 5, 4, 2, 1};

   //std::unique_ptr<TCanvas> 
   std::array<std::unique_ptr<TH2F>, 7> sectorStatusMap, linkStatusMap, channelStatusMap;


  ClassDef(HMPIDTools, 0);
}; // end class HMPIDTools


//} // end namespace o2::hmpidTools


#endif
