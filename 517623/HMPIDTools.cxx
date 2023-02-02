#include "HMPIDTools.h"

ClassImp(HMPIDTools);
//namespace o2::hmpidtools
//{


// chamber [0..6], link 0 = Left; 1 = Right
void HMPIDTools::setLinkStatus(int chamber, int link, bool status)
{
  mempcopy(linkStatus[chamber, link], status);
  //linkStatus[chamber, link] = status;
  const int xMin = sectorWhidth*link;
  const int xMax = sectorWhidth*(link+1);
  const int yMin = 0;
  const int yMax = sectorHeight;
  setPadChannel(chamber, xMin, xMax, yMin, yMax, status);
}        

// chamber [0..6], sector [0..5]
void HMPIDTools::setSectorStatus(int chamber, int sector, bool status)     
{
  sectorStatus[chamber, sector] = status;
  const int xMin = 0;
  const int xMax = sectorWhidth;
  const int yMin = sectorHeight*sector;
  const int yMax = sectorHeight*(sector+1);
  setPadChannel(chamber, xMin, xMax, yMin, yMax, status);
}   


// chamber [0..6], radiator [0..2]
//void HMPIDTools::setRadiatorStatus(int chamber, int radiator);    

// Set A list of Links off
void HMPIDTools::setLinkListStatus(vector<std::array<int, 2>> linkList, bool status)
{
  for(const auto& linkChamber : linkList)
  {
    const auto& chamber = linkChamber[0];
    const auto& link = linkChamber[0];
    setLinkStatus(chamber, link, status);
  }
}      

 // Set A List of Sectors off
void HMPIDTools::setSectorListStatus(vector<std::array<int, 2>> sectorList, bool status)
{
  for(const auto& sectorChamber : sectorList)
  {
    const auto& chamber = sectorChamber[0];
    const auto& sector = sectorChamber[0];
    setSectorStatus(chamber, sector, status);
  }
}

// Set A List of Radiators off
/*void HMPIDTools::setRadiatorListStatus(vector<std::array<int, 2>> radiatorList, bool status)
{
   
} */ 


// Set Pad Channels Off
void HMPIDTools::setPadChannel(int chamber, int xLow, int xHigh, int yLow, int yHigh, bool status)
{
  for (int chamber = 0; chamber < 7; chamber++) {
    for (int x = 0; x < 160; x++) {
      for (int y = 0; y < 144; y++) {
        channelStatus[chamber][x][y] = status;
      }
    }
  }

}


//} // end namespace
