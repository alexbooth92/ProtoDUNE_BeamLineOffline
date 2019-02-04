#include "class_AcquisitionXBTF.h"

AcquisitionXBTF::AcquisitionXBTF()
{
  fFirstTriggerTimestamp = 2.e9;
  fLastTriggerTimestamp   = 0.;
  return;
}

void AcquisitionXBTF::addAcqDetails(double cTimestamp, unsigned int cTimestampCount)
{
  fTimestamp      = cTimestamp;
  fTimestampCount = cTimestampCount;
  return;
}

void AcquisitionXBTF::clearAcquisition()
{
  fFirstTriggerTimestamp = 2.e9;
  fLastTriggerTimestamp   = 0.;
  fTimestamp = 0.;
  fTimestampCount = 0.;
  fVecDataRec.clear();
  fVecDataRecHR.clear();

  fIsEmpty = true;
  return;
}

bool AcquisitionXBTF::isEmpty()
{
  return fIsEmpty;
}

AcquisitionXBTF::EventRecordHR AcquisitionXBTF::decodeRaw(EventRecordRaw const &cEventRecord)
{
  unsigned long seconds         = cEventRecord.fSeconds;
  unsigned long long NSAccurate = cEventRecord.fSeconds*1000000000 + cEventRecord.fCoarse*8; 
  double subSeconds             = cEventRecord.fCoarse*8. + cEventRecord.fFrac/512.;
  double seconds_FullTime       = cEventRecord.fSeconds + (cEventRecord.fCoarse*8 + cEventRecord.fFrac/512.)/1.e9;

  if(seconds_FullTime<fFirstTriggerTimestamp)
  {
    fFirstTriggerTimestamp = seconds_FullTime;
  }
  if(seconds_FullTime>fLastTriggerTimestamp)
  {
    fLastTriggerTimestamp = seconds_FullTime;
  }

  return {seconds, subSeconds, seconds_FullTime, NSAccurate, cEventRecord.fFrac/512.}; 
}

double AcquisitionXBTF::getTimestamp()
{
  return fTimestamp;
}

double AcquisitionXBTF::getFirstTriggerTimestamp()
{
  return fFirstTriggerTimestamp;
}

double AcquisitionXBTF::getLastTriggerTimestamp()
{
  return fLastTriggerTimestamp;
}

unsigned int AcquisitionXBTF::getTimestampCount()
{
  //return fTimestampCount;
  return fVecDataRecHR.size();
}

std::vector<AcquisitionXBTF::EventRecordRaw> AcquisitionXBTF::getDataRaw()
{
  return fVecDataRec;
}

std::vector<AcquisitionXBTF::EventRecordHR> AcquisitionXBTF::getDataHR()
{
  return fVecDataRecHR;
}

void AcquisitionXBTF::addEventRecord(EventRecordRaw cEventRecord)
{
  fIsEmpty = false;
  fVecDataRec  .push_back(cEventRecord);
  fVecDataRecHR.push_back(decodeRaw(cEventRecord));
  return;
}

void AcquisitionXBTF::printRaw()
{
  std::cout.precision(14);
  if(fIsEmpty)
  {
    std::cout << "***********************************************************************************" << std::endl;
    std::cout << "*  TIMESTAMP: " << fTimestamp << " IS EMPTY!" << std::endl;
    std::cout << "***********************************************************************************" << std::endl;
  }
  else
  {
    std::cout << "***********************************************************************************" << std::endl;
    std::cout << "TIMESTAMP: " << fTimestamp << ", NUMBER OF EVENTS RECORDED: " << fTimestampCount          << std::endl;
    std::cout << std::endl;
    for(unsigned int i = 0; i < fVecDataRec.size(); i++)
    {
      std::cout << "EVENT " << i   << ", FRAC: " << fVecDataRec.at(i).fFrac << ", COARSE: " << fVecDataRec.at(i).fCoarse
                << ", SECONDS: " << fVecDataRec.at(i).fSeconds << "." << std::endl;  
    }
    std::cout << "***********************************************************************************" << std::endl;
  }
  return;
}

void AcquisitionXBTF::printHR()
{
  std::cout.precision(20);
  if(fIsEmpty)
  {
    std::cout << "***********************************************************************************" << std::endl;
    std::cout << "*  TIMESTAMP: " << fTimestamp << " IS EMPTY!" << std::endl;
    std::cout << "***********************************************************************************" << std::endl;
  }
  else
  {
    std::cout << "***********************************************************************************" << std::endl;
    std::cout << "TIMESTAMP: " << fTimestamp << ", NUMBER OF EVENTS RECORDED: " << fTimestampCount          << std::endl;
    std::cout << std::endl;
    for(unsigned int i = 0; i < fVecDataRecHR.size(); i++)
    {
      std::cout << "EVENT " << i   << ", SECONDS: " << fVecDataRecHR.at(i).fSeconds << ", SUBSECONDS: " << fVecDataRecHR.at(i).fSubSeconds
                << ", FULL TIME: " << fVecDataRecHR.at(i).fSeconds_FullTime << ", NS, ACCURATE TO NS: " << fVecDataRecHR.at(i).fNSAccurate << std::endl;  
    }
    std::cout << "***********************************************************************************" << std::endl;
  }
  return;
}

void AcquisitionXBTF::setEmpty(bool cIsEmpty)
{
  fIsEmpty = cIsEmpty;
  return;
}

void AcquisitionXBTF::sortByNSAccurate()
{
  std::sort(fVecDataRecHR.begin(), fVecDataRecHR.end(), [](const EventRecordHR &left, const EventRecordHR &right){return ((left.fNSAccurate)<(right.fNSAccurate));}); 
  return;
}
