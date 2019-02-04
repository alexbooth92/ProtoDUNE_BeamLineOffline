#include "class_AcquisitionXBPF.h"

AcquisitionXBPF::AcquisitionXBPF()
{
  fFirstTriggerTimestamp = 2.e9;
  fLastTriggerTimestamp  = 0.;
  fCurrent = 508.4;

  return;
}

void AcquisitionXBPF::addAcqDetails(double cTimestamp, unsigned int cNEventRec)
{
  fTimestamp = cTimestamp;
  fNEventRec = cNEventRec;
  return;
}

void AcquisitionXBPF::clearAcquisition()
{
  fFirstTriggerTimestamp = 2.e9;
  fLastTriggerTimestamp  = 0.;
  fCurrent = 508.4;
  fTimestamp = 0.;
  fNEventRec = 0;
  fVecDataRec.clear();
  fVecDataRecHR.clear();

  fIsEmpty   = true;

  return;
}

bool AcquisitionXBPF::isEmpty()
{
  return fIsEmpty;
}

AcquisitionXBPF::EventRecordHR AcquisitionXBPF::decodeFibres(EventRecordRaw const &cEventRecord)
{
  //LSB IS 8NS PER TICK.
  long long conversion_TickToNs = 8;
  long long conversion_sToNs    = 1000000000;
  long long TTS_LSB = cEventRecord.fTTS_LSB*conversion_TickToNs;
  long long TTS_MSB = cEventRecord.fTTS_MSB*conversion_sToNs;
  long long ETS_LSB = cEventRecord.fETS_LSB*conversion_TickToNs;
  long long ETS_MSB = cEventRecord.fETS_MSB*conversion_sToNs;

  long long triggerTimestamp = TTS_MSB+TTS_LSB;
  long long eventTimestamp   = ETS_MSB+ETS_LSB;

  std::vector<unsigned int> fibrelist;
  std::array<bool, 192> arr_Fibres = {0};
  unsigned int nFibresHit     = 0;
  unsigned int lowestChannel  = 1000;
  unsigned int highestChannel = 1000;
  unsigned int count          = 0;

  for(unsigned int i = 0; i < 6; i++)
  {
    if(cEventRecord.fArrFibres[i]!=0)
    {
      std::bitset<32> currentBitSet{cEventRecord.fArrFibres[i]};
      for(unsigned int j = 0; j < 32; j++)
      {
        arr_Fibres[j+count ] = currentBitSet[j];
        if(currentBitSet[j]==1)
        {
          fibrelist.push_back(j+count);
          nFibresHit++;
          if(lowestChannel==1000)
          {
            lowestChannel = j + count;
          }
          else
          {
            highestChannel = j + count;
          }
        }
      }
    }
    count+=32;
  }
  if(triggerTimestamp/1.e9>fLastTriggerTimestamp)
  {
    fLastTriggerTimestamp = triggerTimestamp/1.e9;
  }
  if(triggerTimestamp/1.e9<fFirstTriggerTimestamp && triggerTimestamp!=0)
  {
    fFirstTriggerTimestamp = triggerTimestamp/1.e9;
  }

  return {triggerTimestamp, eventTimestamp, arr_Fibres, fibrelist, nFibresHit, (nFibresHit > 1 ? (highestChannel - lowestChannel) : 0)}; 
}

AcquisitionXBPF::EventRecordHR AcquisitionXBPF::decodeFibres(VecBoolVar const &cEventRecord)
{
  //LSB IS 8NS PER TICK.
  long long triggerTimestamp = cEventRecord.fTimeS*1000000000 + cEventRecord.fTimeNS;
  long long eventTimestamp   = triggerTimestamp;

  std::vector<unsigned int> fibrelist;
  std::array<bool, 192> arr_Fibres = {0};
  unsigned int nFibresHit     = 0;
  unsigned int lowestChannel  = 1000;
  unsigned int highestChannel = 1000;

  for(unsigned int i = 0; i < cEventRecord.fVar.size(); i++)
  {
    arr_Fibres[i] = cEventRecord.fVar[i];
    if(cEventRecord.fVar[i]!=0)
    {
      fibrelist.push_back(i);
      nFibresHit++;
      if(lowestChannel==1000)
      {
        lowestChannel = i;
      }
      else
      {
        highestChannel = i;
      }
    }
  }

  if(triggerTimestamp/1.e9>fLastTriggerTimestamp)
  {
    fLastTriggerTimestamp = triggerTimestamp/1.e9;
  }
  if(triggerTimestamp/1.e9<fFirstTriggerTimestamp && triggerTimestamp!=0)
  {
    fFirstTriggerTimestamp = triggerTimestamp/1.e9;
  }

  return {triggerTimestamp, eventTimestamp, arr_Fibres, fibrelist, nFibresHit, (nFibresHit > 1 ? (highestChannel - lowestChannel) : 0)}; 
}

double AcquisitionXBPF::getFirstTriggerTimestamp()
{
  return fFirstTriggerTimestamp;
}

double AcquisitionXBPF::getLastTriggerTimestamp()
{
  return fLastTriggerTimestamp;
}

double AcquisitionXBPF::getTimestamp()
{
  return fTimestamp;
}

double AcquisitionXBPF::getCurrent()
{
  return fCurrent;
}

unsigned int AcquisitionXBPF::getNEventRec()
{
  return fNEventRec;
}

std::vector<AcquisitionXBPF::EventRecordRaw> AcquisitionXBPF::getDataRaw()
{
  return fVecDataRec;
}

std::vector<AcquisitionXBPF::EventRecordHR> AcquisitionXBPF::getDataHR()
{
  return fVecDataRecHR;
}

void AcquisitionXBPF::addEventRecord(EventRecordRaw cEventRecord)
{
  fVecDataRec  .push_back(cEventRecord);
  fVecDataRecHR.push_back(decodeFibres(cEventRecord));
  return;
}

void AcquisitionXBPF::addEventRecord(VecBoolVar cVecBoolVar)
{
  fVecDataRecHR.push_back(decodeFibres(cVecBoolVar));
  return;
}

void AcquisitionXBPF::printRaw()
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
    std::cout << "TIMESTAMP: " << fTimestamp << ", NUMBER OF EVENTS RECORDED: " << fNEventRec  << ", NUMBER OF NONZERO EVENTS: " << fNNonZeroEvents << std::endl;
    std::cout << std::endl;
    for(unsigned int i = 0; i < fVecDataRec.size(); i++)
    {
      std::cout << "EVENT " << i   << ", TRIGTS_LSB: "     << fVecDataRec.at(i).fTTS_LSB << ", TRIGTS_MSB: " << fVecDataRec.at(i).fTTS_MSB
                << ", EVTTS_LSB: " << fVecDataRec.at(i).fETS_LSB << ", EVTTS_MSB: " << fVecDataRec.at(i).fETS_MSB  << ". ";  
      std::cout << "DATA: ";
      for(unsigned int j = 0; j < 6; j++)
      {
        std::cout << fVecDataRec.at(i).fArrFibres[j] << " ";
      }
      std::cout << std::endl; 
    }
    std::cout << "***********************************************************************************" << std::endl;
  }
  return;
}

void AcquisitionXBPF::printHR()
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
    std::cout << "TIMESTAMP: " << fTimestamp << ", NUMBER OF EVENTS RECORDED: " << fNEventRec  << ", NUMBER OF NONZERO EVENTS: " << fNNonZeroEvents << std::endl;
    std::cout << std::endl;
    for(unsigned int i = 0; i < fVecDataRecHR.size(); i++)
    {
      std::cout << "EVENT " << i   << ", TRIGTS: " << fVecDataRecHR.at(i).fTriggerTimestamp << ", EVTTS: " << fVecDataRecHR.at(i).fEventTimestamp << "." << std::endl;
      std::cout << "DATA: ";
      for(unsigned int j = 0; j < 192; j++)
      {
        if(j % 32 == 0)
        {
          std::cout << "\nCHANNELS: [" << j << "-" << j+31 << "]" << std::endl; 
        }
        std::cout << fVecDataRecHR.at(i).fArrFibres[j] << " ";
      }
      std::cout << std::endl; 
    }
    std::cout << "***********************************************************************************" << std::endl;
  }
  return;
}

void AcquisitionXBPF::setEmpty(bool cIsEmpty)
{
  fIsEmpty = cIsEmpty;
  return;
}

void AcquisitionXBPF::setCurrent(double cCurrent)
{
  fCurrent = cCurrent;
  return;
}
