#include "class_AcquisitionXBPF.h"
#include "class_Detector.h"

Detector::Detector()
{
  return;
}

Detector::Detector(std::string cDetName, std::string cType, std::string cDetLabel)
{
  fDetName  = cDetName;
  fType     = cType;
  fDetLabel = cDetLabel;
  return;
}

unsigned int Detector::getNAcquisitions()
{
  if(fType=="XBPF")
  {
    return fVecAcquisitionXBPF.size();
  }
  //else if(fType=="XBTF")
  else 
  {
    return fVecAcquisitionXBTF.size();
  }
}

double Detector::getFirstTimestamp()
{
  return fFirstTimestamp;
}

double Detector::getLastTimestamp()
{
  return fLastTimestamp;
}

double Detector::getFirstTriggerTimestamp()
{
  return fDetFirstTriggerTimestamp;
}

double Detector::getLastTriggerTimestamp()
{
  return fDetLastTriggerTimestamp;
}

double Detector::getTimeElapsed()
{
  return (fDetLastTriggerTimestamp-fDetFirstTriggerTimestamp);
}

std::string Detector::getDetName()
{
  return fDetName;
}

std::string Detector::getDetLabel()
{
  return fDetLabel;
}

std::string Detector::getType()
{
  return fType;
}

void Detector::getAcquisitions(std::vector<AcquisitionXBPF> &acq)
{
  acq = fVecAcquisitionXBPF;
  return;
}

void Detector::getAcquisitions(std::vector<AcquisitionXBTF> &acq)
{
  acq = fVecAcquisitionXBTF;
  return;
}

void Detector::addAcquisition(AcquisitionXBPF cAcq)
{
  if(fVecAcquisitionXBPF.size()==0)
  {
    fFirstTimestamp = cAcq.getTimestamp();
    fLastTimestamp  = cAcq.getTimestamp();
    fDetFirstTriggerTimestamp = cAcq.getFirstTriggerTimestamp();
    fDetLastTriggerTimestamp  = cAcq.getLastTriggerTimestamp();
  }

  fVecAcquisitionXBPF.push_back(cAcq);

  if(cAcq.getTimestamp()>fLastTimestamp)
  {
    fLastTimestamp = cAcq.getTimestamp();   
  }
  if(cAcq.getLastTriggerTimestamp()>fDetLastTriggerTimestamp)
  {
    fDetLastTriggerTimestamp = cAcq.getLastTriggerTimestamp();
  }

  return;
}

void Detector::addAcquisition(AcquisitionXBTF cAcq)
{
  if(fVecAcquisitionXBTF.size()==0)
  {
    fFirstTimestamp = cAcq.getTimestamp();
    fLastTimestamp  = cAcq.getTimestamp();
    fDetFirstTriggerTimestamp = cAcq.getFirstTriggerTimestamp();
    fDetLastTriggerTimestamp  = cAcq.getLastTriggerTimestamp();
  }

  fVecAcquisitionXBTF.push_back(cAcq);

  if(cAcq.getTimestamp()>fLastTimestamp)
  {
    fLastTimestamp = cAcq.getTimestamp();   
  }
  if(cAcq.getLastTriggerTimestamp()>fDetLastTriggerTimestamp)
  {
    fDetLastTriggerTimestamp = cAcq.getLastTriggerTimestamp();
  }

  return;
}

void Detector::printDescription()
{
  std::cout << "****************************************************************************************" << std::endl;
  std::cout << fDetName << ": OF TYPE " << fType << ", WITH LABEL " << fDetLabel                          << std::endl;
  std::cout << "****************************************************************************************" << std::endl;

  return;
}

void Detector::printRaw()
{
  std::cout << "-----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "PRINTING DETECTOR " << fDetName << std::endl;
  for(unsigned int i = 0; i < fVecAcquisitionXBPF.size(); i++)
  {
    fVecAcquisitionXBPF.at(i).printRaw();
  }
  std::cout << "-----------------------------------------------------------------------------------------------------------------------------------" << std::endl;

  return;
}

void Detector::printHR()
{
  std::cout << "-----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "PRINTING DETECTOR " << fDetName << std::endl;
  for(unsigned int i = 0; i < fVecAcquisitionXBPF.size(); i++)
  {
    fVecAcquisitionXBPF.at(i).printHR();
  }
  std::cout << "-----------------------------------------------------------------------------------------------------------------------------------" << std::endl;

  return;
}
