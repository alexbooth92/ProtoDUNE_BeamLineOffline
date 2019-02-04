#include "class_DetectorCoincidences.h"
#include "TInterpreter.h"

double Pi       = 3.14159265;
double radToDeg = 180./Pi;

std::vector<TFCoincidenceRecord::TFCoincidence> TFCoincidenceRecord::getTFCoincidencesDegenerate()
{
  if(fVecTFCoincidesUnique.size()==0)
  {
    getTFCoincidencesUnique();
  }

  return fVecTFCoincidesDegenerate;
}

std::vector<TFCoincidenceRecord::TFCoincidence> TFCoincidenceRecord::getTFCoincidencesUnique()
{
  if(fVecTFCoincidesUnique.size()==0)
  {
    for(unsigned int i = 0; i < fVecTFCoincidesDegenerate.size(); i++)
    {
      std::pair<unsigned int, unsigned int> multiplicity = getMultiplicity(fVecTFCoincidesDegenerate[i]);
      fVecTFCoincidesDegenerate[i].setMultiplicity(multiplicity);
      if(multiplicity.first==1 && multiplicity.second==1)
      {
        fVecTFCoincidesDegenerate[i].setUniqueness(true);
        fVecTFCoincidesUnique.push_back(fVecTFCoincidesDegenerate[i]);
      }
      else
      {
       fVecTFCoincidesDegenerate[i].setUniqueness(false);
      }
    }
    return fVecTFCoincidesUnique;
  }
  else
  {
    return fVecTFCoincidesUnique;
  }
}

std::pair<unsigned int,unsigned int> TFCoincidenceRecord::getMultiplicity(TFCoincidenceRecord::TFCoincidence const &cTFCoincidence)
{
  return {fMapUS[{cTFCoincidence.fDetUS, cTFCoincidence.fAcqUS, cTFCoincidence.fEventUS}].size(), fMapDS[{cTFCoincidence.fDetDS, cTFCoincidence.fAcqDS, cTFCoincidence.fEventDS}].size()};
}

unsigned int TFCoincidenceRecord::getMultiplicity(std::string const &cUSOrDS, unsigned int const &cDet, unsigned int const &cAcq, unsigned int const &cEvent)
{
  if(cUSOrDS=="US")
  {
    return fMapUS[{cDet,cAcq,cEvent}].size();
  }
  else
  {
    return fMapDS[{cDet,cAcq,cEvent}].size();
  }
}

void TFCoincidenceRecord::addCoincidence(unsigned int cDetUS, unsigned int cAcqUS, unsigned int cEventUS, unsigned int cDetDS, unsigned int cAcqDS, unsigned int cEventDS, double cTF)
{
  TFCoincidenceRecord::TFCoincidence tofCo = {cDetUS, cAcqUS, cEventUS, cDetDS, cAcqDS, cEventDS, cTF};
  fMapUS[{cDetUS,cAcqUS,cEventUS}].push_back(tofCo);
  fMapDS[{cDetDS,cAcqDS,cEventDS}].push_back(tofCo);
  fVecTFCoincidesDegenerate.push_back(tofCo);

  return;
}

void TFCoincidenceRecord::printDegenerate(std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex)
{
  if(fVecTFCoincidesUnique.size()==0)
  {
    getTFCoincidencesUnique();
  }

  std::vector<std::vector<AcquisitionXBTF>> vec_DetAcq(cMapDetectors.size());
  for(auto det : cMapDetectors)
  {
    if(det.second.getType()=="XBTF")
    {
      det.second.getAcquisitions(vec_DetAcq[cXBTFNameToIndex[det.first]]);     
    }
  }

  std::cout << "**************************************************************************" << std::endl;
  std::cout << "PRINTING DEGENERATE XBTF COINCIDENCES" << std::endl;

  for(unsigned int i = 0; i < fVecTFCoincidesDegenerate.size(); i++)
  {
    TFCoincidence currentTFCO = fVecTFCoincidesDegenerate[i];
    AcquisitionXBTF::EventRecordHR ER_US = vec_DetAcq[currentTFCO.fDetUS][currentTFCO.fAcqUS].getDataHR()[currentTFCO.fEventUS];
    AcquisitionXBTF::EventRecordHR ER_DS = vec_DetAcq[currentTFCO.fDetDS][currentTFCO.fAcqDS].getDataHR()[currentTFCO.fEventDS];
    std::pair<unsigned int, unsigned int> multiplicity = getMultiplicity(currentTFCO);

    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "UPSTREAM DET:  " << currentTFCO.fDetUS << ", ACQ " << currentTFCO.fAcqUS << ", EVENT " << currentTFCO.fEventUS << " TIME: " << ER_US.fSeconds_FullTime << std::endl;
    std::cout << "DOWNSTREM DET: " << currentTFCO.fDetDS << ", ACQ " << currentTFCO.fAcqDS << ", EVENT " << currentTFCO.fEventDS << " TIME: " << ER_DS.fSeconds_FullTime << std::endl;
    std::cout << "UPSTREAM MULTIPLICITY: " << multiplicity.first << ", DOWNSTREAM MULTIPLICITY: " << multiplicity.second << ", TF: " << currentTFCO.fTF << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
  }
  std::cout << "**************************************************************************" << std::endl;

  return;
};

void TFCoincidenceRecord::printUnique(std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex)
{
  if(fVecTFCoincidesUnique.size()==0)
  {
    getTFCoincidencesUnique();
  }

  std::vector<std::vector<AcquisitionXBTF>> vec_DetAcq(cMapDetectors.size());
  for(auto det : cMapDetectors)
  {
    if(det.second.getType()=="XBTF")
    {
      det.second.getAcquisitions(vec_DetAcq[cXBTFNameToIndex[det.first]]);     
    }
  }

  std::cout.precision(25);
  std::cout << "**************************************************************************" << std::endl;
  std::cout << "PRINTING UNIQUE XBTF COINCIDENCES" << std::endl;

  //for(unsigned int i = 0; i < fVecTFCoincidesUnique.size(); i++)
  for(unsigned int i = 0; i < 10; i++)
  {
    TFCoincidence currentTFCO = fVecTFCoincidesUnique[i];
    AcquisitionXBTF::EventRecordHR ER_US = vec_DetAcq[currentTFCO.fDetUS][currentTFCO.fAcqUS].getDataHR()[currentTFCO.fEventUS];
    AcquisitionXBTF::EventRecordHR ER_DS = vec_DetAcq[currentTFCO.fDetDS][currentTFCO.fAcqDS].getDataHR()[currentTFCO.fEventDS];
    std::pair<unsigned int, unsigned int> multiplicity = getMultiplicity(currentTFCO);

    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "UPSTREAM DET:  " << currentTFCO.fDetUS << ", ACQ " << currentTFCO.fAcqUS << ", EVENT " << currentTFCO.fEventUS << " TIME: " << ER_US.fSeconds_FullTime << std::endl;
    std::cout << "DOWNSTREM DET: " << currentTFCO.fDetDS << ", ACQ " << currentTFCO.fAcqDS << ", EVENT " << currentTFCO.fEventDS << " TIME: " << ER_DS.fSeconds_FullTime << std::endl;
    std::cout << "UPSTREAM MULTIPLICITY: " << multiplicity.first << ", DOWNSTREAM MULTIPLICITY: " << multiplicity.second << ", TF: " << currentTFCO.fTF << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
  }
  std::cout << "**************************************************************************" << std::endl;

  return;
};

void TFCoincidenceRecord::dumpDegenerate(std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TString const &cFilePathName)
{
  if(fVecTFCoincidesUnique.size()==0)
  {
    getTFCoincidencesUnique();
  }

  std::map<unsigned int,TString> map_XBTFIndexToName;
  std::vector<std::vector<AcquisitionXBTF>> vec_DetAcq(cMapDetectors.size());
  for(auto det : cMapDetectors)
  {
    if(det.second.getType()=="XBTF")
    {
      det.second.getAcquisitions(vec_DetAcq[cXBTFNameToIndex[det.first]]);     
      map_XBTFIndexToName[cXBTFNameToIndex[det.first]] = det.first;
    }
  }

  TFile *f_DumpFile = new TFile(cFilePathName,"UPDATE");
  TTree *t_DumpTree = new TTree("DumpBL_XBTF_Degen", "DumpBL_XBTF_Degen");

  TString out_US_DetName;               //UPSTREAM DETECTOR NAME.
  TString out_DS_DetName;               //DOWNSTREAM DETECTOR NAME.
  unsigned int out_US_AcqIndex;         //UPSTREAM TF ACQUISITION INDEX.
  unsigned int out_US_EventIndex;       //UPSTREAM TF EVENT INDEX.
  unsigned int out_US_Multiplicity;     //NUMBER OF TIMES THAT THIS EXACT EVENT IN THIS UPSTREAM DETECTOR HAS BEEN MATCHED TO AN EVENT TIME IN DOWNSTREAM DETECTORS.
  double       out_US_TriggerTimestamp; //TRIGGER TIMESTAMP OF THIS EVENT IN THE UPSTREAM DETECTOR. SECONDS.
  unsigned int out_DS_AcqIndex;         //DOWNSTREAM TF ACQUISITION INDEX.
  unsigned int out_DS_EventIndex;       //DOWNSTREAM TF EVENT INDEX.
  unsigned int out_DS_Multiplicity;     //NUMBER OF TIMES THAT THIS EXACT EVENT IN THIS DOWNSTREAM DETECTOR HAS BEEN MATCHED TO AN EVENT TIME IN UPSTREAM DETECTORS.
  double       out_DS_TriggerTimestamp; //TRIGGER TIMESTAMP OF THIS EVENT IN THE DOWNSTREAM DETECTOR. SECONDS.
  bool         out_IsUnique;            //FLAG FOR THE UNIQUENESS OF US/DS TIMESTAMP COMBINATION.
  double       out_TF;                  //DIFFERENCE IN UPSTREAM AND DOWNSTREAM TIMES.

  t_DumpTree->Branch("US_DetName",          &out_US_DetName);
  t_DumpTree->Branch("US_AcqIndex",         &out_US_AcqIndex);
  t_DumpTree->Branch("US_EventIndex",       &out_US_EventIndex);
  t_DumpTree->Branch("US_Multiplicity",     &out_US_Multiplicity);
  t_DumpTree->Branch("US_TriggerTimestamp", &out_US_TriggerTimestamp);
  t_DumpTree->Branch("DS_DetName",          &out_DS_DetName);
  t_DumpTree->Branch("DS_AcqIndex",         &out_DS_AcqIndex);
  t_DumpTree->Branch("DS_EventIndex",       &out_DS_EventIndex);
  t_DumpTree->Branch("DS_Multiplicity",     &out_DS_Multiplicity);
  t_DumpTree->Branch("DS_TriggerTimestamp", &out_DS_TriggerTimestamp);
  t_DumpTree->Branch("IsUnique",            &out_IsUnique);
  t_DumpTree->Branch("TF",                  &out_TF);

  for(unsigned int i = 0; i < fVecTFCoincidesDegenerate.size(); i++)
  {
    TFCoincidence currentTFCO = fVecTFCoincidesDegenerate[i];
    AcquisitionXBTF::EventRecordHR ER_US = vec_DetAcq[currentTFCO.fDetUS][currentTFCO.fAcqUS].getDataHR()[currentTFCO.fEventUS];
    AcquisitionXBTF::EventRecordHR ER_DS = vec_DetAcq[currentTFCO.fDetDS][currentTFCO.fAcqDS].getDataHR()[currentTFCO.fEventDS];
    std::pair<unsigned int, unsigned int> multiplicity = getMultiplicity(currentTFCO);

    out_US_DetName          = map_XBTFIndexToName[currentTFCO.fDetUS];
    out_US_AcqIndex         = currentTFCO.fAcqUS;
    out_US_EventIndex       = currentTFCO.fEventUS;
    out_US_Multiplicity     = multiplicity.first;
    out_US_TriggerTimestamp = ER_US.fSeconds_FullTime;

    out_DS_DetName          = map_XBTFIndexToName[currentTFCO.fDetDS];
    out_DS_AcqIndex         = currentTFCO.fAcqDS;
    out_DS_EventIndex       = currentTFCO.fEventDS;
    out_DS_Multiplicity     = multiplicity.second;
    out_DS_TriggerTimestamp = ER_DS.fSeconds_FullTime;

    out_IsUnique = currentTFCO.fUniqueness;
    out_TF       = currentTFCO.fTF;
    t_DumpTree->Fill();

    out_US_DetName = ""; out_US_AcqIndex = 0; out_US_EventIndex = 0; out_US_Multiplicity = 0; out_US_TriggerTimestamp = 0.; 
    out_DS_DetName = ""; out_DS_AcqIndex = 0; out_DS_EventIndex = 0; out_DS_Multiplicity = 0; out_DS_TriggerTimestamp = 0.; 
    out_IsUnique = 0;
    out_TF = 0.;
  }

  t_DumpTree->Write();

  f_DumpFile->Close();
  delete f_DumpFile;

  return;
}
void TFCoincidenceRecord::dumpUnique(std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TString const &cFilePathName)
{
  if(fVecTFCoincidesUnique.size()==0)
  {
    getTFCoincidencesUnique();
  }

  std::map<unsigned int,TString> map_XBTFIndexToName;
  std::vector<std::vector<AcquisitionXBTF>> vec_DetAcq(cMapDetectors.size());
  for(auto det : cMapDetectors)
  {
    if(det.second.getType()=="XBTF")
    {
      det.second.getAcquisitions(vec_DetAcq[cXBTFNameToIndex[det.first]]);     
      map_XBTFIndexToName[cXBTFNameToIndex[det.first]] = det.first;
    }
  }

  TFile *f_DumpFile = new TFile(cFilePathName,"UPDATE");
  TTree *t_DumpTree = new TTree("DumpBL_XBTF_Unique", "DumpBL_XBTF_Unique");

  TString out_US_DetName;
  TString out_DS_DetName;
  unsigned int out_US_AcqIndex;
  unsigned int out_US_EventIndex;
  unsigned int out_US_Multiplicity;
  double       out_US_TriggerTimestamp;
  unsigned int out_DS_AcqIndex;
  unsigned int out_DS_EventIndex;
  unsigned int out_DS_Multiplicity;
  double       out_DS_TriggerTimestamp;
  bool         out_IsUnique;
  double       out_TF;

  t_DumpTree->Branch("US_DetName",          &out_US_DetName);           //UPSTREAM DETECTOR NAME.
  t_DumpTree->Branch("DS_DetName",          &out_DS_DetName);           //DOWNSTREAM DETECTOR NAME.
  t_DumpTree->Branch("US_AcqIndex",         &out_US_AcqIndex);          //UPSTREAM TF ACQUISITION INDEX.
  t_DumpTree->Branch("US_EventIndex",       &out_US_EventIndex);        //UPSTREAM TF EVENT INDEX.
  t_DumpTree->Branch("US_Multiplicity",     &out_US_Multiplicity);      //NUMBER OF TIMES THAT THIS EXACT EVENT IN THIS UPSTREAM DETECTOR HAS BEEN MATCHED TO AN EVENT TIME IN DOWNSTREAM DETECTORS.
  t_DumpTree->Branch("US_TriggerTimestamp", &out_US_TriggerTimestamp);  //TRIGGER TIMESTAMP OF THIS EVENT IN THE UPSTREAM DETECTOR. SECONDS.
  t_DumpTree->Branch("DS_AcqIndex",         &out_DS_AcqIndex);          //DOWNSTREAM TF ACQUISITION INDEX.
  t_DumpTree->Branch("DS_EventIndex",       &out_DS_EventIndex);        //DOWNSTREAM TF EVENT INDEX.
  t_DumpTree->Branch("DS_Multiplicity",     &out_DS_Multiplicity);      //NUMBER OF TIMES THAT THIS EXACT EVENT IN THIS DOWNSTREAM DETECTOR HAS BEEN MATCHED TO AN EVENT TIME IN UPSTREAM DETECTORS.
  t_DumpTree->Branch("DS_TriggerTimestamp", &out_DS_TriggerTimestamp);  //TRIGGER TIMESTAMP OF THIS EVENT IN THE DOWNSTREAM DETECTOR. SECONDS.
  t_DumpTree->Branch("IsUnique",            &out_IsUnique);             //FLAG FOR THE UNIQUENESS OF US/DS TIMESTAMP COMBINATION.
  t_DumpTree->Branch("TF",                  &out_TF);                   //DIFFERENCE IN UPSTREAM AND DOWNSTREAM TIMES.

  for(unsigned int i = 0; i < fVecTFCoincidesUnique.size(); i++)
  {
    TFCoincidenceRecord::TFCoincidence currentTFCO = fVecTFCoincidesUnique[i];
    AcquisitionXBTF::EventRecordHR ER_US = vec_DetAcq[currentTFCO.fDetUS][currentTFCO.fAcqUS].getDataHR()[currentTFCO.fEventUS];
    AcquisitionXBTF::EventRecordHR ER_DS = vec_DetAcq[currentTFCO.fDetDS][currentTFCO.fAcqDS].getDataHR()[currentTFCO.fEventDS];
    std::pair<unsigned int, unsigned int> multiplicity = getMultiplicity(currentTFCO);

    out_US_DetName          = map_XBTFIndexToName[currentTFCO.fDetUS];
    out_US_AcqIndex         = currentTFCO.fAcqUS;
    out_US_EventIndex       = currentTFCO.fEventUS;
    out_US_Multiplicity     = multiplicity.first;
    out_US_TriggerTimestamp = ER_US.fSeconds_FullTime;

    out_DS_DetName          = map_XBTFIndexToName[currentTFCO.fDetDS];
    out_DS_AcqIndex         = currentTFCO.fAcqDS;
    out_DS_EventIndex       = currentTFCO.fEventDS;
    out_DS_Multiplicity     = multiplicity.second;
    out_DS_TriggerTimestamp = ER_DS.fSeconds_FullTime;

    out_IsUnique = currentTFCO.fUniqueness;
    out_TF       = currentTFCO.fTF;
    t_DumpTree->Fill();

    out_US_DetName = ""; out_US_AcqIndex = 0; out_US_EventIndex = 0; out_US_Multiplicity = 0; out_US_TriggerTimestamp = 0.; 
    out_DS_DetName = ""; out_DS_AcqIndex = 0; out_DS_EventIndex = 0; out_DS_Multiplicity = 0; out_DS_TriggerTimestamp = 0.; 
    out_IsUnique = 0;
    out_TF = 0.;
  }

  t_DumpTree->Write();

  f_DumpFile->Close();
  delete f_DumpFile;

  return;
}


std::vector<PROFCoincidenceRecord::PROFCoincidence> PROFCoincidenceRecord::getPROFCoincidencesDegenerate()
{
  if(fVecPROFCoincidesUnique.size()==0)
  {
    getPROFCoincidencesUnique();
  }
  return fVecPROFCoincidesDegenerate;
}

std::vector<PROFCoincidenceRecord::PROFCoincidence> PROFCoincidenceRecord::getPROFCoincidencesUnique()
{
  if(fVecPROFCoincidesUnique.size()==0)
  {
    for(unsigned int i = 0; i < fVecPROFCoincidesDegenerate.size(); i++)
    {
      std::array<long unsigned int,3> multiplicity = getMultiplicity(fVecPROFCoincidesDegenerate[i]);
      fVecPROFCoincidesDegenerate[i].setMultiplicity(multiplicity);
      if(multiplicity[0]==1 && multiplicity[1]==1 && multiplicity[2]==1)
      {
        fVecPROFCoincidesDegenerate[i].setUniqueness(true);
        fVecPROFCoincidesUnique.push_back(fVecPROFCoincidesDegenerate[i]);
      }
      else
      {
        fVecPROFCoincidesDegenerate[i].setUniqueness(false);
      }
    }
    return fVecPROFCoincidesUnique; 
  }
  else
  {
    return fVecPROFCoincidesUnique;
  }
}

std::array<long unsigned int,3> PROFCoincidenceRecord::getMultiplicity(PROFCoincidenceRecord::PROFCoincidence const &cPROFCoincidence)
{
  return {fMapPROF1[{cPROFCoincidence.fAcqPROF1, cPROFCoincidence.fEventPROF1}].size(),fMapPROF2[{cPROFCoincidence.fAcqPROF2, cPROFCoincidence.fEventPROF2}].size(),
          fMapPROF3[{cPROFCoincidence.fAcqPROF1, cPROFCoincidence.fEventPROF3}].size()};
}

unsigned int PROFCoincidenceRecord::getMultiplicity(unsigned int const &cDet, unsigned int const &cAcq, unsigned int const &cEvent)
{
  if(cDet==1)
  {
    return fMapPROF1[{cAcq,cEvent}].size();
  }
  else if(cDet==2)
  {
    return fMapPROF2[{cAcq,cEvent}].size();
  }
  else
  {
    return fMapPROF3[{cAcq,cEvent}].size();
  }
}

void PROFCoincidenceRecord::addCoincidence(unsigned int cAcqPROF1, unsigned int cEventPROF1, unsigned int cAcqPROF2, unsigned int cEventPROF2, unsigned int cAcqPROF3, unsigned int cEventPROF3,
                                           std::vector<double> cCosTheta, std::vector<double> cTheta, std::vector<double> cMomentum)
{
  PROFCoincidenceRecord::PROFCoincidence profCo(cAcqPROF1, cEventPROF1, cAcqPROF2, cEventPROF2, cAcqPROF3, cEventPROF3);
  profCo.addDeflectionInfo(cCosTheta, cTheta, cMomentum);
  fMapPROF1[{cAcqPROF1,cEventPROF1}].push_back(profCo);
  fMapPROF2[{cAcqPROF2,cEventPROF2}].push_back(profCo);
  fMapPROF3[{cAcqPROF3,cEventPROF3}].push_back(profCo);
  fVecPROFCoincidesDegenerate.push_back(profCo);

  return;
}

void PROFCoincidenceRecord::printDegenerate(std::map<std::string,Detector> &cMapDetectors)
{
  if(fVecPROFCoincidesUnique.size()==0)
  {
    getPROFCoincidencesUnique();
  }

  std::vector<AcquisitionXBPF> vec_DetAcqPROF1;
  cMapDetectors["XBPF022697"].getAcquisitions(vec_DetAcqPROF1);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF2;
  cMapDetectors["XBPF022701"].getAcquisitions(vec_DetAcqPROF2);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF3;
  cMapDetectors["XBPF022702"].getAcquisitions(vec_DetAcqPROF3);

  std::cout << "**************************************************************************" << std::endl;
  std::cout << "PRINTING DEGENERATE PROF COINCIDENCES" << std::endl;

  for(unsigned int i = 0; i < fVecPROFCoincidesDegenerate.size(); i++)
  {
    PROFCoincidence currentPROFCO = fVecPROFCoincidesDegenerate[i];
    AcquisitionXBPF::EventRecordHR ER_PROF1 = vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getDataHR()[currentPROFCO.fEventPROF1];
    AcquisitionXBPF::EventRecordHR ER_PROF2 = vec_DetAcqPROF2[currentPROFCO.fAcqPROF2].getDataHR()[currentPROFCO.fEventPROF2];
    AcquisitionXBPF::EventRecordHR ER_PROF3 = vec_DetAcqPROF3[currentPROFCO.fAcqPROF3].getDataHR()[currentPROFCO.fEventPROF3];

    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "PROF1: ACQ " << currentPROFCO.fAcqPROF1 << ", EVENT " << currentPROFCO.fEventPROF1 << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[0]  << " TIME: " << ER_PROF1.fTriggerTimestamp
              << "\n    FIBRES ";
    for(unsigned int j = 0; j < ER_PROF1.fFibresList.size(); j++)
    {
      std::cout << ER_PROF1.fFibresList[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "PROF2: ACQ " << currentPROFCO.fAcqPROF2 << ", EVENT " << currentPROFCO.fEventPROF2  << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[1] << " TIME: " << ER_PROF2.fTriggerTimestamp
              << "\n    FIBRES ";
    for(unsigned int j = 0; j < ER_PROF2.fFibresList.size(); j++)
    {
      std::cout << ER_PROF2.fFibresList[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "PROF3: ACQ " << currentPROFCO.fAcqPROF3 << ", EVENT " << currentPROFCO.fEventPROF3  << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[2] << " TIME: " << ER_PROF3.fTriggerTimestamp
              << "\n    FIBRES ";
    for(unsigned int j = 0; j < ER_PROF3.fFibresList.size(); j++)
    {
      std::cout << ER_PROF3.fFibresList[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "    NUMBER OF POSSIBLE FIBRE COMBINATIONS: " << currentPROFCO.fMomentum.size() << std::endl;
    for(unsigned int j = 0; j < currentPROFCO.getCosTheta().size(); j++)
    {
      std::cout << "    COSTHETA: " << currentPROFCO.getCosTheta()[j] << ", THETA (DEGREES): " << currentPROFCO.getTheta()[j]*radToDeg << ", MOMENTUM (GeV): " << currentPROFCO.getMomentum()[j] << std::endl;
    }
    std::cout << "-------------------------------------------------------------------------" << std::endl;
  }
  std::cout << "**************************************************************************" << std::endl;

  return;
};

void PROFCoincidenceRecord::printUnique(std::map<std::string,Detector> &cMapDetectors)
{
  if(fVecPROFCoincidesUnique.size()==0)
  {
    getPROFCoincidencesUnique();
  }

  std::vector<AcquisitionXBPF> vec_DetAcqPROF1;
  cMapDetectors["XBPF022697"].getAcquisitions(vec_DetAcqPROF1);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF2;
  cMapDetectors["XBPF022701"].getAcquisitions(vec_DetAcqPROF2);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF3;
  cMapDetectors["XBPF022702"].getAcquisitions(vec_DetAcqPROF3);

  std::cout << "**************************************************************************" << std::endl;
  std::cout << "PRINTING UNIQUE PROF COINCIDENCES" << std::endl;

  for(unsigned int i = 0; i < 10; i++)
  //for(unsigned int i = 0; i < fVecPROFCoincidesUnique.size(); i++)
  {
    PROFCoincidence currentPROFCO = fVecPROFCoincidesUnique[i];
    AcquisitionXBPF::EventRecordHR ER_PROF1 = vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getDataHR()[currentPROFCO.fEventPROF1];
    AcquisitionXBPF::EventRecordHR ER_PROF2 = vec_DetAcqPROF2[currentPROFCO.fAcqPROF2].getDataHR()[currentPROFCO.fEventPROF2];
    AcquisitionXBPF::EventRecordHR ER_PROF3 = vec_DetAcqPROF3[currentPROFCO.fAcqPROF3].getDataHR()[currentPROFCO.fEventPROF3];

    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "PROF1:   ACQ " << currentPROFCO.fAcqPROF1 << ", EVENT " << currentPROFCO.fEventPROF1 << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[0] <<  " TIME: " << ER_PROF1.fTriggerTimestamp
              << "\n    FIBRES ";
    for(unsigned int j = 0; j < ER_PROF1.fFibresList.size(); j++)
    {
      std::cout << ER_PROF1.fFibresList[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "PROF2:   ACQ " << currentPROFCO.fAcqPROF2 << ", EVENT " << currentPROFCO.fEventPROF2  << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[1] << " TIME: " << ER_PROF2.fTriggerTimestamp
              << "\n    FIBRES ";
    for(unsigned int j = 0; j < ER_PROF2.fFibresList.size(); j++)
    {
      std::cout << ER_PROF2.fFibresList[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "PROF3:   ACQ " << currentPROFCO.fAcqPROF3 << ", EVENT " << currentPROFCO.fEventPROF3  << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[2] << " TIME: " << ER_PROF3.fTriggerTimestamp
              << "\n    FIBRES ";
    for(unsigned int j = 0; j < ER_PROF3.fFibresList.size(); j++)
    {
      std::cout << ER_PROF3.fFibresList[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "    NUMBER OF POSSIBLE FIBRE COMBINATIONS: " << currentPROFCO.getNFibreCombinations() << std::endl;
    for(unsigned int j = 0; j < currentPROFCO.getCosTheta().size(); j++)
    {
      std::cout << "    COSTHETA: " << currentPROFCO.getCosTheta()[j] << ", THETA (DEGREES): " << currentPROFCO.getTheta()[j]*radToDeg << ", MOMENTUM (GeV): " << currentPROFCO.getMomentum()[j] << std::endl;
    }
    std::cout << "-------------------------------------------------------------------------" << std::endl;
  }
  std::cout << "**************************************************************************" << std::endl;

  return;
};

void PROFCoincidenceRecord::dumpDegenerate(std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName)
{
  if(fVecPROFCoincidesUnique.size()==0)
  {
    getPROFCoincidencesUnique();
  }

  //IN THIS TREE, EACH ENTRY IS A POSSIBLE TRIPLE COINCIDENCE. FOR EACH POSSIBLE FIBRE COMBINATION, A MOMENTUM IS CALCULATED.
  TFile *f_DumpFile = new TFile(cFilePathName,"UPDATE");
  TTree *t_DumpTree = new TTree("DumpBL_XBPF_Degen", "DumpBL_XBPF_Degen");

  std::vector<AcquisitionXBPF> vec_DetAcqPROF1;
  cMapDetectors["XBPF022697"].getAcquisitions(vec_DetAcqPROF1);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF2;
  cMapDetectors["XBPF022701"].getAcquisitions(vec_DetAcqPROF2);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF3;
  cMapDetectors["XBPF022702"].getAcquisitions(vec_DetAcqPROF3);

  unsigned int out_PROF1_AcqIndex;               //ACQUISITION INDEX OF PROF1.
  unsigned int out_PROF1_EventIndex;             //EVENT INDEX OF PROF1.
  unsigned int out_PROF1_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF2 AND PROF3 TIMING COINCIDENCES.
  long long    out_PROF1_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF1. NANOSECONDS.
  unsigned int out_PROF1_NFibresHit;             //NUMBER OF PROF1 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_PROF1_FibresHit; //LIST OF PROF2 FIBRES HIT IN THIS EVENT.
  unsigned int out_PROF2_AcqIndex;               //ACQUISITION INDEX OF PROF2.
  unsigned int out_PROF2_EventIndex;             //EVENT INDEX OF PROF2.
  unsigned int out_PROF2_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF1 AND PROF3 TIMING COINCIDENCES.
  long long    out_PROF2_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF2. NANOSECONDS.
  unsigned int out_PROF2_NFibresHit;             //NUMBER OF PROF2 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_PROF2_FibresHit; //LIST OF PROF2 FIBRES HIT IN THIS EVENT.
  unsigned int out_PROF3_AcqIndex;               //ACQUISITION INDEX OF PROF3. 
  unsigned int out_PROF3_EventIndex;             //EVENT INDEX OF PROF3.
  unsigned int out_PROF3_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF1 AND PROF2 TIMING COINCIDENCES.
  long long    out_PROF3_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF3. NANOSECONDS.
  unsigned int out_PROF3_NFibresHit;             //NUMBER OF PROF3 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_PROF3_FibresHit; //LIST OF PROF3 FIBRES HIT IN THIS EVENT.
  unsigned int out_NFibreCombinations;           //NUMBER OF POSSIBLE FIBRES COMBINATIONS ACROSS THE THREE PROFILERS.
  double       out_Current;                      //CURRENT OF THE BENDING MAGNET USED IN THE CALCULATION OF BFIELD, AMPS.
  bool         out_IsUnique;                     //FLAG AS TO WHEN THIS COINCIDENCE IS TIME UNIQUE.
  std::vector<double> out_CosTheta;              //CALCULATED COSTHETA.
  std::vector<double> out_Theta;                 //CALCULATED THETA, RADIANS.
  std::vector<double> out_Momentum;              //CALCULATED MOMENTUM, GeV.

  t_DumpTree->Branch("PROF1_AcqIndex",         &out_PROF1_AcqIndex);
  t_DumpTree->Branch("PROF1_EventIndex",       &out_PROF1_EventIndex);
  t_DumpTree->Branch("PROF1_Multiplicity",     &out_PROF1_Multiplicity);
  t_DumpTree->Branch("PROF1_TriggerTimestamp", &out_PROF1_TriggerTimestamp);
  t_DumpTree->Branch("PROF1_NFibresHit",       &out_PROF1_NFibresHit);
  t_DumpTree->Branch("PROF1_FibresHit",        &out_PROF1_FibresHit);
  t_DumpTree->Branch("PROF2_AcqIndex",         &out_PROF2_AcqIndex);
  t_DumpTree->Branch("PROF2_EventIndex",       &out_PROF2_EventIndex);
  t_DumpTree->Branch("PROF2_Multiplicity",     &out_PROF2_Multiplicity);
  t_DumpTree->Branch("PROF2_TriggerTimestamp", &out_PROF2_TriggerTimestamp);
  t_DumpTree->Branch("PROF2_NFibresHit",       &out_PROF2_NFibresHit);
  t_DumpTree->Branch("PROF2_FibresHit",        &out_PROF2_FibresHit);
  t_DumpTree->Branch("PROF3_AcqIndex",         &out_PROF3_AcqIndex);
  t_DumpTree->Branch("PROF3_EventIndex",       &out_PROF3_EventIndex);
  t_DumpTree->Branch("PROF3_Multiplicity",     &out_PROF3_Multiplicity);
  t_DumpTree->Branch("PROF3_TriggerTimestamp", &out_PROF3_TriggerTimestamp);
  t_DumpTree->Branch("PROF3_NFibresHit",       &out_PROF3_NFibresHit);
  t_DumpTree->Branch("PROF3_FibresHit",        &out_PROF3_FibresHit);
  t_DumpTree->Branch("NFibreCombinations",     &out_NFibreCombinations);
  t_DumpTree->Branch("IsUnqiue",               &out_IsUnique);
  t_DumpTree->Branch("Current",                &out_Current);
  t_DumpTree->Branch("CosTheta",               &out_CosTheta);
  t_DumpTree->Branch("Theta",                  &out_Theta);
  t_DumpTree->Branch("Momentum",               &out_Momentum);

  for(unsigned int i = 0; i < fVecPROFCoincidesDegenerate.size(); i++)
  {
    PROFCoincidence currentPROFCO = fVecPROFCoincidesDegenerate[i];
    AcquisitionXBPF::EventRecordHR ER_PROF1 = vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getDataHR()[currentPROFCO.fEventPROF1];
    AcquisitionXBPF::EventRecordHR ER_PROF2 = vec_DetAcqPROF2[currentPROFCO.fAcqPROF2].getDataHR()[currentPROFCO.fEventPROF2];
    AcquisitionXBPF::EventRecordHR ER_PROF3 = vec_DetAcqPROF3[currentPROFCO.fAcqPROF3].getDataHR()[currentPROFCO.fEventPROF3];

    out_PROF1_AcqIndex = currentPROFCO.fAcqPROF1; out_PROF1_EventIndex = currentPROFCO.fEventPROF1; out_PROF1_Multiplicity = currentPROFCO.fMultiplicity[0]; 
    out_PROF1_TriggerTimestamp = ER_PROF1.fTriggerTimestamp; out_PROF1_NFibresHit = ER_PROF1.fFibresList.size();
    for(unsigned int j = 0; j < ER_PROF1.fFibresList.size(); j++)
    {
      out_PROF1_FibresHit.push_back(ER_PROF1.fFibresList[j]);
    }

    out_PROF2_AcqIndex = currentPROFCO.fAcqPROF2; out_PROF2_EventIndex = currentPROFCO.fEventPROF2; out_PROF2_Multiplicity = currentPROFCO.fMultiplicity[1]; 
    out_PROF2_TriggerTimestamp = ER_PROF2.fTriggerTimestamp; out_PROF2_NFibresHit = ER_PROF2.fFibresList.size();
    for(unsigned int j = 0; j < ER_PROF2.fFibresList.size(); j++)
    {
      out_PROF2_FibresHit.push_back(ER_PROF2.fFibresList[j]);
    }

    out_PROF3_AcqIndex = currentPROFCO.fAcqPROF3; out_PROF3_EventIndex = currentPROFCO.fEventPROF3; out_PROF3_Multiplicity = currentPROFCO.fMultiplicity[2]; 
    out_PROF3_TriggerTimestamp = ER_PROF3.fTriggerTimestamp; out_PROF3_NFibresHit = ER_PROF3.fFibresList.size();
    for(unsigned int j = 0; j < ER_PROF3.fFibresList.size(); j++)
    {
      out_PROF3_FibresHit.push_back(ER_PROF3.fFibresList[j]);
    }

    out_NFibreCombinations = currentPROFCO.getNFibreCombinations();
    out_IsUnique = currentPROFCO.fUniqueness;
    out_Current = vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getCurrent();
    for(unsigned int j = 0; j < currentPROFCO.getCosTheta().size(); j++)
    {
      out_CosTheta.push_back(currentPROFCO.getCosTheta()[j]);
      out_Theta   .push_back(currentPROFCO.getTheta()[j]);
      out_Momentum.push_back(currentPROFCO.getMomentum()[j]);
    }
    t_DumpTree->Fill();
    out_PROF1_AcqIndex = 0; out_PROF1_EventIndex = 0; out_PROF1_Multiplicity = 0; out_PROF1_TriggerTimestamp = 0.; out_PROF1_NFibresHit = 0; out_PROF1_FibresHit.clear();
    out_PROF2_AcqIndex = 0; out_PROF2_EventIndex = 0; out_PROF2_Multiplicity = 0; out_PROF2_TriggerTimestamp = 0.; out_PROF2_NFibresHit = 0; out_PROF2_FibresHit.clear();
    out_PROF3_AcqIndex = 0; out_PROF3_EventIndex = 0; out_PROF3_Multiplicity = 0; out_PROF3_TriggerTimestamp = 0.; out_PROF3_NFibresHit = 0; out_PROF3_FibresHit.clear();
    out_IsUnique = 0.; out_Current = 0.; out_NFibreCombinations = 0; out_CosTheta.clear(); out_Theta.clear(); out_Momentum.clear();
  }

  t_DumpTree->Write();

  f_DumpFile->Close();
  delete f_DumpFile;

  return;
}

void PROFCoincidenceRecord::dumpUnique(std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName)
{
  if(fVecPROFCoincidesUnique.size()==0)
  {
    getPROFCoincidencesUnique();
  }

  TFile *f_DumpFile = new TFile(cFilePathName,"UPDATE");
  TTree *t_DumpTree = new TTree("DumpBL_XBPF_Unique", "DumpBL_XBPF_Unique");

  std::vector<AcquisitionXBPF> vec_DetAcqPROF1;
  cMapDetectors["XBPF022697"].getAcquisitions(vec_DetAcqPROF1);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF2;
  cMapDetectors["XBPF022701"].getAcquisitions(vec_DetAcqPROF2);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF3;
  cMapDetectors["XBPF022702"].getAcquisitions(vec_DetAcqPROF3);

  unsigned int out_PROF1_AcqIndex;               //ACQUISITION INDEX OF PROF1.
  unsigned int out_PROF1_EventIndex;             //EVENT INDEX OF PROF1.
  unsigned int out_PROF1_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF2 AND PROF3 TIMING COINCIDENCES.
  long long    out_PROF1_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF1. NANOSECONDS.
  unsigned int out_PROF1_NFibresHit;             //NUMBER OF PROF1 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_PROF1_FibresHit; //LIST OF PROF2 FIBRES HIT IN THIS EVENT.
  unsigned int out_PROF2_AcqIndex;               //ACQUISITION INDEX OF PROF2.
  unsigned int out_PROF2_EventIndex;             //EVENT INDEX OF PROF2.
  unsigned int out_PROF2_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF1 AND PROF3 TIMING COINCIDENCES.
  long long    out_PROF2_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF2. NANOSECONDS.
  unsigned int out_PROF2_NFibresHit;             //NUMBER OF PROF2 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_PROF2_FibresHit; //LIST OF PROF2 FIBRES HIT IN THIS EVENT.
  unsigned int out_PROF3_AcqIndex;               //ACQUISITION INDEX OF PROF3. 
  unsigned int out_PROF3_EventIndex;             //EVENT INDEX OF PROF3.
  unsigned int out_PROF3_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF1 AND PROF2 TIMING COINCIDENCES.
  long long    out_PROF3_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF3. NANOSECONDS.
  unsigned int out_PROF3_NFibresHit;             //NUMBER OF PROF3 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_PROF3_FibresHit; //LIST OF PROF3 FIBRES HIT IN THIS EVENT.
  unsigned int out_NFibreCombinations;           //NUMBER OF POSSIBLE FIBRES COMBINATIONS ACROSS THE THREE PROFILERS.
  bool         out_IsUnique;                     //CURRENT OF THE BENDING MAGNET USED IN THE CALCULATION OF BFIELD, AMPS.
  double       out_Current;                      //FLAG AS TO WHEN THIS COINCIDENCE IS TIME UNIQUE.
  std::vector<double> out_CosTheta;              //CALCULATED COSTHETA.
  std::vector<double> out_Theta;                 //CALCULATED THETA, RADIANS.
  std::vector<double> out_Momentum;              //CALCULATED MOMENTUM, GeV.

  t_DumpTree->Branch("PROF1_AcqIndex",         &out_PROF1_AcqIndex);
  t_DumpTree->Branch("PROF1_EventIndex",       &out_PROF1_EventIndex);
  t_DumpTree->Branch("PROF1_Multiplicity",     &out_PROF1_Multiplicity);
  t_DumpTree->Branch("PROF1_TriggerTimestamp", &out_PROF1_TriggerTimestamp);
  t_DumpTree->Branch("PROF1_NFibresHit",       &out_PROF1_NFibresHit);
  t_DumpTree->Branch("PROF1_FibresHit",        &out_PROF1_FibresHit);
  t_DumpTree->Branch("PROF2_AcqIndex",         &out_PROF2_AcqIndex);
  t_DumpTree->Branch("PROF2_EventIndex",       &out_PROF2_EventIndex);
  t_DumpTree->Branch("PROF2_Multiplicity",     &out_PROF2_Multiplicity);
  t_DumpTree->Branch("PROF2_TriggerTimestamp", &out_PROF2_TriggerTimestamp);
  t_DumpTree->Branch("PROF2_NFibresHit",       &out_PROF2_NFibresHit);
  t_DumpTree->Branch("PROF2_FibresHit",        &out_PROF2_FibresHit);
  t_DumpTree->Branch("PROF3_AcqIndex",         &out_PROF3_AcqIndex);
  t_DumpTree->Branch("PROF3_EventIndex",       &out_PROF3_EventIndex);
  t_DumpTree->Branch("PROF3_Multiplicity",     &out_PROF3_Multiplicity);
  t_DumpTree->Branch("PROF3_TriggerTimestamp", &out_PROF3_TriggerTimestamp);
  t_DumpTree->Branch("PROF3_NFibresHit",       &out_PROF3_NFibresHit);
  t_DumpTree->Branch("PROF3_FibresHit",        &out_PROF3_FibresHit);
  t_DumpTree->Branch("NFibreCombinations",     &out_NFibreCombinations);
  t_DumpTree->Branch("IsUnique",               &out_IsUnique);
  t_DumpTree->Branch("Current",                &out_Current);
  t_DumpTree->Branch("CosTheta",               &out_CosTheta);
  t_DumpTree->Branch("Theta",                  &out_Theta);
  t_DumpTree->Branch("Momentum",               &out_Momentum);

  for(unsigned int i = 0; i < fVecPROFCoincidesUnique.size(); i++)
  {
    PROFCoincidence currentPROFCO = fVecPROFCoincidesUnique[i];
    AcquisitionXBPF::EventRecordHR ER_PROF1 = vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getDataHR()[currentPROFCO.fEventPROF1];
    AcquisitionXBPF::EventRecordHR ER_PROF2 = vec_DetAcqPROF2[currentPROFCO.fAcqPROF2].getDataHR()[currentPROFCO.fEventPROF2];
    AcquisitionXBPF::EventRecordHR ER_PROF3 = vec_DetAcqPROF3[currentPROFCO.fAcqPROF3].getDataHR()[currentPROFCO.fEventPROF3];

    out_PROF1_AcqIndex = currentPROFCO.fAcqPROF1; out_PROF1_EventIndex = currentPROFCO.fEventPROF1; out_PROF1_Multiplicity = currentPROFCO.fMultiplicity[0]; 
    out_PROF1_TriggerTimestamp = ER_PROF1.fTriggerTimestamp; out_PROF1_NFibresHit = ER_PROF1.fFibresList.size();
    for(unsigned int j = 0; j < ER_PROF1.fFibresList.size(); j++)
    {
      out_PROF1_FibresHit.push_back(ER_PROF1.fFibresList[j]);
    }

    out_PROF2_AcqIndex = currentPROFCO.fAcqPROF2; out_PROF2_EventIndex = currentPROFCO.fEventPROF2; out_PROF2_Multiplicity = currentPROFCO.fMultiplicity[1]; 
    out_PROF2_TriggerTimestamp = ER_PROF2.fTriggerTimestamp; out_PROF2_NFibresHit = ER_PROF2.fFibresList.size();
    for(unsigned int j = 0; j < ER_PROF2.fFibresList.size(); j++)
    {
      out_PROF2_FibresHit.push_back(ER_PROF2.fFibresList[j]);
    }

    out_PROF3_AcqIndex = currentPROFCO.fAcqPROF3; out_PROF3_EventIndex = currentPROFCO.fEventPROF3; out_PROF3_Multiplicity = currentPROFCO.fMultiplicity[2]; 
    out_PROF3_TriggerTimestamp = ER_PROF3.fTriggerTimestamp; out_PROF3_NFibresHit = ER_PROF3.fFibresList.size();
    for(unsigned int j = 0; j < ER_PROF3.fFibresList.size(); j++)
    {
      out_PROF3_FibresHit.push_back(ER_PROF3.fFibresList[j]);
    }

    out_NFibreCombinations = currentPROFCO.getNFibreCombinations();
    out_IsUnique = currentPROFCO.fUniqueness;
    out_Current = vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getCurrent();
    for(unsigned int j = 0; j < currentPROFCO.getCosTheta().size(); j++)
    {
      out_CosTheta.push_back(currentPROFCO.getCosTheta()[j]);
      out_Theta   .push_back(currentPROFCO.getTheta()[j]);
      out_Momentum.push_back(currentPROFCO.getMomentum()[j]);
    }
    t_DumpTree->Fill();
    out_PROF1_AcqIndex = 0; out_PROF1_EventIndex = 0; out_PROF1_Multiplicity = 0; out_PROF1_TriggerTimestamp = 0.; out_PROF1_NFibresHit = 0; out_PROF1_FibresHit.clear();
    out_PROF2_AcqIndex = 0; out_PROF2_EventIndex = 0; out_PROF2_Multiplicity = 0; out_PROF2_TriggerTimestamp = 0.; out_PROF2_NFibresHit = 0; out_PROF2_FibresHit.clear();
    out_PROF3_AcqIndex = 0; out_PROF3_EventIndex = 0; out_PROF3_Multiplicity = 0; out_PROF3_TriggerTimestamp = 0.; out_PROF3_NFibresHit = 0; out_PROF3_FibresHit.clear();
    out_Current = 0.; out_NFibreCombinations = 0; out_CosTheta.clear(); out_Theta.clear(); out_Momentum.clear();
  }

  t_DumpTree->Write();

  f_DumpFile->Close();
  delete f_DumpFile;

  return;
}

std::vector<CombinedCoincidenceRecord::TFPROFCoincidence> CombinedCoincidenceRecord::getTFPROFCoincidencesDegenerate()
{
  if(fTFPROFCoincidencesUnique.size()==0)
  {
    getTFPROFCoincidencesUnique();
  }

  return fTFPROFCoincidencesDegenerate;
}

std::vector<CombinedCoincidenceRecord::TFPROFCoincidence> CombinedCoincidenceRecord::getTFPROFCoincidencesUnique()
{
  if(fTFPROFCoincidencesUnique.size()==0)
  {
    for(unsigned int i = 0; i < fTFPROFCoincidencesDegenerate.size(); i++)
    {
      if(fTFPROFCoincidencesDegenerate[i].getMultiplicity()==1)
      {
        fTFPROFCoincidencesUnique.push_back(fTFPROFCoincidencesDegenerate[i]);
      }
    }
  }

  return fTFPROFCoincidencesUnique;
}

unsigned int CombinedCoincidenceRecord::getMultiplicity(PROFCoincidenceRecord::PROFCoincidence const & cPROFCoincidence)
{
  return fMapPROFCount[{cPROFCoincidence.fAcqPROF1, cPROFCoincidence.fEventPROF1, cPROFCoincidence.fAcqPROF2, cPROFCoincidence.fEventPROF2, cPROFCoincidence.fAcqPROF3, cPROFCoincidence.fEventPROF3}];
}

void CombinedCoincidenceRecord::addTFPROFCoincidence(CombinedCoincidenceRecord::TFPROFCoincidence cTFPROFCoincidence)
{
  fTFPROFCoincidencesDegenerate.push_back(cTFPROFCoincidence);
  return;
}

void CombinedCoincidenceRecord::addMultiplicityInfo(std::array<unsigned int,6> &cPROFCoincidenceIncides)
{
  fMapPROFCount[cPROFCoincidenceIncides]++;
  return;
}

void CombinedCoincidenceRecord::dumpDegenerate(std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TString const &cFilePathName,
                                               TFCoincidenceRecord &cTFCoincidenceRecord)
{
  gInterpreter->GenerateDictionary("vector<vector<UInt_t> >","vector");
  
  if(fTFPROFCoincidencesUnique.size()==0)
  {
    getTFPROFCoincidencesUnique();  
  }

  std::map<unsigned int,TString> map_XBTFIndexToName;
  std::vector<std::vector<AcquisitionXBTF>> vec_DetAcq(cMapDetectors.size());
  for(auto det : cMapDetectors)
  {
    if(det.second.getType()=="XBTF")
    {
      det.second.getAcquisitions(vec_DetAcq[cXBTFNameToIndex[det.first]]);     
      map_XBTFIndexToName[cXBTFNameToIndex[det.first]] = det.first;
    }
  }

  std::vector<AcquisitionXBPF> vec_DetAcqPROF1;
  cMapDetectors["XBPF022697"].getAcquisitions(vec_DetAcqPROF1);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF2;
  cMapDetectors["XBPF022701"].getAcquisitions(vec_DetAcqPROF2);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF3;
  cMapDetectors["XBPF022702"].getAcquisitions(vec_DetAcqPROF3);

  TFile *f_DumpFile = new TFile(cFilePathName,"UPDATE");
  TTree *t_DumpTree = new TTree("DumpBL_XBTFXBPF_Degen", "DumpBL_XBTFXBPF_Degen");

  TString out_US_DetName;                                     
  TString out_DS_DetName;
  unsigned int out_US_AcqIndex;
  unsigned int out_US_EventIndex;
  unsigned int out_US_Multiplicity;
  double       out_US_TriggerTimestamp;
  unsigned int out_DS_AcqIndex;
  unsigned int out_DS_EventIndex;
  unsigned int out_DS_Multiplicity;
  double       out_DS_TriggerTimestamp;
  bool         out_IsTFUnique;
  double       out_TF;
  unsigned int out_Multiplicty;
  std::vector<unsigned int> out_PROFCoMultiplicity;
  std::vector<unsigned int> out_PROF1_AcqIndex;               //ACQUISITION INDEX OF PROF1.
  std::vector<unsigned int> out_PROF1_EventIndex;             //EVENT INDEX OF PROF1.
  std::vector<unsigned int> out_PROF1_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF2 AND PROF3 TIMING COINCIDENCES.
  std::vector<long long>    out_PROF1_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF1. NANOSECONDS.
  std::vector<unsigned int> out_PROF1_NFibresHit;             //NUMBER OF PROF1 FIBRES HIT IN THIS EVENT.
  std::vector<std::vector<unsigned int>> out_PROF1_FibresHit; //LIST OF PROF1 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_PROF2_AcqIndex;               //ACQUISITION INDEX OF PROF2.
  std::vector<unsigned int> out_PROF2_EventIndex;             //EVENT INDEX OF PROF2.
  std::vector<unsigned int> out_PROF2_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF1 AND PROF3 TIMING COINCIDENCES.
  std::vector<long long>    out_PROF2_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF2. NANOSECONDS.
  std::vector<unsigned int> out_PROF2_NFibresHit;             //NUMBER OF PROF2 FIBRES HIT IN THIS EVENT.
  std::vector<std::vector<unsigned int>> out_PROF2_FibresHit; //LIST OF PROF2 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_PROF3_AcqIndex;               //ACQUISITION INDEX OF PROF3. 
  std::vector<unsigned int> out_PROF3_EventIndex;             //EVENT INDEX OF PROF3.
  std::vector<unsigned int> out_PROF3_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF1 AND PROF2 TIMING COINCIDENCES.
  std::vector<long long>    out_PROF3_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF3. NANOSECONDS.
  std::vector<unsigned int> out_PROF3_NFibresHit;             //NUMBER OF PROF3 FIBRES HIT IN THIS EVENT.
  std::vector<std::vector<unsigned int>> out_PROF3_FibresHit; //LIST OF PROF3 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_NFibreCombinations;           //NUMBER OF POSSIBLE FIBRES COMBINATIONS ACROSS THE THREE PROFILERS.
  std::vector<bool>         out_IsPROFUnique;                 //CURRENT OF THE BENDING MAGNET USED IN THE CALCULATION OF BFIELD, AMPS.
  std::vector<double>       out_Current;                      //FLAG AS TO WHEN THIS COINCIDENCE IS TIME UNIQUE.
  std::vector<std::vector<double>> out_CosTheta;              //CALCULATED COSTHETA.
  std::vector<std::vector<double>> out_Theta;                 //CALCULATED THETA, RADIANS.
  std::vector<std::vector<double>> out_Momentum;              //CALCULATED MOMENTUM, GeV.
  std::vector<std::vector<double>> out_Mass;                  //CALCULATED MASS.

  t_DumpTree->Branch("US_DetName",             &out_US_DetName);           
  t_DumpTree->Branch("DS_DetName",             &out_DS_DetName);           
  t_DumpTree->Branch("US_AcqIndex",            &out_US_AcqIndex);          
  t_DumpTree->Branch("US_EventIndex",          &out_US_EventIndex);        
  t_DumpTree->Branch("US_Multiplicity",        &out_US_Multiplicity);      
  t_DumpTree->Branch("US_TriggerTimestamp",    &out_US_TriggerTimestamp);  
  t_DumpTree->Branch("DS_AcqIndex",            &out_DS_AcqIndex);          
  t_DumpTree->Branch("DS_EventIndex",          &out_DS_EventIndex);        
  t_DumpTree->Branch("DS_Multiplicity",        &out_DS_Multiplicity);      
  t_DumpTree->Branch("DS_TriggerTimestamp",    &out_DS_TriggerTimestamp);  
  t_DumpTree->Branch("IsTFUnique",             &out_IsTFUnique);             
  t_DumpTree->Branch("TF",                     &out_TF);                   
  t_DumpTree->Branch("Multiplicity",           &out_Multiplicty);
  t_DumpTree->Branch("PROF1_AcqIndex",         &out_PROF1_AcqIndex);
  t_DumpTree->Branch("PROF1_EventIndex",       &out_PROF1_EventIndex);
  t_DumpTree->Branch("PROF1_Multiplicity",     &out_PROF1_Multiplicity);
  t_DumpTree->Branch("PROF1_TriggerTimestamp", &out_PROF1_TriggerTimestamp);
  t_DumpTree->Branch("PROF1_NFibresHit",       &out_PROF1_NFibresHit);
  t_DumpTree->Branch("PROF1_FibresHit",        &out_PROF1_FibresHit);
  t_DumpTree->Branch("PROF2_AcqIndex",         &out_PROF2_AcqIndex);
  t_DumpTree->Branch("PROF2_EventIndex",       &out_PROF2_EventIndex);
  t_DumpTree->Branch("PROF2_Multiplicity",     &out_PROF2_Multiplicity);
  t_DumpTree->Branch("PROF2_TriggerTimestamp", &out_PROF2_TriggerTimestamp);
  t_DumpTree->Branch("PROF2_NFibresHit",       &out_PROF2_NFibresHit);
  t_DumpTree->Branch("PROF2_FibresHit",        &out_PROF2_FibresHit);
  t_DumpTree->Branch("PROF3_AcqIndex",         &out_PROF3_AcqIndex);
  t_DumpTree->Branch("PROF3_EventIndex",       &out_PROF3_EventIndex);
  t_DumpTree->Branch("PROF3_Multiplicity",     &out_PROF3_Multiplicity);
  t_DumpTree->Branch("PROF3_TriggerTimestamp", &out_PROF3_TriggerTimestamp);
  t_DumpTree->Branch("PROF3_NFibresHit",       &out_PROF3_NFibresHit);
  t_DumpTree->Branch("PROF3_FibresHit",        &out_PROF3_FibresHit);
  t_DumpTree->Branch("NFibreCombinations",     &out_NFibreCombinations);
  t_DumpTree->Branch("IsPROFUnique",           &out_IsPROFUnique);
  t_DumpTree->Branch("Current",                &out_Current);
  t_DumpTree->Branch("CosTheta",               &out_CosTheta);
  t_DumpTree->Branch("Theta",                  &out_Theta);
  t_DumpTree->Branch("Momentum",               &out_Momentum);
  t_DumpTree->Branch("Mass",                   &out_Mass);                   

  for(unsigned int i = 0; i < fTFPROFCoincidencesDegenerate.size(); i++)
  {
    TFCoincidenceRecord::TFCoincidence currentTFCO = fTFPROFCoincidencesDegenerate[i].getTFCoincidence();
    std::vector<PROFCoincidenceRecord::PROFCoincidence> currentPROFCOs = fTFPROFCoincidencesDegenerate[i].getPROFCoincidences();

    AcquisitionXBTF::EventRecordHR ER_US = vec_DetAcq[currentTFCO.fDetUS][currentTFCO.fAcqUS].getDataHR()[currentTFCO.fEventUS];
    AcquisitionXBTF::EventRecordHR ER_DS = vec_DetAcq[currentTFCO.fDetDS][currentTFCO.fAcqDS].getDataHR()[currentTFCO.fEventDS];
    std::pair<unsigned int, unsigned int> multiplicity = cTFCoincidenceRecord.getMultiplicity(currentTFCO);

    out_US_DetName          = map_XBTFIndexToName[currentTFCO.fDetUS];
    out_US_AcqIndex         = currentTFCO.fAcqUS;
    out_US_EventIndex       = currentTFCO.fEventUS;
    out_US_Multiplicity     = multiplicity.first;
    out_US_TriggerTimestamp = ER_US.fSeconds_FullTime;

    out_DS_DetName          = map_XBTFIndexToName[currentTFCO.fDetDS];
    out_DS_AcqIndex         = currentTFCO.fAcqDS;
    out_DS_EventIndex       = currentTFCO.fEventDS;
    out_DS_Multiplicity     = multiplicity.second;
    out_DS_TriggerTimestamp = ER_DS.fSeconds_FullTime;

    out_IsTFUnique = currentTFCO.fUniqueness;
    out_TF          = currentTFCO.fTF;
    out_Multiplicty = fTFPROFCoincidencesDegenerate[i].getMultiplicity(); 

    for(unsigned int j = 0; j < currentPROFCOs.size(); j++)
    {
      PROFCoincidenceRecord::PROFCoincidence currentPROFCO = currentPROFCOs[j];
      AcquisitionXBPF::EventRecordHR ER_PROF1 = vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getDataHR()[currentPROFCO.fEventPROF1];
      AcquisitionXBPF::EventRecordHR ER_PROF2 = vec_DetAcqPROF2[currentPROFCO.fAcqPROF2].getDataHR()[currentPROFCO.fEventPROF2];
      AcquisitionXBPF::EventRecordHR ER_PROF3 = vec_DetAcqPROF3[currentPROFCO.fAcqPROF3].getDataHR()[currentPROFCO.fEventPROF3];

      out_PROF1_AcqIndex.push_back(currentPROFCO.fAcqPROF1); out_PROF1_EventIndex.push_back(currentPROFCO.fEventPROF1); out_PROF1_Multiplicity.push_back(currentPROFCO.fMultiplicity[0]); 
      out_PROF1_TriggerTimestamp.push_back(ER_PROF1.fTriggerTimestamp); out_PROF1_NFibresHit.push_back(ER_PROF1.fFibresList.size());
      std::vector<unsigned int> temp1FibresHit;
      for(unsigned int k = 0; k < ER_PROF1.fFibresList.size(); k++)
      {
        temp1FibresHit.push_back(ER_PROF1.fFibresList[k]);
      }
      out_PROF1_FibresHit.push_back(temp1FibresHit);

      out_PROF2_AcqIndex.push_back(currentPROFCO.fAcqPROF2); out_PROF2_EventIndex.push_back(currentPROFCO.fEventPROF2); out_PROF2_Multiplicity.push_back(currentPROFCO.fMultiplicity[0]); 
      out_PROF2_TriggerTimestamp.push_back(ER_PROF2.fTriggerTimestamp); out_PROF2_NFibresHit.push_back(ER_PROF2.fFibresList.size());
      std::vector<unsigned int> temp2FibresHit;
      for(unsigned int k = 0; k < ER_PROF2.fFibresList.size(); k++)
      {
        temp2FibresHit.push_back(ER_PROF2.fFibresList[k]);
      }
      out_PROF2_FibresHit.push_back(temp2FibresHit);

      out_PROF3_AcqIndex.push_back(currentPROFCO.fAcqPROF3); out_PROF3_EventIndex.push_back(currentPROFCO.fEventPROF3); out_PROF3_Multiplicity.push_back(currentPROFCO.fMultiplicity[0]); 
      out_PROF3_TriggerTimestamp.push_back(ER_PROF3.fTriggerTimestamp); out_PROF3_NFibresHit.push_back(ER_PROF3.fFibresList.size());
      std::vector<unsigned int> temp3FibresHit;
      for(unsigned int k = 0; k < ER_PROF3.fFibresList.size(); k++)
      {
        temp3FibresHit.push_back(ER_PROF3.fFibresList[k]);
      }
      out_PROF3_FibresHit.push_back(temp3FibresHit);

      out_NFibreCombinations.push_back(currentPROFCO.getNFibreCombinations());
      out_IsPROFUnique.push_back(currentPROFCO.fUniqueness);
      out_Current.push_back(vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getCurrent());
      std::vector<double> tempCosTheta;
      std::vector<double> tempTheta;
      std::vector<double> tempMomentum;
      std::vector<double> tempMass;
      for(unsigned int k = 0; k < currentPROFCO.getCosTheta().size(); k++)
      {
        tempCosTheta.push_back(currentPROFCO.getCosTheta()[k]);
        tempTheta   .push_back(currentPROFCO.getTheta()[k]);
        tempMomentum.push_back(currentPROFCO.getMomentum()[k]);
        tempMass.push_back(fTFPROFCoincidencesDegenerate[i].fMass[j][k]);
      }
      out_CosTheta.push_back(tempCosTheta);
      out_Theta   .push_back(tempTheta);
      out_Momentum.push_back(tempMomentum);
      out_Mass.push_back(tempMass);
    }

    t_DumpTree->Fill();
    out_US_DetName = ""; out_US_AcqIndex = 0; out_US_EventIndex = 0; out_US_Multiplicity = 0; out_US_TriggerTimestamp = 0.; 
    out_DS_DetName = ""; out_DS_AcqIndex = 0; out_DS_EventIndex = 0; out_DS_Multiplicity = 0; out_DS_TriggerTimestamp = 0.; 
    out_IsTFUnique = 0;
    out_TF = 0.; out_Multiplicty = 0;
    out_PROF1_AcqIndex.clear(); out_PROF1_EventIndex.clear(); out_PROF1_Multiplicity.clear(); out_PROF1_TriggerTimestamp.clear(); out_PROF1_NFibresHit.clear(); out_PROF1_FibresHit.clear();
    out_PROF2_AcqIndex.clear(); out_PROF2_EventIndex.clear(); out_PROF2_Multiplicity.clear(); out_PROF2_TriggerTimestamp.clear(); out_PROF2_NFibresHit.clear(); out_PROF2_FibresHit.clear();
    out_PROF3_AcqIndex.clear(); out_PROF3_EventIndex.clear(); out_PROF3_Multiplicity.clear(); out_PROF3_TriggerTimestamp.clear(); out_PROF3_NFibresHit.clear(); out_PROF3_FibresHit.clear();
    out_IsPROFUnique.clear(); out_Current.clear(); out_NFibreCombinations.clear(); out_CosTheta.clear(); out_Theta.clear(); out_Momentum.clear(); out_Mass.clear();
  }

  t_DumpTree->Write();

  f_DumpFile->Close();
  delete f_DumpFile;
  return;
}

void CombinedCoincidenceRecord::dumpUnique(std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TString const &cFilePathName,
                                           TFCoincidenceRecord &cTFCoincidenceRecord)
{
  gInterpreter->GenerateDictionary("vector<vector<UInt_t> >","vector");

  if(fTFPROFCoincidencesUnique.size()==0)
  {
    getTFPROFCoincidencesUnique();  
  }

  std::map<unsigned int,TString> map_XBTFIndexToName;
  std::vector<std::vector<AcquisitionXBTF>> vec_DetAcq(cMapDetectors.size());
  for(auto det : cMapDetectors)
  {
    if(det.second.getType()=="XBTF")
    {
      det.second.getAcquisitions(vec_DetAcq[cXBTFNameToIndex[det.first]]);     
      map_XBTFIndexToName[cXBTFNameToIndex[det.first]] = det.first;
    }
  }

  std::vector<AcquisitionXBPF> vec_DetAcqPROF1;
  cMapDetectors["XBPF022697"].getAcquisitions(vec_DetAcqPROF1);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF2;
  cMapDetectors["XBPF022701"].getAcquisitions(vec_DetAcqPROF2);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF3;
  cMapDetectors["XBPF022702"].getAcquisitions(vec_DetAcqPROF3);

  TFile *f_DumpFile = new TFile(cFilePathName,"UPDATE");
  TTree *t_DumpTree = new TTree("DumpBL_XBTFXBPF_Unique", "DumpBL_XBTFXBPF_Unique");

  TString out_US_DetName;                                     
  TString out_DS_DetName;
  unsigned int out_US_AcqIndex;
  unsigned int out_US_EventIndex;
  unsigned int out_US_Multiplicity;
  double       out_US_TriggerTimestamp;
  unsigned int out_DS_AcqIndex;
  unsigned int out_DS_EventIndex;
  unsigned int out_DS_Multiplicity;
  double       out_DS_TriggerTimestamp;
  bool         out_IsTFUnique;
  double       out_TF;
  unsigned int out_Multiplicty;
  std::vector<unsigned int> out_PROFCoMultiplicity;
  std::vector<unsigned int> out_PROF1_AcqIndex;               //ACQUISITION INDEX OF PROF1.
  std::vector<unsigned int> out_PROF1_EventIndex;             //EVENT INDEX OF PROF1.
  std::vector<unsigned int> out_PROF1_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF2 AND PROF3 TIMING COINCIDENCES.
  std::vector<long long>    out_PROF1_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF1. NANOSECONDS.
  std::vector<unsigned int> out_PROF1_NFibresHit;             //NUMBER OF PROF1 FIBRES HIT IN THIS EVENT.
  std::vector<std::vector<unsigned int>> out_PROF1_FibresHit; //LIST OF PROF2 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_PROF2_AcqIndex;               //ACQUISITION INDEX OF PROF2.
  std::vector<unsigned int> out_PROF2_EventIndex;             //EVENT INDEX OF PROF2.
  std::vector<unsigned int> out_PROF2_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF1 AND PROF3 TIMING COINCIDENCES.
  std::vector<long long>    out_PROF2_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF2. NANOSECONDS.
  std::vector<unsigned int> out_PROF2_NFibresHit;             //NUMBER OF PROF2 FIBRES HIT IN THIS EVENT.
  std::vector<std::vector<unsigned int>> out_PROF2_FibresHit; //LIST OF PROF2 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_PROF3_AcqIndex;               //ACQUISITION INDEX OF PROF3. 
  std::vector<unsigned int> out_PROF3_EventIndex;             //EVENT INDEX OF PROF3.
  std::vector<unsigned int> out_PROF3_Multiplicity;           //NUMBER OF TIMES THIS EXACT EVENT IS MATCHED WITH OTHER POSSIBLE PROF1 AND PROF2 TIMING COINCIDENCES.
  std::vector<long long>    out_PROF3_TriggerTimestamp;       //TRIGGER TIMESTAMP OF EVENT IN PROF3. NANOSECONDS.
  std::vector<unsigned int> out_PROF3_NFibresHit;             //NUMBER OF PROF3 FIBRES HIT IN THIS EVENT.
  std::vector<std::vector<unsigned int>> out_PROF3_FibresHit; //LIST OF PROF3 FIBRES HIT IN THIS EVENT.
  std::vector<unsigned int> out_NFibreCombinations;           //NUMBER OF POSSIBLE FIBRES COMBINATIONS ACROSS THE THREE PROFILERS.
  std::vector<bool>         out_IsPROFUnique;                 //CURRENT OF THE BENDING MAGNET USED IN THE CALCULATION OF BFIELD, AMPS.
  std::vector<double>       out_Current;                      //FLAG AS TO WHEN THIS COINCIDENCE IS TIME UNIQUE.
  std::vector<std::vector<double>> out_CosTheta;              //CALCULATED COSTHETA.
  std::vector<std::vector<double>> out_Theta;                 //CALCULATED THETA, RADIANS.
  std::vector<std::vector<double>> out_Momentum;              //CALCULATED MOMENTUM, GeV.
  std::vector<std::vector<double>> out_Mass;                  //CALCULATED MASS.

  t_DumpTree->Branch("US_DetName",             &out_US_DetName);           
  t_DumpTree->Branch("DS_DetName",             &out_DS_DetName);           
  t_DumpTree->Branch("US_AcqIndex",            &out_US_AcqIndex);          
  t_DumpTree->Branch("US_EventIndex",          &out_US_EventIndex);        
  t_DumpTree->Branch("US_Multiplicity",        &out_US_Multiplicity);      
  t_DumpTree->Branch("US_TriggerTimestamp",    &out_US_TriggerTimestamp);  
  t_DumpTree->Branch("DS_AcqIndex",            &out_DS_AcqIndex);          
  t_DumpTree->Branch("DS_EventIndex",          &out_DS_EventIndex);        
  t_DumpTree->Branch("DS_Multiplicity",        &out_DS_Multiplicity);      
  t_DumpTree->Branch("DS_TriggerTimestamp",    &out_DS_TriggerTimestamp);  
  t_DumpTree->Branch("IsTFUnique",             &out_IsTFUnique);             
  t_DumpTree->Branch("TF",                     &out_TF);                   
  t_DumpTree->Branch("Multiplicity",           &out_Multiplicty);
  t_DumpTree->Branch("PROF1_AcqIndex",         &out_PROF1_AcqIndex);
  t_DumpTree->Branch("PROF1_EventIndex",       &out_PROF1_EventIndex);
  t_DumpTree->Branch("PROF1_Multiplicity",     &out_PROF1_Multiplicity);
  t_DumpTree->Branch("PROF1_TriggerTimestamp", &out_PROF1_TriggerTimestamp);
  t_DumpTree->Branch("PROF1_NFibresHit",       &out_PROF1_NFibresHit);
  t_DumpTree->Branch("PROF1_FibresHit",        &out_PROF1_FibresHit);
  t_DumpTree->Branch("PROF2_AcqIndex",         &out_PROF2_AcqIndex);
  t_DumpTree->Branch("PROF2_EventIndex",       &out_PROF2_EventIndex);
  t_DumpTree->Branch("PROF2_Multiplicity",     &out_PROF2_Multiplicity);
  t_DumpTree->Branch("PROF2_TriggerTimestamp", &out_PROF2_TriggerTimestamp);
  t_DumpTree->Branch("PROF2_NFibresHit",       &out_PROF2_NFibresHit);
  t_DumpTree->Branch("PROF2_FibresHit",        &out_PROF2_FibresHit);
  t_DumpTree->Branch("PROF3_AcqIndex",         &out_PROF3_AcqIndex);
  t_DumpTree->Branch("PROF3_EventIndex",       &out_PROF3_EventIndex);
  t_DumpTree->Branch("PROF3_Multiplicity",     &out_PROF3_Multiplicity);
  t_DumpTree->Branch("PROF3_TriggerTimestamp", &out_PROF3_TriggerTimestamp);
  t_DumpTree->Branch("PROF3_NFibresHit",       &out_PROF3_NFibresHit);
  t_DumpTree->Branch("PROF3_FibresHit",        &out_PROF3_FibresHit);
  t_DumpTree->Branch("NFibreCombinations",     &out_NFibreCombinations);
  t_DumpTree->Branch("IsPROFUnique",           &out_IsPROFUnique);
  t_DumpTree->Branch("Current",                &out_Current);
  t_DumpTree->Branch("CosTheta",               &out_CosTheta);
  t_DumpTree->Branch("Theta",                  &out_Theta);
  t_DumpTree->Branch("Momentum",               &out_Momentum);
  t_DumpTree->Branch("Mass",                   &out_Mass);                   

  for(unsigned int i = 0; i < fTFPROFCoincidencesUnique.size(); i++)
  {
    TFCoincidenceRecord::TFCoincidence currentTFCO = fTFPROFCoincidencesUnique[i].getTFCoincidence();
    std::vector<PROFCoincidenceRecord::PROFCoincidence> currentPROFCOs = fTFPROFCoincidencesUnique[i].getPROFCoincidences();

    AcquisitionXBTF::EventRecordHR ER_US = vec_DetAcq[currentTFCO.fDetUS][currentTFCO.fAcqUS].getDataHR()[currentTFCO.fEventUS];
    AcquisitionXBTF::EventRecordHR ER_DS = vec_DetAcq[currentTFCO.fDetDS][currentTFCO.fAcqDS].getDataHR()[currentTFCO.fEventDS];
    std::pair<unsigned int, unsigned int> multiplicity = cTFCoincidenceRecord.getMultiplicity(currentTFCO);

    out_US_DetName          = map_XBTFIndexToName[currentTFCO.fDetUS];
    out_US_AcqIndex         = currentTFCO.fAcqUS;
    out_US_EventIndex       = currentTFCO.fEventUS;
    out_US_Multiplicity     = multiplicity.first;
    out_US_TriggerTimestamp = ER_US.fSeconds_FullTime;

    out_DS_DetName          = map_XBTFIndexToName[currentTFCO.fDetDS];
    out_DS_AcqIndex         = currentTFCO.fAcqDS;
    out_DS_EventIndex       = currentTFCO.fEventDS;
    out_DS_Multiplicity     = multiplicity.second;
    out_DS_TriggerTimestamp = ER_DS.fSeconds_FullTime;

    out_IsTFUnique  = currentTFCO.fUniqueness;
    out_TF          = currentTFCO.fTF;
    out_Multiplicty = fTFPROFCoincidencesUnique[i].getMultiplicity();

    for(unsigned int j = 0; j < currentPROFCOs.size(); j++)
    {
      PROFCoincidenceRecord::PROFCoincidence currentPROFCO = currentPROFCOs[j];
      AcquisitionXBPF::EventRecordHR ER_PROF1 = vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getDataHR()[currentPROFCO.fEventPROF1];
      AcquisitionXBPF::EventRecordHR ER_PROF2 = vec_DetAcqPROF2[currentPROFCO.fAcqPROF2].getDataHR()[currentPROFCO.fEventPROF2];
      AcquisitionXBPF::EventRecordHR ER_PROF3 = vec_DetAcqPROF3[currentPROFCO.fAcqPROF3].getDataHR()[currentPROFCO.fEventPROF3];

      out_PROF1_AcqIndex.push_back(currentPROFCO.fAcqPROF1); out_PROF1_EventIndex.push_back(currentPROFCO.fEventPROF1); out_PROF1_Multiplicity.push_back(currentPROFCO.fMultiplicity[0]); 
      out_PROF1_TriggerTimestamp.push_back(ER_PROF1.fTriggerTimestamp); out_PROF1_NFibresHit.push_back(ER_PROF1.fFibresList.size());
      std::vector<unsigned int> temp1FibresHit;
      for(unsigned int k = 0; k < ER_PROF1.fFibresList.size(); k++)
      {
        temp1FibresHit.push_back(ER_PROF1.fFibresList[k]);
      }
      out_PROF1_FibresHit.push_back(temp1FibresHit);

      out_PROF2_AcqIndex.push_back(currentPROFCO.fAcqPROF2); out_PROF2_EventIndex.push_back(currentPROFCO.fEventPROF2); out_PROF2_Multiplicity.push_back(currentPROFCO.fMultiplicity[0]); 
      out_PROF2_TriggerTimestamp.push_back(ER_PROF2.fTriggerTimestamp); out_PROF2_NFibresHit.push_back(ER_PROF2.fFibresList.size());
      std::vector<unsigned int> temp2FibresHit;
      for(unsigned int k = 0; k < ER_PROF2.fFibresList.size(); k++)
      {
        temp2FibresHit.push_back(ER_PROF2.fFibresList[k]);
      }
      out_PROF2_FibresHit.push_back(temp2FibresHit);

      out_PROF3_AcqIndex.push_back(currentPROFCO.fAcqPROF3); out_PROF3_EventIndex.push_back(currentPROFCO.fEventPROF3); out_PROF3_Multiplicity.push_back(currentPROFCO.fMultiplicity[0]); 
      out_PROF3_TriggerTimestamp.push_back(ER_PROF3.fTriggerTimestamp); out_PROF3_NFibresHit.push_back(ER_PROF3.fFibresList.size());
      std::vector<unsigned int> temp3FibresHit;
      for(unsigned int k = 0; k < ER_PROF3.fFibresList.size(); k++)
      {
        temp3FibresHit.push_back(ER_PROF3.fFibresList[k]);
      }
      out_PROF3_FibresHit.push_back(temp3FibresHit);

      out_NFibreCombinations.push_back(currentPROFCO.getNFibreCombinations());
      out_IsPROFUnique.push_back(currentPROFCO.fUniqueness);
      out_Current.push_back(vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getCurrent());
      std::vector<double> tempCosTheta;
      std::vector<double> tempTheta;
      std::vector<double> tempMomentum;
      std::vector<double> tempMass;
      for(unsigned int k = 0; k < currentPROFCO.getCosTheta().size(); k++)
      {
        tempCosTheta.push_back(currentPROFCO.getCosTheta()[k]);
        tempTheta   .push_back(currentPROFCO.getTheta()[k]);
        tempMomentum.push_back(currentPROFCO.getMomentum()[k]);
        tempMass.push_back(fTFPROFCoincidencesUnique[i].fMass[j][k]);
      }
      out_CosTheta.push_back(tempCosTheta);
      out_Theta   .push_back(tempTheta);
      out_Momentum.push_back(tempMomentum);
      out_Mass.push_back(tempMass);
    }

    t_DumpTree->Fill();
    out_US_DetName = ""; out_US_AcqIndex = 0; out_US_EventIndex = 0; out_US_Multiplicity = 0; out_US_TriggerTimestamp = 0.; 
    out_DS_DetName = ""; out_DS_AcqIndex = 0; out_DS_EventIndex = 0; out_DS_Multiplicity = 0; out_DS_TriggerTimestamp = 0.; 
    out_IsTFUnique = 0;
    out_TF = 0.; out_Multiplicty = 0;
    out_PROF1_AcqIndex.clear(); out_PROF1_EventIndex.clear(); out_PROF1_Multiplicity.clear(); out_PROF1_TriggerTimestamp.clear(); out_PROF1_NFibresHit.clear(); out_PROF1_FibresHit.clear();
    out_PROF2_AcqIndex.clear(); out_PROF2_EventIndex.clear(); out_PROF2_Multiplicity.clear(); out_PROF2_TriggerTimestamp.clear(); out_PROF2_NFibresHit.clear(); out_PROF2_FibresHit.clear();
    out_PROF3_AcqIndex.clear(); out_PROF3_EventIndex.clear(); out_PROF3_Multiplicity.clear(); out_PROF3_TriggerTimestamp.clear(); out_PROF3_NFibresHit.clear(); out_PROF3_FibresHit.clear();
    out_IsPROFUnique.clear(); out_Current.clear(); out_NFibreCombinations.clear(); out_CosTheta.clear(); out_Theta.clear(); out_Momentum.clear(); out_Mass.clear();
  }

  t_DumpTree->Write();

  f_DumpFile->Close();
  delete f_DumpFile;

  return;
}

void CombinedCoincidenceRecord::printDegenerate(std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex,
                                                TFCoincidenceRecord &cTFCoincidenceRecord)
{
  /*
  if(fTFPROFCoincidencesUnique.size()==0)
  {
    getTFPROFCoincidencesUnique();  
  }*/

  std::vector<std::vector<AcquisitionXBTF>> vec_DetAcq(cMapDetectors.size());
  for(auto det : cMapDetectors)
  {
    if(det.second.getType()=="XBTF")
    {
      det.second.getAcquisitions(vec_DetAcq[cXBTFNameToIndex[det.first]]);     
    }
  }

  std::vector<AcquisitionXBPF> vec_DetAcqPROF1;
  cMapDetectors["XBPF022697"].getAcquisitions(vec_DetAcqPROF1);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF2;
  cMapDetectors["XBPF022701"].getAcquisitions(vec_DetAcqPROF2);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF3;
  cMapDetectors["XBPF022702"].getAcquisitions(vec_DetAcqPROF3);

  std::cout << "**************************************************************************" << std::endl;
  std::cout << "PRINTING DEGENERATE TF PROF COINCIDENCES" << std::endl;

  for(unsigned int i = 0; i < fTFPROFCoincidencesDegenerate.size(); i++)
  {
    TFCoincidenceRecord::TFCoincidence currentTFCO = fTFPROFCoincidencesDegenerate[i].getTFCoincidence();
    std::vector<PROFCoincidenceRecord::PROFCoincidence> currentPROFCOs = fTFPROFCoincidencesDegenerate[i].getPROFCoincidences();

    AcquisitionXBTF::EventRecordHR ER_US = vec_DetAcq[currentTFCO.fDetUS][currentTFCO.fAcqUS].getDataHR()[currentTFCO.fEventUS];
    AcquisitionXBTF::EventRecordHR ER_DS = vec_DetAcq[currentTFCO.fDetDS][currentTFCO.fAcqDS].getDataHR()[currentTFCO.fEventDS];
    std::pair<unsigned int, unsigned int> multiplicity = cTFCoincidenceRecord.getMultiplicity(currentTFCO);

    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "-----------------MATCHED THE FOLLOWING TF COINCIDENCE--------------------" << std::endl;
    std::cout << "UPSTREAM DET:  " << currentTFCO.fDetUS << ", ACQ " << currentTFCO.fAcqUS << ", EVENT " << currentTFCO.fEventUS << " TIME: " << ER_US.fSeconds_FullTime << std::endl;
    std::cout << "DOWNSTREM DET: " << currentTFCO.fDetDS << ", ACQ " << currentTFCO.fAcqDS << ", EVENT " << currentTFCO.fEventDS << " TIME: " << ER_DS.fSeconds_FullTime << std::endl;
    std::cout << "UPSTREAM MULTIPLICITY: " << multiplicity.first << ", DOWNSTREAM MULTIPLICITY: " << multiplicity.second << ", TF: " << currentTFCO.fTF << std::endl;

    std::cout << "---------------WITH THE FOLLOWING PROF COINCIDENCES-------------------" << std::endl;
    std::cout << "    THE MULTIPLICITY IS: " << fTFPROFCoincidencesDegenerate[i].getMultiplicity() << std::endl;

    for(unsigned int j = 0; j < currentPROFCOs.size(); j++)
    {
      PROFCoincidenceRecord::PROFCoincidence currentPROFCO = currentPROFCOs[j];
      AcquisitionXBPF::EventRecordHR ER_PROF1 = vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getDataHR()[currentPROFCO.fEventPROF1];
      AcquisitionXBPF::EventRecordHR ER_PROF2 = vec_DetAcqPROF2[currentPROFCO.fAcqPROF2].getDataHR()[currentPROFCO.fEventPROF2];
      AcquisitionXBPF::EventRecordHR ER_PROF3 = vec_DetAcqPROF3[currentPROFCO.fAcqPROF3].getDataHR()[currentPROFCO.fEventPROF3];

      std::cout << "PROF1:   ACQ " << currentPROFCO.fAcqPROF1 << ", EVENT " << currentPROFCO.fEventPROF1 << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[0] <<  " TIME: " << ER_PROF1.fTriggerTimestamp
        << "\n    FIBRES ";
      for(unsigned int k = 0; k < ER_PROF1.fFibresList.size(); k++)
      {
        std::cout << ER_PROF1.fFibresList[k] << " ";
      }
      std::cout << std::endl;
      std::cout << "PROF2:   ACQ " << currentPROFCO.fAcqPROF2 << ", EVENT " << currentPROFCO.fEventPROF2  << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[1] << " TIME: " << ER_PROF2.fTriggerTimestamp
        << "\n    FIBRES ";
      for(unsigned int k = 0; k < ER_PROF2.fFibresList.size(); k++)
      {
        std::cout << ER_PROF2.fFibresList[k] << " ";
      }
      std::cout << std::endl;
      std::cout << "PROF3:   ACQ " << currentPROFCO.fAcqPROF3 << ", EVENT " << currentPROFCO.fEventPROF3  << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[2] << " TIME: " << ER_PROF3.fTriggerTimestamp
        << "\n    FIBRES ";
      for(unsigned int k = 0; k < ER_PROF3.fFibresList.size(); k++)
      {
        std::cout << ER_PROF3.fFibresList[k] << " ";
      }
      std::cout << std::endl;
      std::cout << "    NUMBER OF POSSIBLE FIBRE COMBINATIONS: " << currentPROFCO.getNFibreCombinations() << std::endl;
      for(unsigned int k = 0; k < currentPROFCO.getCosTheta().size(); k++)
      {
        std::cout << "    COSTHETA: " << currentPROFCO.getCosTheta()[k] << ", THETA (DEGREES): " << currentPROFCO.getTheta()[k]*radToDeg << ", MOMENTUM (GeV): " << currentPROFCO.getMomentum()[k] << std::endl;
        std::cout << "    MASS (GeV): " << fTFPROFCoincidencesDegenerate[i].fMass[j][k] << std::endl; 
      }
      std::cout << "-------------------------------------------------------------------------" << std::endl;
    }
  }
  std::cout << "**************************************************************************" << std::endl;

  return;
}

void CombinedCoincidenceRecord::printUnique(std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TFCoincidenceRecord &cTFCoincidenceRecord)
{
  if(fTFPROFCoincidencesUnique.size()==0)
  {
    getTFPROFCoincidencesUnique();  
  }

  std::vector<std::vector<AcquisitionXBTF>> vec_DetAcq(cMapDetectors.size());
  for(auto det : cMapDetectors)
  {
    if(det.second.getType()=="XBTF")
    {
      det.second.getAcquisitions(vec_DetAcq[cXBTFNameToIndex[det.first]]);     
    }
  }

  std::vector<AcquisitionXBPF> vec_DetAcqPROF1;
  cMapDetectors["XBPF022697"].getAcquisitions(vec_DetAcqPROF1);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF2;
  cMapDetectors["XBPF022701"].getAcquisitions(vec_DetAcqPROF2);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF3;
  cMapDetectors["XBPF022702"].getAcquisitions(vec_DetAcqPROF3);

  std::cout << "**************************************************************************" << std::endl;
  std::cout << "PRINTING DEGENERATE TF PROF COINCIDENCES" << std::endl;

  for(unsigned int i = 0; i < 10; i++)
  //for(unsigned int i = 0; i < fTFPROFCoincidencesUnique.size(); i++)
  {
    TFCoincidenceRecord::TFCoincidence currentTFCO = fTFPROFCoincidencesUnique[i].getTFCoincidence();
    std::vector<PROFCoincidenceRecord::PROFCoincidence> currentPROFCOs = fTFPROFCoincidencesDegenerate[i].getPROFCoincidences();

    AcquisitionXBTF::EventRecordHR ER_US = vec_DetAcq[currentTFCO.fDetUS][currentTFCO.fAcqUS].getDataHR()[currentTFCO.fEventUS];
    AcquisitionXBTF::EventRecordHR ER_DS = vec_DetAcq[currentTFCO.fDetDS][currentTFCO.fAcqDS].getDataHR()[currentTFCO.fEventDS];
    std::pair<unsigned int, unsigned int> multiplicity = cTFCoincidenceRecord.getMultiplicity(currentTFCO);

    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "-----------------MATCHED THE FOLLOWING TF COINCIDENCE--------------------" << std::endl;
    std::cout << "UPSTREAM DET:  " << currentTFCO.fDetUS << ", ACQ " << currentTFCO.fAcqUS << ", EVENT " << currentTFCO.fEventUS << " TIME: " << ER_US.fSeconds_FullTime << std::endl;
    std::cout << "DOWNSTREM DET: " << currentTFCO.fDetDS << ", ACQ " << currentTFCO.fAcqDS << ", EVENT " << currentTFCO.fEventDS << " TIME: " << ER_DS.fSeconds_FullTime << std::endl;
    std::cout << "UPSTREAM MULTIPLICITY: " << multiplicity.first << ", DOWNSTREAM MULTIPLICITY: " << multiplicity.second << ", TF: " << currentTFCO.fTF << std::endl;

    std::cout << "---------------WITH THE FOLLOWING PROF COINCIDENCES-------------------" << std::endl;
    std::cout << "    THE MULTIPLICITY IS: " << fTFPROFCoincidencesUnique[i].getMultiplicity() << std::endl;

    for(unsigned int j = 0; j < currentPROFCOs.size(); j++)
    {
      PROFCoincidenceRecord::PROFCoincidence currentPROFCO = currentPROFCOs[j];
      AcquisitionXBPF::EventRecordHR ER_PROF1 = vec_DetAcqPROF1[currentPROFCO.fAcqPROF1].getDataHR()[currentPROFCO.fEventPROF1];
      AcquisitionXBPF::EventRecordHR ER_PROF2 = vec_DetAcqPROF2[currentPROFCO.fAcqPROF2].getDataHR()[currentPROFCO.fEventPROF2];
      AcquisitionXBPF::EventRecordHR ER_PROF3 = vec_DetAcqPROF3[currentPROFCO.fAcqPROF3].getDataHR()[currentPROFCO.fEventPROF3];

      std::cout << "-------------------------------------------------------------------------" << std::endl;
      std::cout << "PROF1:   ACQ " << currentPROFCO.fAcqPROF1 << ", EVENT " << currentPROFCO.fEventPROF1 << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[0] <<  " TIME: " << ER_PROF1.fTriggerTimestamp
        << "\n    FIBRES ";
      for(unsigned int k = 0; k < ER_PROF1.fFibresList.size(); k++)
      {
        std::cout << ER_PROF1.fFibresList[k] << " ";
      }
      std::cout << std::endl;
      std::cout << "PROF2:   ACQ " << currentPROFCO.fAcqPROF2 << ", EVENT " << currentPROFCO.fEventPROF2  << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[1] << " TIME: " << ER_PROF2.fTriggerTimestamp
        << "\n    FIBRES ";
      for(unsigned int k = 0; k < ER_PROF2.fFibresList.size(); k++)
      {
        std::cout << ER_PROF2.fFibresList[k] << " ";
      }
      std::cout << std::endl;
      std::cout << "PROF3:   ACQ " << currentPROFCO.fAcqPROF3 << ", EVENT " << currentPROFCO.fEventPROF3  << ", MULTIPLICITY " << currentPROFCO.fMultiplicity[2] << " TIME: " << ER_PROF3.fTriggerTimestamp
        << "\n    FIBRES ";
      for(unsigned int k = 0; k < ER_PROF3.fFibresList.size(); k++)
      {
        std::cout << ER_PROF3.fFibresList[k] << " ";
      }
      std::cout << std::endl;
      std::cout << "    NUMBER OF POSSIBLE FIBRE COMBINATIONS: " << currentPROFCO.getNFibreCombinations() << std::endl;
      for(unsigned int k = 0; k < currentPROFCO.getCosTheta().size(); k++)
      {
        std::cout << "    COSTHETA: " << currentPROFCO.getCosTheta()[k] << ", THETA (DEGREES): " << currentPROFCO.getTheta()[k]*radToDeg << ", MOMENTUM (GeV): " << currentPROFCO.getMomentum()[k] << std::endl;
        std::cout << "    MASS (GeV): " << fTFPROFCoincidencesUnique[i].fMass[j][k] << std::endl; 
      }
      std::cout << "-------------------------------------------------------------------------" << std::endl;
    }
  }
  std::cout << "**************************************************************************" << std::endl;

  return;
}


