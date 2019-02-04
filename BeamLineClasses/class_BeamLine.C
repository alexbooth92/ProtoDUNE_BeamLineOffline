#include "class_BeamLine.h"
#include "class_Detector.h"
#include <TInterpreter.h>


BeamLine::BeamLine()
{
  return;
}

double BeamLine::getCoordinate(unsigned int const &cDet, unsigned int const &cFibre)
{
  double systShift702 = -1.45;
  if(cDet==1)
  {
    return -1.*FibreToCoord_XBPF022697[cFibre]*1e-3;
  }
  else if(cDet==2)
  {
    return -1.*FibreToCoord_XBPF022701[cFibre]*1e-3;
  }
  else if(cDet==3)
  {
    return -1.*(FibreToCoord_XBPF022702[cFibre]+systShift702)*1e-3;
  }

  return 0.;
}

std::vector<double> BeamLine::getCoordinates(std::string const &cDet, std::vector<unsigned int> const &cFibre)
{
  std::vector<double> fibreCoords;
  for(unsigned int i = 0; i < cFibre.size(); i++)
  {
    //ORIENTATION IE THE MINUS IS CONSISTENT WITH FIBRE 0 ON JURA SIDE AND FIBRE 0 ON TOP.
    if     (cDet.find("697")!=std::string::npos){fibreCoords.push_back(-1.*FibreToCoord_XBPF022697[cFibre[i]]*1e-3);}
    else if(cDet.find("698")!=std::string::npos){fibreCoords.push_back(-1.*FibreToCoord_XBPF022701[cFibre[i]]*1e-3);}
    else if(cDet.find("701")!=std::string::npos){fibreCoords.push_back(-1.*FibreToCoord_XBPF022701[cFibre[i]]*1e-3);}
    else if(cDet.find("702")!=std::string::npos){fibreCoords.push_back(-1.*FibreToCoord_XBPF022702[cFibre[i]]*1e-3);}
    else if(cDet.find("707")!=std::string::npos){fibreCoords.push_back(-1.*FibreToCoord_XBPF022701[cFibre[i]]*1e-3);}
    else if(cDet.find("708")!=std::string::npos){fibreCoords.push_back(-1.*FibreToCoord_XBPF022702[cFibre[i]]*1e-3);}
    else if(cDet.find("716")!=std::string::npos){fibreCoords.push_back(-1.*FibreToCoord_XBPF022701[cFibre[i]]*1e-3);}
    else if(cDet.find("717")!=std::string::npos){fibreCoords.push_back(-1.*FibreToCoord_XBPF022702[cFibre[i]]*1e-3);}
  }

  return fibreCoords;
}


double BeamLine::getCosTheta(unsigned int const &cPROF1, unsigned const &cPROF2, unsigned int const &cPROF3)
{
  double CosTheta0 = std::cos(Theta0);
  double TanTheta0 = std::tan(Theta0);
  double XPROF1    = getCoordinate(1, cPROF1);
  double XPROF2    = getCoordinate(2, cPROF2);
  double XPROF3    = getCoordinate(3, cPROF3);

  /*double a    = (XPROF2*Offset_XBPF022702-XPROF3*Offset_XBPF022701)/(Offset_XBPF022702-Offset_XBPF022701); 
  double num  = ((a-XPROF1)*(XPROF3-XPROF2)*CosTheta0 - (Offset_XBPF022701-Offset_XBPF022702)*(Offset_XBPF022697+(a-XPROF1)*TanTheta0));
  double den1 = std::sqrt(Offset_XBPF022697*Offset_XBPF022697+(a-XPROF1)*(a-XPROF1));
  double den2 = std::sqrt(((Offset_XBPF022701-Offset_XBPF022702)*std::tan(Theta0)-(XPROF3-XPROF2)*std::cos(Theta0))*((Offset_XBPF022701-Offset_XBPF022702)*std::tan(Theta0)-(XPROF3-XPROF2)*std::cos(Theta0))
                          +(Offset_XBPF022701-Offset_XBPF022702)*(Offset_XBPF022701-Offset_XBPF022702));*/
  double a    = (XPROF3*Offset_XBPF022701-XPROF2*Offset_XBPF022702)*CosTheta0/(Offset_XBPF022702-Offset_XBPF022701); 
  double num  = (a+XPROF1)*((Offset_XBPF022702-Offset_XBPF022701)*TanTheta0+(XPROF2-XPROF3)*CosTheta0)+Offset_XBPF022697*(Offset_XBPF022702-Offset_XBPF022701);
  double den1 = std::sqrt(Offset_XBPF022697*Offset_XBPF022697+(a+XPROF1)*(a+XPROF1));
  double den2 = std::sqrt(((Offset_XBPF022702-Offset_XBPF022701)*TanTheta0+(XPROF2-XPROF3)*CosTheta0)*((Offset_XBPF022702-Offset_XBPF022701)*TanTheta0+(XPROF2-XPROF3)*CosTheta0)+
                           (Offset_XBPF022702-Offset_XBPF022701)*(Offset_XBPF022702-Offset_XBPF022701));

  double cosTheta = num/(den1*den2);

  return cosTheta;
}

double BeamLine::getMomentum(double const &cCosTheta, double const &cCurrent)
{
  double LB     = mag_P1*cCurrent;
  double deltaI = cCurrent - mag_P4; 
  if(deltaI>0)
  {
    LB+= mag_P3*deltaI*deltaI;
  }

  return (299792458*LB/(1.e9*(std::acos(cCosTheta)))); 
}

void BeamLine::considerMomenta(std::vector<unsigned int> const &cPROF1Fibs, std::vector<unsigned int> const &cPROF2Fibs, std::vector<unsigned int> const &cPROF3Fibs, double const &cCurrent,
                               std::vector<double> &cCosTheta, std::vector<double> &cTheta, std::vector<double> &cMomentum)
{
  for(unsigned int i = 0; i < cPROF1Fibs.size(); i++)
  {
    for(unsigned int j = 0; j < cPROF2Fibs.size(); j++)
    {
      for(unsigned int k = 0; k < cPROF3Fibs.size(); k++)
      {
        double cosTheta = getCosTheta(cPROF1Fibs[i], cPROF2Fibs[j], cPROF3Fibs[k]);
        cCosTheta.push_back(cosTheta);
        cTheta.push_back(std::acos(cosTheta));
        cMomentum.push_back(getMomentum(cosTheta, cCurrent));
      }
    }
  }

  return;
}

std::vector<TFCoincidenceRecord::TFCoincidence> BeamLine::getTFCoincidencesDegenerate()
{
  return fTFCoincidenceRecord.getTFCoincidencesDegenerate();
}

std::vector<TFCoincidenceRecord::TFCoincidence> BeamLine::getTFCoincidencesUnique()
{
  return fTFCoincidenceRecord.getTFCoincidencesUnique();
}

std::vector<PROFCoincidenceRecord::PROFCoincidence> BeamLine::getPROFCoincidencesDegenerate()
{
  return fPROFCoincdenceRecord.getPROFCoincidencesDegenerate();
}

std::vector<PROFCoincidenceRecord::PROFCoincidence> BeamLine::getPROFCoincidencesUnique()
{
  return fPROFCoincdenceRecord.getPROFCoincidencesUnique();
}

std::vector<CombinedCoincidenceRecord::TFPROFCoincidence> BeamLine::getTFPROFCoincidencesDegenerate()
{
  return fCombinedCoincidenceRecord.getTFPROFCoincidencesDegenerate();
}

std::vector<CombinedCoincidenceRecord::TFPROFCoincidence> BeamLine::getTFPROFCoincidencesUnique()
{
  return fCombinedCoincidenceRecord.getTFPROFCoincidencesUnique();
}

void BeamLine::findTFCoincidences(std::map<std::string,Detector> &cMapDetectors)
{
  std::vector<Detector> vec_TFDet;
  unsigned int index = 0;
  for(auto detector : cMapDetectors)
  {
    //if(!detector.first.find("XBTF"))
    if(detector.second.getType()=="XBTF")
    {
      vec_TFDet.push_back(detector.second);
      fXBTFDetNameToIndex[detector.first] = index;
      fXBTFIndexToDetName[index] = detector.first;
      index++;
    }
  }

  //LOOP I OVER ALL POSSIBLE US DETECTORS, LOOP J OVER ALL POSSIBLE DS DETECTORS.
  for(unsigned int i = 0; i < vec_TFDet.size(); i++)
  {
    for(unsigned int j = 0; j < vec_TFDet.size(); j++)
    {
      std::cout << fXBTFIndexToDetName[i] << std::endl;
      //ONLY COMPARE ONE US AND DS PAIR ONCE.
      //if(i < j)
      if((fXBTFIndexToDetName[i]=="XBTF022687A" && fXBTFIndexToDetName[j]=="XBTF022716A")||
         (fXBTFIndexToDetName[i]=="XBTF022687A" && fXBTFIndexToDetName[j]=="XBTF022716B")||
         (fXBTFIndexToDetName[i]=="XBTF022687B" && fXBTFIndexToDetName[j]=="XBTF022716A")||
         (fXBTFIndexToDetName[i]=="XBTF022687B" && fXBTFIndexToDetName[j]=="XBTF022716B"))
      {
        std::cout << "CALCULATING TF BETWEEN: " << vec_TFDet[i].getDetName() << " AND " << vec_TFDet[j].getDetName() << std::endl;

        std::vector<AcquisitionXBTF> vec_DetAcqUS;
        vec_TFDet[i].getAcquisitions(vec_DetAcqUS);
        std::vector<AcquisitionXBTF> vec_DetAcqDS;
        vec_TFDet[j].getAcquisitions(vec_DetAcqDS);

        //LOOP K OVER US TF ACQUISITIONS.
        for(unsigned int k = 0; k < vec_DetAcqUS.size(); k++)
        {
          //LOOP L OVER DS TF ACQUISITIONS.
          for(unsigned int l = 0; l < vec_DetAcqDS.size(); l++)
          {
            std::vector<AcquisitionXBTF::EventRecordHR> vec_CurrentRecordsUS = vec_DetAcqUS[k].getDataHR(); 
            std::vector<AcquisitionXBTF::EventRecordHR> vec_CurrentRecordsDS = vec_DetAcqDS[l].getDataHR(); 
            if(vec_CurrentRecordsDS.size()!=0 && vec_CurrentRecordsUS.size()!=0)
            {
              //IF THE EARLIEST UPSTREAM TIME IS LATER THAN THE LATEST DOWNSTREAM TIME, DONT ACCEPT.
              if(!((((double)vec_CurrentRecordsUS[0].fSeconds   -(double)vec_CurrentRecordsDS[vec_CurrentRecordsDS.size()-1].fSeconds)
                           +(vec_CurrentRecordsUS[0].fSubSeconds-vec_CurrentRecordsDS[vec_CurrentRecordsDS.size()-1].fSubSeconds)/1.e9)*1.e9>0))
              {
                //IF THE LATERST UPSTREAM TIME IS MORE THAN THE TOLERANCE AWAY FROM THE FIRST DS TIME, DONT ACCEPT.
                if(!((((double)vec_CurrentRecordsUS[vec_CurrentRecordsUS.size()-1].fSeconds-(double)vec_CurrentRecordsDS[0].fSeconds)
                             +(vec_CurrentRecordsUS[vec_CurrentRecordsUS.size()-1].fSubSeconds     -vec_CurrentRecordsDS[0].fSubSeconds)/1.e9)*1.e9<-1.*toleranceTF))
                {
                  unsigned int newStartValue = 0;
                  //LOOP OVER THE EVENTS IN US ACQUISITION K.
                  for(unsigned int m = 0; m < (vec_CurrentRecordsUS.size() < eventLimit ? vec_CurrentRecordsUS.size() : eventLimit); m++)
                  //for(unsigned int m = 0; m < vec_CurrentRecordsUS.size(); m++)
                  {
                    bool beenInRegion = false;
                    if(!((((double)vec_CurrentRecordsUS[m].fSeconds   -(double)vec_CurrentRecordsDS[vec_CurrentRecordsDS.size()-1].fSeconds)
                                 +(vec_CurrentRecordsUS[m].fSubSeconds-vec_CurrentRecordsDS[vec_CurrentRecordsDS.size()-1].fSubSeconds)/1.e9)*1.e9>0))
                    {
                      //LOOP OVER THE EVENTS IN DS ACQUISITION L.
                      for(unsigned int n = newStartValue; n < (vec_CurrentRecordsDS.size() < eventLimit ? vec_CurrentRecordsDS.size() : eventLimit); n++)
                      //for(unsigned int n = newStartValue; n < vec_CurrentRecordsDS.size(); n++)
                      {
                        double deltaT = (((double)vec_CurrentRecordsUS[m].fSeconds   -(double)vec_CurrentRecordsDS[n].fSeconds)
                                                +(vec_CurrentRecordsUS[m].fSubSeconds-vec_CurrentRecordsDS[n].fSubSeconds)/1.e9)*1.e9;
                        if(deltaT<-1.*toleranceTF)
                        {
                          break;  
                        }
                        else if(deltaT>=(-1*toleranceTF) && deltaT<=0)
                        {
                          fTFCoincidenceRecord.addCoincidence(i, k, m, j, l, n, -1*deltaT);
                          if(beenInRegion==false)
                          {
                            beenInRegion  = true;
                            if(n>0)
                            {
                              newStartValue = n - 1; 
                            }
                            else
                            {
                              newStartValue = 0;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return;
}

void BeamLine::findPROFCoincidences(std::map<std::string,Detector> &cMapDetectors) 
{
  std::map<std::string,Detector> map_PFDet;
  for(auto detector : cMapDetectors)
  {
    std::string s_detName = detector.second.getDetName();
    if(s_detName=="XBPF022697"|| s_detName=="XBPF022701" || s_detName=="XBPF022702")
    {
      map_PFDet[s_detName] = detector.second;
    }
  }
  if(map_PFDet.size()<3)
  {
    std::cout << "NOT ENOUGH XBPFs TO MAKE DEFLECTION CALCULATION" << std::endl;
  }

  std::vector<AcquisitionXBPF> vec_DetAcqPROF1;
  map_PFDet["XBPF022697"].getAcquisitions(vec_DetAcqPROF1);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF2;
  map_PFDet["XBPF022701"].getAcquisitions(vec_DetAcqPROF2);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF3;
  map_PFDet["XBPF022702"].getAcquisitions(vec_DetAcqPROF3);

  //LOOPS I AND J RUN OVER THE ACQUISITIONS RECORDED FOR EACH DETECTOR SINCE THE LAST TIME THE SCRIPT WAS RUN.
  for(unsigned int i = 0; i < vec_DetAcqPROF1.size(); i++)
  {
    for(unsigned int j = 0; j < vec_DetAcqPROF3.size(); j++)
    {
      //GET THE UNPACKED EVENTS FROM THE I AND J ACQUISITIONS.
      std::vector<AcquisitionXBPF::EventRecordHR> vec_CurrentRecordsPROF1 = vec_DetAcqPROF1[i].getDataHR(); 
      std::vector<AcquisitionXBPF::EventRecordHR> vec_CurrentRecordsPROF3 = vec_DetAcqPROF3[j].getDataHR(); 
      //LOOP K OVER THE EVENTS IN ACQUISITION I FROM PROF1.
      if(vec_DetAcqPROF1[i].getNNonZeroEvents()!=0 && vec_DetAcqPROF3[j].getNNonZeroEvents()!=0)
      //if(vec_CurrentRecordsPROF1.size()!=0 && vec_CurrentRecordsPROF3.size()!=0)
      {
        //IF THE EARLIEST UPSTREAM TIME IS LATER THAN THE LATEST DOWNSTREAM TIME, DONT ACCEPT.
        if(!(vec_CurrentRecordsPROF1[0].fTriggerTimestamp-vec_CurrentRecordsPROF3[vec_DetAcqPROF3[j].getNNonZeroEvents()-1].fTriggerTimestamp>0))
        {
          //IF THE LATERST UPSTREAM TIME IS MORE THAN THE TOLERANCE AWAY FROM THE FIRST DS TIME, DONT ACCEPT.
          if(!(vec_CurrentRecordsPROF1[vec_DetAcqPROF1[i].getNNonZeroEvents()-1].fTriggerTimestamp-vec_CurrentRecordsPROF3[0].fTriggerTimestamp<-1*tolerancePROF))
          {
            unsigned int newStartValue = 0;
            for(unsigned int k = 0; k < (vec_DetAcqPROF1[i].getNNonZeroEvents()<eventLimit ? vec_DetAcqPROF1[i].getNNonZeroEvents() : eventLimit); k++)
            //for(unsigned int k = 0; k < vec_DetAcqPROF1[i].getNNonZeroEvents(); k++)
            {
              bool beenInRegion = false;
              //IF THE CURRENT UPSTREAM TIME IS LATER THAN THE LATEST DOWNSTREAM TIME, DONT ACCEPT.
              if(!(vec_CurrentRecordsPROF1[k].fTriggerTimestamp-vec_CurrentRecordsPROF3[vec_DetAcqPROF3[j].getNNonZeroEvents()-1].fTriggerTimestamp>0)) 
              {
                //LOOP L OVER THE EVENTS IN ACQUISITION J FROM PROF3.
                for(unsigned int l = newStartValue; l < (vec_DetAcqPROF1[j].getNNonZeroEvents()<eventLimit ? vec_DetAcqPROF1[j].getNNonZeroEvents() : eventLimit); l++)
                //for(unsigned int l = newStartValue; l < vec_DetAcqPROF3[j].getNNonZeroEvents(); l++)
                {
                  long long deltaT = vec_CurrentRecordsPROF1[k].fTriggerTimestamp-vec_CurrentRecordsPROF3[l].fTriggerTimestamp;
                  //long long deltaT = vec_CurrentRecordsPROF1[k].fTriggerTimestamp-(vec_CurrentRecordsPROF3[l].fTriggerTimestamp+40);
                  if(deltaT<-1.*tolerancePROF)
                  {
                    break;  
                  }
                  //else if(deltaT>=(-1*tolerancePROF) && deltaT<=0)
                  else if(std::abs(deltaT)<=tolerancePROF)
                  {
                    //LOOP M OVER THE MIDDLE PROFILER'S ACQUISITIONS.
                    for(unsigned int m = 0; m < vec_DetAcqPROF2.size(); m++)
                    {
                      std::vector<AcquisitionXBPF::EventRecordHR> vec_CurrentRecordsPROF2 = vec_DetAcqPROF2[m].getDataHR(); 
                      unsigned int newStartValue2 = 0;
                      bool         beenInRegion2  = false;
                      if(vec_CurrentRecordsPROF2.size()!=0)
                      {
                        //IF THE LATEST MIDTIME IS EARLIER THAN THE CURRENT US TIME OR THE EARLIERST MIDTIME IS LATER THAN CURRENT DS TIME, DONT ACCEPT. 
                        if(vec_CurrentRecordsPROF2[vec_DetAcqPROF2[m].getNNonZeroEvents()-1].fTriggerTimestamp>=vec_CurrentRecordsPROF1[k].fTriggerTimestamp && 
                           vec_CurrentRecordsPROF2[0].fTriggerTimestamp<=vec_CurrentRecordsPROF3[l].fTriggerTimestamp)
                           //vec_CurrentRecordsPROF2[0].fTriggerTimestamp<=vec_CurrentRecordsPROF3[l].fTriggerTimestamp+40)
                        {
                          //LOOP N OVER THE MIDDLE PROFILER'S EVENTS IN ACQUISITION M.
                          for(unsigned int n = newStartValue2; n < vec_DetAcqPROF2[m].getNNonZeroEvents(); n++)
                          {
                            //if(vec_CurrentRecordsPROF2[n].fTriggerTimestamp>vec_CurrentRecordsPROF3[l].fTriggerTimestamp+40)
                            //if(vec_CurrentRecordsPROF2[n].fTriggerTimestamp>vec_CurrentRecordsPROF3[l].fTriggerTimestamp)
                            if(vec_CurrentRecordsPROF2[n].fTriggerTimestamp-vec_CurrentRecordsPROF3[l].fTriggerTimestamp>tolerancePROF)
                            {
                              break;
                            }
                            else if(std::abs(vec_CurrentRecordsPROF2[n].fTriggerTimestamp-vec_CurrentRecordsPROF3[l].fTriggerTimestamp)<=tolerancePROF &&
                                    std::abs(vec_CurrentRecordsPROF2[n].fTriggerTimestamp-vec_CurrentRecordsPROF1[k].fTriggerTimestamp)<=tolerancePROF )
                            //else if(vec_CurrentRecordsPROF2[n].fTriggerTimestamp<=vec_CurrentRecordsPROF3[l].fTriggerTimestamp+40 &&
                            //else if(vec_CurrentRecordsPROF2[n].fTriggerTimestamp<=vec_CurrentRecordsPROF3[l].fTriggerTimestamp &&
                              //  vec_CurrentRecordsPROF2[n].fTriggerTimestamp>=vec_CurrentRecordsPROF1[k].fTriggerTimestamp )
                            {
                              if(beenInRegion2==false)
                              {
                                beenInRegion2 = true;
                                if(n>0)
                                {
                                  newStartValue2 = n - 1;
                                }
                              }
                              std::vector<double> cosTheta; std::vector<double> theta; std::vector<double> momentum;
                              considerMomenta(vec_CurrentRecordsPROF1[k].fFibresList, vec_CurrentRecordsPROF2[n].fFibresList, vec_CurrentRecordsPROF3[l].fFibresList, vec_DetAcqPROF1[i].getCurrent(),
                                              cosTheta, theta, momentum);
                              fPROFCoincdenceRecord.addCoincidence(i, k, m, n, j, l, cosTheta, theta, momentum);
                            }
                          }
                        }
                      }
                    }
                    if(beenInRegion==false)
                    {
                      beenInRegion  = true;
                      if(l>0)
                      {
                        newStartValue = l - 1; 
                      }
                      else
                      {
                        newStartValue = l;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      } 
    }
  }

  return;
}

void BeamLine::printTFCoincidencesDegenerate(std::map<std::string,Detector> &cMapDetectors)
{
  fTFCoincidenceRecord.printDegenerate(cMapDetectors, fXBTFDetNameToIndex);
  return;
}

void BeamLine::printTFCoincidencesUnique(std::map<std::string,Detector> &cMapDetectors)
{
  fTFCoincidenceRecord.printUnique(cMapDetectors, fXBTFDetNameToIndex);
  return;
}

void BeamLine::printPROFCoincidencesDegenerate(std::map<std::string,Detector> &cMapDetectors)
{
  fPROFCoincdenceRecord.printDegenerate(cMapDetectors);
  return;
}

void BeamLine::printPROFCoincidencesUnique(std::map<std::string,Detector> &cMapDetectors)
{
  fPROFCoincdenceRecord.printUnique(cMapDetectors);
  return;
}

void BeamLine::dumpBIData(std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName)
{
  gInterpreter->GenerateDictionary("vector<vector<int> >","vector");

  TFile *f_DumpFile    = new TFile(cFilePathName,"UPDATE");
  TTree *t_DumpBI_XBPF = new TTree("DumpBI_XBPF", "DumpBI_XBPF");

  std::string  out_XBPF_DetName;                             //FULL NAME OF DETECTOR.
  unsigned int out_XBPF_AcquisitionIndex;                    //INDEX OF THE ACQUISITION WITHIN THE DATA DUMP PERIOD.
  double       out_XBPF_AcquisitionTimestamp;                //TIMESTAMP OF THE ACQUISITION, SECONDS.
  unsigned int out_XBPF_NEventsRec;                          //NUMBER OF EVENTS RECORDED IN CURRENT ACQUISITION, COUNTSREC.
  double       out_XBPF_FirstTriggerTimestamp;               //FIRST TRIGGER TIMESTAMP IN THE ACQUISITION, SECONDS.
  double       out_XBPF_LastTriggerTimestamp;                //LAST TRIGGER TIMESTAMP IN THE ACQUISITION, SECONDS.
  double       out_XBPF_Current;                             //CURRENT OF THE BENDING MAGNET, AMPERE.
  std::vector<unsigned int> out_XBPF_EventIndex;             //EVENT INDEX WITHIN THE ACQUISITION.
  std::vector<long long>    out_XBPF_TriggerTimestamp;       //TRIGGER TIMESTAMP OF THE EVENT, NANOSECONDS.
  std::vector<long long>    out_XBPF_EventTimestamp;         //EVENTS TIMESTAMP OF THE EVENT, NANOSECONDS.
  std::vector<unsigned int> out_XBPF_Span;                   //DISTANCE BETWEEN FURTHEST FIBRES HIT IN EVENT, NUMBER OF FIBRES. ZERO WHEN SINGLE FIBRE IS HIT.
  std::vector<unsigned int> out_XBPF_NFibresHit;             //NUMBER OF FIBRES HIT IN AN EVENT.
  std::vector<std::vector<int>> out_XBPF_FibresHit;          //LIST OF FIBRES HIT IN THE EVENT.

  t_DumpBI_XBPF->Branch("DetName",               &out_XBPF_DetName);
  t_DumpBI_XBPF->Branch("AcquisitionIndex",      &out_XBPF_AcquisitionIndex);
  t_DumpBI_XBPF->Branch("AcquisitionTimestamp",  &out_XBPF_AcquisitionTimestamp);
  t_DumpBI_XBPF->Branch("NEventsRec",            &out_XBPF_NEventsRec);
  t_DumpBI_XBPF->Branch("FirstTriggerTimestamp", &out_XBPF_FirstTriggerTimestamp);
  t_DumpBI_XBPF->Branch("LastTriggerTimestamp",  &out_XBPF_LastTriggerTimestamp);
  t_DumpBI_XBPF->Branch("Current",               &out_XBPF_Current);
  t_DumpBI_XBPF->Branch("EventIndex",            &out_XBPF_EventIndex);
  t_DumpBI_XBPF->Branch("TriggerTimestamp",      &out_XBPF_TriggerTimestamp);
  t_DumpBI_XBPF->Branch("EventTimestamp",        &out_XBPF_EventTimestamp);
  t_DumpBI_XBPF->Branch("Span",                  &out_XBPF_Span);
  t_DumpBI_XBPF->Branch("NFibresHit",            &out_XBPF_NFibresHit);
  t_DumpBI_XBPF->Branch("FibresHit",             &out_XBPF_FibresHit);

  for(auto detector : cMapDetectors)
  {
    if(detector.second.getType()=="XBPF")
    {
      std::vector<AcquisitionXBPF> vec_Acq;
      detector.second.getAcquisitions(vec_Acq);
      for(unsigned int i = 0; i < vec_Acq.size(); i++)
      {
        out_XBPF_DetName               = (TString)detector.first;
        out_XBPF_AcquisitionIndex      = i;
        out_XBPF_AcquisitionTimestamp  = vec_Acq[i].getTimestamp();
        out_XBPF_NEventsRec            = vec_Acq[i].getNEventRec();
        out_XBPF_FirstTriggerTimestamp = vec_Acq[i].getFirstTriggerTimestamp();
        out_XBPF_LastTriggerTimestamp  = vec_Acq[i].getLastTriggerTimestamp();
        out_XBPF_Current               = vec_Acq[i].getCurrent();

        std::vector<AcquisitionXBPF::EventRecordHR> vec_currentRecord = vec_Acq[i].getDataHR();
        for(unsigned int j = 0; j < vec_currentRecord.size(); j++)
        {
          out_XBPF_EventIndex      .push_back(j);
          out_XBPF_TriggerTimestamp.push_back(vec_currentRecord[j].fTriggerTimestamp);
          out_XBPF_EventTimestamp  .push_back(vec_currentRecord[j].fEventTimestamp);
          out_XBPF_Span            .push_back(vec_currentRecord[j].fSpan);
          out_XBPF_NFibresHit      .push_back(vec_currentRecord[j].fNFibresHit);
          std::vector<int> temp_Fibres;
          for(unsigned int k = 0; k < vec_currentRecord[j].fFibresList.size(); k++)
          {
            temp_Fibres.push_back(vec_currentRecord[j].fFibresList[k]);
          }
          out_XBPF_FibresHit.push_back(temp_Fibres);
        }
        t_DumpBI_XBPF->Fill();
        out_XBPF_DetName = ""; out_XBPF_AcquisitionIndex = 0; out_XBPF_AcquisitionTimestamp = 0.; out_XBPF_NEventsRec = 0; out_XBPF_FirstTriggerTimestamp = 0.; out_XBPF_LastTriggerTimestamp = 0.;
        out_XBPF_Current = 0.; 
        out_XBPF_EventIndex.clear(); out_XBPF_TriggerTimestamp.clear(); out_XBPF_EventTimestamp.clear(); out_XBPF_Span.clear(); out_XBPF_NFibresHit.clear(); out_XBPF_FibresHit.clear();
      }
    }
  }

  t_DumpBI_XBPF->Write();

  TTree *t_DumpBI_XBTF = new TTree("DumpBI_XBTF", "DumpBI_XBTF");

  std::string  out_XBTF_DetName;                       //FULL NAME OF DETECTOR.
  unsigned int out_XBTF_AcquisitionIndex;              //INDEX OF THE ACQUISITION WITHIN THE DATA DUMP PERIOD.
  double       out_XBTF_AcquisitionTimestamp;          //TIMESTAMP OF THE ACQUISITION, SECONDS.
  unsigned int out_XBTF_NTimestampCount;               //NUMBER OF EVENTS RECORDED IN CURRENT ACQUISITION, COUNTSREC.
  double       out_XBTF_FirstTriggerTimestamp;         //FIRST EVENT TRIGGER TIMESTAMP IN THE ACQUISITION, SECONDS.
  double       out_XBTF_LastTriggerTimestamp;          //LAST EVENT TRIGGER TIMESTAMP IN THE ACQUISITION, SECONDS.
  std::vector<unsigned int>  out_XBTF_EventIndex;      //EVENT INDEX WITHIN THE ACQUISITION.
  std::vector<unsigned long> out_XBTF_TriggerSeconds;  //SECONDS PART OF THE EVENT TRIGGER TIMESTAMP, SECONDS.
  std::vector<double>        out_XBTF_TriggerSubSeconds;    //SUBSECOND PART OF EVENT TRIGGER TIMESTAMP, NANOSECONDS.
  std::vector<double>        out_XBTF_TriggerFullTimestamp; //FULL EVENT TRIGGER TIMESTAMP OF EVENT, SECONDS.

  t_DumpBI_XBTF->Branch("DetName",               &out_XBTF_DetName);
  t_DumpBI_XBTF->Branch("AcquisitionIndex",      &out_XBTF_AcquisitionIndex);
  t_DumpBI_XBTF->Branch("AcquisitionTimestamp",  &out_XBTF_AcquisitionTimestamp);
  t_DumpBI_XBTF->Branch("NTimestampCount",       &out_XBTF_NTimestampCount);
  t_DumpBI_XBTF->Branch("FirstTriggerTimestamp", &out_XBTF_FirstTriggerTimestamp);
  t_DumpBI_XBTF->Branch("LastTriggerTimestamp",  &out_XBTF_LastTriggerTimestamp);
  t_DumpBI_XBTF->Branch("EventIndex",            &out_XBTF_EventIndex);
  t_DumpBI_XBTF->Branch("TriggerSeconds",        &out_XBTF_TriggerSeconds);
  t_DumpBI_XBTF->Branch("TriggerSubSeconds",     &out_XBTF_TriggerSubSeconds);
  t_DumpBI_XBTF->Branch("TriggerFullTimestamp",  &out_XBTF_TriggerFullTimestamp);

  for(auto detector : cMapDetectors)
  {
    if(detector.second.getType()=="XBTF")
    {
      std::vector<AcquisitionXBTF> vec_Acq;
      detector.second.getAcquisitions(vec_Acq);
      for(unsigned int i = 0; i < vec_Acq.size(); i++)
      {
        out_XBTF_DetName               = (TString)detector.first;
        out_XBPF_AcquisitionIndex      = i;
        out_XBTF_AcquisitionTimestamp  = vec_Acq[i].getTimestamp();
        out_XBTF_NTimestampCount       = vec_Acq[i].getTimestampCount();
        out_XBTF_FirstTriggerTimestamp = vec_Acq[i].getFirstTriggerTimestamp();
        out_XBTF_LastTriggerTimestamp  = vec_Acq[i].getLastTriggerTimestamp();

        std::vector<AcquisitionXBTF::EventRecordHR> vec_currentRecord = vec_Acq[i].getDataHR();
        for(unsigned int j = 0; j < vec_currentRecord.size(); j++)
        {
          out_XBTF_EventIndex     .push_back(j);
          out_XBTF_TriggerSeconds .push_back(vec_currentRecord[j].fSeconds);
          out_XBTF_TriggerSubSeconds   .push_back(vec_currentRecord[j].fSubSeconds);
          out_XBTF_TriggerFullTimestamp.push_back(vec_currentRecord[j].fSeconds_FullTime);
        }
        t_DumpBI_XBTF->Fill();
        out_XBTF_DetName = ""; out_XBTF_AcquisitionIndex = 0; out_XBTF_AcquisitionTimestamp = 0.; out_XBTF_NTimestampCount = 0; out_XBTF_FirstTriggerTimestamp = 0.; out_XBTF_LastTriggerTimestamp = 0.;
        out_XBTF_EventIndex.clear(); out_XBTF_TriggerSeconds.clear(); out_XBTF_TriggerSubSeconds.clear(); out_XBTF_TriggerFullTimestamp.clear();
      }
    }
  }

  t_DumpBI_XBTF->Write();

  f_DumpFile->Close();
  delete f_DumpFile;

  return;
}

void BeamLine::dumpTFCoincidencesDegenerate(std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName)
{
  fTFCoincidenceRecord.dumpDegenerate(cMapDetectors, fXBTFDetNameToIndex, cFilePathName);
  return;
}

void BeamLine::dumpTFCoincidencesUnique(std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName)
{
  fTFCoincidenceRecord.dumpUnique(cMapDetectors, fXBTFDetNameToIndex, cFilePathName);
  return;
}

void BeamLine::dumpPROFCoincidencesDegenerate(std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName)
{
  fPROFCoincdenceRecord.dumpDegenerate(cMapDetectors, cFilePathName);
  return;
}

void BeamLine::dumpPROFCoincidencesUnique(std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName)
{
  fPROFCoincdenceRecord.dumpUnique(cMapDetectors, cFilePathName);
  return;
}

void BeamLine::dumpTFPROFCoincidencesDegenerate (std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName)
{
  fCombinedCoincidenceRecord.dumpDegenerate(cMapDetectors, fXBTFDetNameToIndex, cFilePathName, fTFCoincidenceRecord);
  return;
}

void BeamLine::dumpTFPROFCoincidencesUnique(std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName)
{
  fCombinedCoincidenceRecord.dumpUnique(cMapDetectors, fXBTFDetNameToIndex, cFilePathName, fTFCoincidenceRecord);
  return;
}

void BeamLine::findTFPROFCoincidences(std::map<std::string,Detector> &cMapDetectors, bool const &cUseUniqueCoincidencesOnly)
{
  std::vector<TFCoincidenceRecord::TFCoincidence>     vec_TFCo;  
  std::vector<PROFCoincidenceRecord::PROFCoincidence> vec_PROFCo;

  if(cUseUniqueCoincidencesOnly)
  {
    vec_TFCo  = fTFCoincidenceRecord.getTFCoincidencesUnique();
    vec_PROFCo = fPROFCoincdenceRecord.getPROFCoincidencesUnique();
  }
  else
  {
    vec_TFCo  = fTFCoincidenceRecord.getTFCoincidencesDegenerate();
    vec_PROFCo = fPROFCoincdenceRecord.getPROFCoincidencesDegenerate();
  }

  std::vector<AcquisitionXBPF> vec_DetAcqPROF1;
  cMapDetectors["XBPF022697"].getAcquisitions(vec_DetAcqPROF1);
  std::vector<AcquisitionXBPF> vec_DetAcqPROF3;
  cMapDetectors["XBPF022702"].getAcquisitions(vec_DetAcqPROF3);

  unsigned int newStartValue = 0;
  //for(unsigned int i = 0; i < vec_TFCo.size(); i++)
  for(unsigned int i = 0; i < (vec_TFCo.size()<eventLimit ? vec_TFCo.size() : eventLimit); i++)
  {
    if((fXBTFIndexToDetName[vec_TFCo[i].fDetUS]=="XBTF022687A" && fXBTFIndexToDetName[vec_TFCo[i].fDetDS]=="XBTF022716A") ||
       (fXBTFIndexToDetName[vec_TFCo[i].fDetUS]=="XBTF022687A" && fXBTFIndexToDetName[vec_TFCo[i].fDetDS]=="XBTF022716B") ||
       (fXBTFIndexToDetName[vec_TFCo[i].fDetUS]=="XBTF022687B" && fXBTFIndexToDetName[vec_TFCo[i].fDetDS]=="XBTF022716A") ||
       (fXBTFIndexToDetName[vec_TFCo[i].fDetUS]=="XBTF022687B" && fXBTFIndexToDetName[vec_TFCo[i].fDetDS]=="XBTF022716B"))
    {
      std::vector<AcquisitionXBTF> vec_AcqUS;
      cMapDetectors[fXBTFIndexToDetName[vec_TFCo[i].fDetUS]].getAcquisitions(vec_AcqUS);
      std::vector<AcquisitionXBTF> vec_AcqDS;
      cMapDetectors[fXBTFIndexToDetName[vec_TFCo[i].fDetDS]].getAcquisitions(vec_AcqDS);

      bool beenInRegion = false; 
      std::vector<PROFCoincidenceRecord::PROFCoincidence> vec_MatchedPROFCo;
      //for(unsigned int j = newStartValue; j < vec_PROFCo.size(); j++)
      for(unsigned int j = newStartValue; j < (vec_PROFCo.size() < eventLimit ? vec_PROFCo.size() : eventLimit); j++)
      //for(unsigned int j = 0; j < vec_PROFCo.size(); j++)
      {
        if(std::abs((double)vec_DetAcqPROF1[vec_PROFCo[j].fAcqPROF1].getDataHR()[vec_PROFCo[j].fEventPROF1].fTriggerTimestamp
                           -vec_AcqUS[vec_TFCo[i].fAcqUS].getDataHR()[vec_TFCo[i].fEventUS].fSeconds_FullTime*1e9 - TFPROFTimingOffset) < toleranceTFPROF &&
           std::abs((double)vec_DetAcqPROF1[vec_PROFCo[j].fAcqPROF3].getDataHR()[vec_PROFCo[j].fEventPROF3].fTriggerTimestamp
                           -vec_AcqDS[vec_TFCo[i].fAcqDS].getDataHR()[vec_TFCo[i].fEventDS].fSeconds_FullTime*1e9 - TFPROFTimingOffset) < toleranceTFPROF)
        {
          if(beenInRegion==false)
          {
            beenInRegion = true;
            if(j>0)
            {
              newStartValue = j - 1; 
              if(j==vec_PROFCo.size()-1)
              {
                newStartValue = 0;
              }
            }
            else
            {
              newStartValue = 0;
            }
          }
          if(j==vec_PROFCo.size()-1)
          {
            newStartValue = 0;
          }
          vec_MatchedPROFCo.push_back(vec_PROFCo[j]);
          std::array<unsigned int,6> coincidenceIndices = {vec_PROFCo[j].fAcqPROF1, vec_PROFCo[j].fEventPROF1,
                                                           vec_PROFCo[j].fAcqPROF2, vec_PROFCo[j].fEventPROF2,
                                                           vec_PROFCo[j].fAcqPROF3, vec_PROFCo[j].fEventPROF3};
          fCombinedCoincidenceRecord.addMultiplicityInfo(coincidenceIndices);
        }
        else if(beenInRegion==true)
        {
          break;
        }
        if(j==vec_PROFCo.size()-1)
        {
          newStartValue = 0;
          break;
        }
      }
      CombinedCoincidenceRecord::TFPROFCoincidence tfprofCo(vec_TFCo[i], vec_MatchedPROFCo);
      fCombinedCoincidenceRecord.addTFPROFCoincidence(tfprofCo);
    }
  }

  return;
}

void BeamLine::printTFPROFCoincidencesDegenerate(std::map<std::string,Detector> &cMapDetectors)
{
  fCombinedCoincidenceRecord.printDegenerate(cMapDetectors, fXBTFDetNameToIndex, fTFCoincidenceRecord);
  return;
}

void BeamLine::printTFPROFCoincidencesUnique(std::map<std::string,Detector> &cMapDetectors)
{
  fCombinedCoincidenceRecord.printUnique(cMapDetectors, fXBTFDetNameToIndex, fTFCoincidenceRecord);
  return;
}
