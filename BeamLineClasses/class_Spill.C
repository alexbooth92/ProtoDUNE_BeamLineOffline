#include "class_Spill.h"
#include "class_BeamLine.h"

Spill::Spill()
{
  return;
}

Spill::Spill(unsigned long long cTimeS, unsigned long long cTimeNS)
{
  fTimeS  = cTimeS;
  fTimeNS = cTimeNS;
  fFullTimestampNS = (double)fTimeS*1.e9 + (double)fTimeNS;
  fMapTDCTolerance["XBH4.XTDC.022.713"] = 500; fMapTDCTolerance["XBH4.XTDC.022.716"] = 500; fMapTDCTolerance["XBTF.022.687A"] = 500;
  fMapTDCTolerance["XBTF.022.687B"] = 500;   fMapTDCTolerance["XBTF.022.716A"] = 500;   fMapTDCTolerance["XBTF.022.716B"] = 500;

  return;
}

bool Spill::compareTimestamps(BasicDoubleVar const &basicDoubleVar)
{
  bool   isBetterMatch    = false;
  double testTimestamp    = (double)basicDoubleVar.fTimeS;
  double currentTimestamp = (double)fMapBasicDoubleVars[basicDoubleVar.fBranchTitle].fTimeS;

  if(std::abs(testTimestamp-fTimeS)<std::abs(currentTimestamp-fTimeS))
  {
    isBetterMatch = true;
  }

  return isBetterMatch;
}

bool Spill::compareTimestamps(VecBoolVar const &vecBoolVar)
{
  bool   isBetterMatch    = false;
  double testTimestamp    = (double)vecBoolVar.fTimeS;
  double currentTimestamp = fMapXBPFAcquisitions[vecBoolVar.fBranchTitle].fTimestamp;

  if(std::abs(testTimestamp-fTimeS)<std::abs(currentTimestamp-fTimeS))
  {
    isBetterMatch = true;
  }

  return isBetterMatch;
}

bool Spill::compareTimestamps(std::string const &s_Hardware, VecLongVar const &vecLongVar)
{
  bool   isBetterMatch    = false;
  double testTimestamp    = (double)vecLongVar.fTimeS;
  double currentTimestamp = fMapTDCAcquisitions[s_Hardware].fTimestamp;

  if(std::abs(testTimestamp-fTimeS)<std::abs(currentTimestamp-fTimeS))
  {
    isBetterMatch = true;
  }

  return isBetterMatch;
}

bool Spill::isComplete()
{
  if(fMapBasicDoubleVars.size()==fListSpillVars.size())
  {
    return true;
  }
  else
  {
    return false;
  }
}

double Spill::getFractionComplete()
{
  return (double)(fMapBasicDoubleVars.size()+fMapXBPFAcquisitions.size()+fMapTDCAcquisitions.size())/(double)fNVarsPossible;
}

double Spill::getSpillVar(std::string const &s_BranchTitle)
{
  return fMapBasicDoubleVars[s_BranchTitle].fVar;
}

void Spill::addIfBetterMatch(BasicDoubleVar const &basicDoubleVar)
{
  if(fMapBasicDoubleVars.count(basicDoubleVar.fBranchTitle))
  {
    bool isBetterMatch = compareTimestamps(basicDoubleVar);
    if(isBetterMatch)
    {
      fMapBasicDoubleVars[basicDoubleVar.fBranchTitle] = basicDoubleVar;
    }
  }
  else
  {
    fMapBasicDoubleVars[basicDoubleVar.fBranchTitle] = basicDoubleVar;
  }

  return;
}

void Spill::addIfBetterMatch(std::vector<VecBoolVar> const &vec_VecBoolVar)
{
  if(vec_VecBoolVar.size()!=0)
  {
    if(fMapXBPFAcquisitions.count(vec_VecBoolVar[0].fBranchTitle))
    {
      bool isBetterMatch = compareTimestamps(vec_VecBoolVar[0]);
      if(isBetterMatch)
      {
        fMapXBPFAcquisitions[vec_VecBoolVar[0].fBranchTitle].clearAcquisition();
        for(unsigned int i = 0; i < vec_VecBoolVar.size(); i++)
        {
          if(fMapXBPFAcquisitions[vec_VecBoolVar[i].fBranchTitle].fTimestamp==0.)
          {
            fMapXBPFAcquisitions[vec_VecBoolVar[i].fBranchTitle].addAcqDetails(vec_VecBoolVar[0].fTimeS, vec_VecBoolVar.size());
          }
          fMapXBPFAcquisitions[vec_VecBoolVar[i].fBranchTitle].addEventRecord(vec_VecBoolVar[i]);
        }
        fMapXBPFAcquisitions[vec_VecBoolVar[0].fBranchTitle].setCurrent(fMapBasicDoubleVars["XBH4.BEND.022.699:I_MEAS"].fVar);
      }
    }
    else
    {
      if(std::abs((double)fTimeS-(double)vec_VecBoolVar[0].fTimeS) < 15.)
      {
        AcquisitionXBPF acqXBPF;
        acqXBPF.addAcqDetails(vec_VecBoolVar[0].fTimeS, vec_VecBoolVar.size());
        for(unsigned int i = 0; i < vec_VecBoolVar.size(); i++)
        {
          acqXBPF.addEventRecord(vec_VecBoolVar[i]);
        }
        acqXBPF.setCurrent(fMapBasicDoubleVars["XBH4.BEND.022.699:I_MEAS"].fVar);
        fMapXBPFAcquisitions[vec_VecBoolVar[0].fBranchTitle] = acqXBPF;
      }
    }
  }

  return;
}

void Spill::addIfBetterMatch(std::string const &s_Hardware, VecLongVar const &vec_Frac, VecLongVar const &vec_Coarse, VecLongVar const &vec_Seconds, VecLongVar const &vec_TimestampCount)
{
  if(vec_Frac.fVar.size()!=0 && vec_TimestampCount.fVar.size()!=0)
  {
    if(fMapTDCAcquisitions.count(s_Hardware))
    {
      bool isBetterMatch = compareTimestamps(s_Hardware, vec_Frac);
      if(isBetterMatch)
      {
        fMapTDCAcquisitions[s_Hardware].clearAcquisition();
        for(unsigned int i = 0; i < vec_TimestampCount.fVar[0]; i++)
        {
          if(fMapTDCAcquisitions[s_Hardware].fTimestamp==0.)
          {
            fMapTDCAcquisitions[s_Hardware].addAcqDetails(vec_Seconds.fVar[0], vec_TimestampCount.fVar[0]);
          }
          fMapTDCAcquisitions[s_Hardware].addEventRecord({(unsigned int)vec_Frac.fVar[i], vec_Coarse.fVar[i], vec_Seconds.fVar[i]});
        }
        fMapTDCAcquisitions[s_Hardware].sortByNSAccurate();
      }
    }
    else
    {
      if(std::abs((double)fTimeS-(double)vec_Frac.fTimeS) < 15.)
      {
        AcquisitionXBTF acqXBTF;
        acqXBTF.addAcqDetails(vec_Seconds.fVar[0], vec_TimestampCount.fVar[0]);
        for(unsigned int i = 0; i < vec_TimestampCount.fVar[0]; i++)
        {
          acqXBTF.addEventRecord({(unsigned int)vec_Frac.fVar[i], vec_Coarse.fVar[i], vec_Seconds.fVar[i]});
        }
        acqXBTF.sortByNSAccurate();
        fMapTDCAcquisitions[s_Hardware] = acqXBTF;
      }
    }
  }

  return;
}

void Spill::findXBPFEvents()
{
  BeamLine beamline;
  if(fMapXBPFAcquisitions.size()==8)
  {
    if(fMapXBPFAcquisitions["XBH4.XBPF.022.701:EVENTSDATA"].getNEventRec()!=0)
    {
      std::vector<AcquisitionXBPF::EventRecordHR> d697 = fMapXBPFAcquisitions["XBH4.XBPF.022.697:EVENTSDATA"].getDataHR();
      std::vector<AcquisitionXBPF::EventRecordHR> d698 = fMapXBPFAcquisitions["XBH4.XBPF.022.698:EVENTSDATA"].getDataHR();
      std::vector<AcquisitionXBPF::EventRecordHR> d701 = fMapXBPFAcquisitions["XBH4.XBPF.022.701:EVENTSDATA"].getDataHR();
      std::vector<AcquisitionXBPF::EventRecordHR> d702 = fMapXBPFAcquisitions["XBH4.XBPF.022.702:EVENTSDATA"].getDataHR();
      std::vector<AcquisitionXBPF::EventRecordHR> d707 = fMapXBPFAcquisitions["XBH4.XBPF.022.707:EVENTSDATA"].getDataHR();
      std::vector<AcquisitionXBPF::EventRecordHR> d708 = fMapXBPFAcquisitions["XBH4.XBPF.022.708:EVENTSDATA"].getDataHR();
      std::vector<AcquisitionXBPF::EventRecordHR> d716 = fMapXBPFAcquisitions["XBH4.XBPF.022.716:EVENTSDATA"].getDataHR();
      std::vector<AcquisitionXBPF::EventRecordHR> d717 = fMapXBPFAcquisitions["XBH4.XBPF.022.717:EVENTSDATA"].getDataHR();
      for(unsigned int i = 0; i < d697.size(); i++)
      {
        std::vector<double> momenta; std::vector<double> cosTheta; std::vector<double> theta;
        beamline.considerMomenta(d697[i].fFibresList, d701[i].fFibresList, d702[i].fFibresList, fMapBasicDoubleVars["XBH4.BEND.022.699:I_MEAS"].fVar, cosTheta, theta, momenta);
        XBPFCoincidence coincidence = {
          {"XBH4.XBPF.022.697", d697[i].fTriggerTimestamp, d697[i].fEventTimestamp, d697[i].fNFibresHit, d697[i].fSpan, d697[i].fFibresList, beamline.getCoordinates("XBH4.XBPF.022.697", d697[i].fFibresList)},
          {"XBH4.XBPF.022.698", d698[i].fTriggerTimestamp, d698[i].fEventTimestamp, d698[i].fNFibresHit, d698[i].fSpan, d698[i].fFibresList, beamline.getCoordinates("XBH4.XBPF.022.698", d698[i].fFibresList)},
          {"XBH4.XBPF.022.701", d701[i].fTriggerTimestamp, d701[i].fEventTimestamp, d701[i].fNFibresHit, d701[i].fSpan, d701[i].fFibresList, beamline.getCoordinates("XBH4.XBPF.022.698", d701[i].fFibresList)},
          {"XBH4.XBPF.022.702", d702[i].fTriggerTimestamp, d702[i].fEventTimestamp, d702[i].fNFibresHit, d702[i].fSpan, d702[i].fFibresList, beamline.getCoordinates("XBH4.XBPF.022.702", d702[i].fFibresList)},
          {"XBH4.XBPF.022.707", d707[i].fTriggerTimestamp, d707[i].fEventTimestamp, d707[i].fNFibresHit, d707[i].fSpan, d707[i].fFibresList, beamline.getCoordinates("XBH4.XBPF.022.707", d707[i].fFibresList)},
          {"XBH4.XBPF.022.708", d708[i].fTriggerTimestamp, d708[i].fEventTimestamp, d708[i].fNFibresHit, d708[i].fSpan, d708[i].fFibresList, beamline.getCoordinates("XBH4.XBPF.022.708", d708[i].fFibresList)},
          {"XBH4.XBPF.022.716", d716[i].fTriggerTimestamp, d716[i].fEventTimestamp, d716[i].fNFibresHit, d716[i].fSpan, d716[i].fFibresList, beamline.getCoordinates("XBH4.XBPF.022.716", d716[i].fFibresList)},
          {"XBH4.XBPF.022.717", d717[i].fTriggerTimestamp, d717[i].fEventTimestamp, d717[i].fNFibresHit, d717[i].fSpan, d717[i].fFibresList, beamline.getCoordinates("XBH4.XBPF.022.717", d717[i].fFibresList)},
        (unsigned long long)d697[i].fTriggerTimestamp, theta, momenta};
        fVecXBPFCoincidence.push_back(coincidence);
        //coincidence.print();
      }
    }
  }

  return;
}

void Spill::matchXTDCToGeneralTrigger()
{
  std::vector<AcquisitionXBTF::EventRecordHR> gtRecord = fMapTDCAcquisitions["XBH4.GENERALTRIGGER"].getDataHR();
  for(auto it : fMapTDCAcquisitions)
  {
    if(it.first!="XBH4.GENERALTRIGGER")
    //if(it.first!="XBH4.GENERALTRIGGER" && it.first=="XBH4.XTDC.022.713")
    {
      std::vector<AcquisitionXBTF::EventRecordHR> tdcRecord = it.second.getDataHR();
      long long tolerance = fMapTDCTolerance[it.first];
      unsigned int newStartValue = 0;
      for(unsigned int i = 0; i < gtRecord.size(); i++)
      {
        std::vector<unsigned int> vec_MatchEvents;
        bool beenInRegion  = false;
        bool firstInRegion = true;
        for(unsigned int j = newStartValue; j < tdcRecord.size(); j++)
        {
          long long deltaT = (long long)gtRecord[i].fNSAccurate - (long long)tdcRecord[j].fNSAccurate;
          /*
          if(i==0 && j==0)
          {
            std::cout << i << " " << j << " " << gtRecord[i].fNSAccurate << " " <<  tdcRecord[j].fNSAccurate << " " << deltaT << " " << it.first << " BLUD " << (double)tolerance << " " << std::abs((double)deltaT) << std::endl;
          }*/
          //IF THE TDC FIRES WITHIN 500NS OF THE GENERAL TRIGGER, MATCH IT.
          if(std::abs((double)deltaT) < (double)tolerance)
          {
            vec_MatchEvents.push_back(j);
            beenInRegion = true;
            if(firstInRegion && j > 0)
            {
              newStartValue = j - 1;
              firstInRegion = false;
            }
            //std::cout << "accepted, newStartValue: " << newStartValue << std::endl;
          }
          else if(beenInRegion)
          {
            break;
          }
        }
        fMapGenTrigEventsToTDCEvents[{i, it.first}] = vec_MatchEvents;
      }
    }
  }
  
  return;
}

void Spill::matchXBPFToGTNSAccurate()
{
  std::vector<AcquisitionXBTF::EventRecordHR> gtRecord = fMapTDCAcquisitions["XBH4.GENERALTRIGGER"].getDataHR();
  long long tolerance = 500;
  unsigned int newStartValue = 0;
  for(unsigned int i = 0; i < gtRecord.size(); i++)
  {
    bool beenInRegion  = false;
    bool firstInRegion = true;
    for(unsigned int j = newStartValue; j < fVecXBPFCoincidence.size(); j++)
    {
      long long deltaT = (long long)gtRecord[i].fNSAccurate - (long long)fVecXBPFCoincidence[j].fUpstreamTime;
      /*
      if(i<10 && j < 10)
      {
        std::cout << i << " " << j << " " << gtRecord[i].fNSAccurate << " " <<  fVecXBPFCoincidence[j].fUpstreamTime << " " << deltaT << " " << (double)tolerance << " " << std::abs((double)deltaT) << std::endl;
      }*/
      //IF THE UPSTREAM XBPF FIRES WITHIN 500NS OF THE GENERAL TRIGGER, MATCH IT.
      if(std::abs((double)deltaT) < (double)tolerance)
      {
        fMapGenTrigEventToXBPFCoincidence[i].push_back(fVecXBPFCoincidence[j]);
        beenInRegion = true;
        if(firstInRegion && j > 0)
        {
          newStartValue = j - 1;
          firstInRegion = false;
        }
        //std::cout << "accepted, newStartValue: " << newStartValue << std::endl;
      }
      else if(beenInRegion)
      {
        break;
      }
    }
  }

  return;
}

void Spill::matchXBPFToGT()
{
  bool firstIteration      = true;
  bool consistentXBPFCount = true;
  unsigned int nEvents     = 0;
  for(auto it : fMapXBPFAcquisitions)
  {
    if(firstIteration)
    {
      nEvents = it.second.getNEventRec();
      firstIteration = false;
    }
    else
    {
      if(it.second.getNEventRec()!=nEvents)
      {
        consistentXBPFCount = false;
        break;
      }
    }
  }
  fConsistentXBPF = consistentXBPFCount;
  if(consistentXBPFCount)
  {
    findXBPFEvents();
    matchXBPFToGTNSAccurate();
  }

  return;
}

void Spill::getSpillVarList(std::vector<std::string> &vec_SpillVarList)
{
  vec_SpillVarList = fListSpillVars;
  return;
}

void Spill::constructTimeOfFlight()
{
  std::vector<AcquisitionXBTF::EventRecordHR> tdcRecord_687A = fMapTDCAcquisitions["XBTF.022.687A"].getDataHR();
  std::vector<AcquisitionXBTF::EventRecordHR> tdcRecord_687B = fMapTDCAcquisitions["XBTF.022.687B"].getDataHR();
  std::vector<AcquisitionXBTF::EventRecordHR> tdcRecord_716A = fMapTDCAcquisitions["XBTF.022.716A"].getDataHR();
  std::vector<AcquisitionXBTF::EventRecordHR> tdcRecord_716B = fMapTDCAcquisitions["XBTF.022.716B"].getDataHR();

  for(unsigned int i = 0; i < fMapTDCAcquisitions["XBH4.GENERALTRIGGER"].getTimestampCount(); i++)
  {
    std::vector<TOF> vec_TOF;
    if(fMapGenTrigEventsToTDCEvents.count({i,"XBTF.022.687A"}) && fMapGenTrigEventsToTDCEvents.count({i,"XBTF.022.716A"}))
    {
      for(unsigned int j = 0; j < fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.687A"}].size(); j++)
      {
        for(unsigned int k = 0; k < fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.716A"}].size(); k++)
        {
          AcquisitionXBTF::EventRecordHR USEvent = tdcRecord_687A[fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.687A"}][j]];
          AcquisitionXBTF::EventRecordHR DSEvent = tdcRecord_716A[fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.716A"}][k]];
          double deltaSeconds    = (double)DSEvent.fSeconds - (double)USEvent.fSeconds;
          double deltaSubSeconds = DSEvent.fSubSeconds - USEvent.fSubSeconds;
          double tof = deltaSeconds*1e9 + deltaSubSeconds;
          if(tof>0)
          {
            vec_TOF.push_back({"AA", tof}); 
          }
        }
      }
    }
    if(fMapGenTrigEventsToTDCEvents.count({i,"XBTF.022.687A"}) && fMapGenTrigEventsToTDCEvents.count({i,"XBTF.022.716B"}))
    {
      for(unsigned int j = 0; j < fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.687A"}].size(); j++)
      {
        for(unsigned int k = 0; k < fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.716B"}].size(); k++)
        {
          AcquisitionXBTF::EventRecordHR USEvent = tdcRecord_687A[fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.687A"}][j]];
          AcquisitionXBTF::EventRecordHR DSEvent = tdcRecord_716B[fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.716B"}][k]];
          double deltaSeconds    = (double)DSEvent.fSeconds - (double)USEvent.fSeconds;
          double deltaSubSeconds = DSEvent.fSubSeconds - USEvent.fSubSeconds;
          double tof = deltaSeconds*1e9 + deltaSubSeconds;
          if(tof>0)
          {
            vec_TOF.push_back({"AB", tof}); 
          }
        }
      }
    }
    if(fMapGenTrigEventsToTDCEvents.count({i,"XBTF.022.687B"}) && fMapGenTrigEventsToTDCEvents.count({i,"XBTF.022.716A"}))
    {
      for(unsigned int j = 0; j < fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.687B"}].size(); j++)
      {
        for(unsigned int k = 0; k < fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.716A"}].size(); k++)
        {
          AcquisitionXBTF::EventRecordHR USEvent = tdcRecord_687B[fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.687B"}][j]];
          AcquisitionXBTF::EventRecordHR DSEvent = tdcRecord_716A[fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.716A"}][k]];
          double deltaSeconds    = (double)DSEvent.fSeconds - (double)USEvent.fSeconds;
          double deltaSubSeconds = DSEvent.fSubSeconds - USEvent.fSubSeconds;
          double tof = deltaSeconds*1e9 + deltaSubSeconds;
          if(tof>0)
          {
            vec_TOF.push_back({"BA", tof}); 
          }
        }
      }
    }
    if(fMapGenTrigEventsToTDCEvents.count({i,"XBTF.022.687B"}) && fMapGenTrigEventsToTDCEvents.count({i,"XBTF.022.716B"}))
    {
      for(unsigned int j = 0; j < fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.687B"}].size(); j++)
      {
        for(unsigned int k = 0; k < fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.716B"}].size(); k++)
        {
          AcquisitionXBTF::EventRecordHR USEvent = tdcRecord_687B[fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.687B"}][j]];
          AcquisitionXBTF::EventRecordHR DSEvent = tdcRecord_716B[fMapGenTrigEventsToTDCEvents[{i,"XBTF.022.716B"}][k]];
          double deltaSeconds    = (double)DSEvent.fSeconds - (double)USEvent.fSeconds;
          double deltaSubSeconds = DSEvent.fSubSeconds - USEvent.fSubSeconds;
          double tof = deltaSeconds*1e9 + deltaSubSeconds;
          if(tof>0)
          {
            vec_TOF.push_back({"BB", tof}); 
          }
        }
      }
    }
    fMapGenTrigEventToTOFs[i] = vec_TOF;
  }

  return;
}

void Spill::printTDCMatching()
{
  std::vector<AcquisitionXBTF::EventRecordHR> gtRecord = fMapTDCAcquisitions["XBH4.GENERALTRIGGER"].getDataHR();
  for(auto it : fMapGenTrigEventsToTDCEvents)
  {
    std::vector<AcquisitionXBTF::EventRecordHR> tdcRecord = fMapTDCAcquisitions[it.first.second].getDataHR();
    AcquisitionXBTF::EventRecordHR gtEvent = gtRecord[it.first.first];
    if(it.second.size()!=0 && it.first.second.find("XBH4")==std::string::npos)
    {
      std::cout << "GENERAL TRIGGER: " << it.first.first << " MATCHED WITH: " << it.first.second << ", THERE ARE " << it.second.size() << " POSSIBLE EVENTS. " << std::endl; 
      for(unsigned int i = 0; i < it.second.size(); i++)
      {
        std::cout << "EVENT: " << it.second[i] << " AT TIME: " << tdcRecord[it.second[i]].fNSAccurate << " ";
      }
      std::cout << std::endl;
    }
  }

  return;
}

void Spill::printSpill()
{
  std::cout << "**************************************************************************************************************************" << std::endl;
  std::cout << "*  SPILL INFO; TIMESTAMP S: " << fTimeS << ", NS: " << fTimeNS 
            << ", NVARS ASSIGNED: "           << (fMapBasicDoubleVars.size()+fMapXBPFAcquisitions.size()+fMapTDCAcquisitions.size()) << "/" << fNVarsPossible << "." << std::endl;
  std::cout << "**************************************************************************************************************************" << std::endl;
  
  unsigned int count = 0;
  std::cout << "--------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "- SPILL VARS: " << std::endl; 
  for(auto it : fMapBasicDoubleVars)
  {
    //std::cout << it.first << ": " << it.second.fVar << ". ";  
    std::cout << it.first << ": " << it.second.fVar << ", AT TIME " << it.second.fTimeS << "s " << it.second.fTimeNS << "ns. " << std::endl; 
  }
  std::cout << "\n--------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "- XBPFs:" << std::endl;
  for(auto it : fMapXBPFAcquisitions)
  {
    std::cout << it.first << std::endl;
    it.second.printHR();
  }
  std::cout << "\n--------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "- XTDCs:" << std::endl;
  for(auto it : fMapTDCAcquisitions)
  {
    std::cout << it.first << std::endl;
    it.second.printHR();
  }
  std::cout << "**************************************************************************************************************************" << std::endl;
  std::cout << std::endl;

  return;
}

void Spill::printSpillConcise()
{
  std::cout.precision(25);
  std::cout << "**************************************************************************************************************************" << std::endl;
  std::cout << "*  SPILL INFO; TIMESTAMP S: " << fTimeS << ", NS: " << fTimeNS 
            << ", NVARS ASSIGNED: "           << (fMapBasicDoubleVars.size()+fMapXBPFAcquisitions.size()+fMapTDCAcquisitions.size()) << "/" << fNVarsPossible << "." << std::endl;
  std::cout << "**************************************************************************************************************************" << std::endl;
  
  unsigned int count = 0;
  std::cout << "--------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "- SPILL VARS: " << std::endl; 
  for(auto it : fMapBasicDoubleVars)
  {
    //std::cout << it.first << ": " << it.second.fVar << ". ";  
    std::cout << it.first << ": " << it.second.fVar << ", AT TIME " << it.second.fTimeS << "s " << it.second.fTimeNS << "ns. " << std::endl; 
  }
  std::cout << "\n--------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "- XBPFs:" << std::endl;
  for(auto it : fMapXBPFAcquisitions)
  {
    std::cout << it.first << " AT TIME: " << it.second.fTimestamp << std::endl;
  }
  std::cout << "\n--------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "- XTDCs:" << std::endl;
  for(auto it : fMapTDCAcquisitions)
  {
    std::cout << it.first << " AT TIME: " << it.second.fTimestamp << std::endl;
  }
  std::cout << "**************************************************************************************************************************" << std::endl;
  std::cout << std::endl;

  return;
}
