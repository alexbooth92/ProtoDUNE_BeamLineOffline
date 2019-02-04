#ifndef CLASS_SPILL_H
#define CLASS_SPILL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include "class_AcquisitionXBPF.h"
#include "class_AcquisitionXBTF.h"

class AcquisitionXBPF;
class AcquisitionXBTF;
class BeamLine;

struct BasicIntVar{
    std::string        fBranchTitle;
    unsigned long long fTimeS;	
    unsigned long long fTimeNS;	
    int fVar;
};

struct BasicDoubleVar{
    std::string        fBranchTitle;
    unsigned long long fTimeS;	
    unsigned long long fTimeNS;	
    double fVar;
};

struct VecBoolVar{
    std::string fBranchTitle;
    unsigned long long fTimeS;
    unsigned long long fTimeNS;
    std::vector<bool>  fVar;
};

struct VecLongVar{
    std::string fBranchTitle;
    unsigned long long fTimeS;
    unsigned long long fTimeNS;
    std::vector<unsigned long long>  fVar;
};

struct XBPFEvent{
    std::string fDetector;
    long long fTriggerTimestamp;
    long long fEventTimestamp;
    unsigned int fFibresHit;
    unsigned int fSpan;
    std::vector<unsigned int> fFibreList;
    std::vector<double> fPhysicalCoordinates;
    void print(){std::cout << fDetector << " TRIGGERED AT TIME: " << fTriggerTimestamp << ", N FIBRES HIT: " << fFibresHit << ", WITH SPAN: " << fSpan << std::endl;}
};

struct XBPFCoincidence{
  XBPFEvent f697;
  XBPFEvent f698;
  XBPFEvent f701;
  XBPFEvent f702;
  XBPFEvent f707;
  XBPFEvent f708;
  XBPFEvent f716;
  XBPFEvent f717;
  unsigned long long fUpstreamTime;
  std::vector<double> fDeflectionAngle;
  std::vector<double> fMomentum;
  void print(){
    std::cout << "\nCOINCIDENCE IN XBPFs, UPSTREAM TIME: " << fUpstreamTime << std::endl;
    f697.print(); f698.print(); f701.print(); f702.print();
    f707.print(); f708.print(); f716.print(); f717.print();
    std::cout << "THE PARTICLE HAD THE FOLLOWING POSSIBLE DEFLECTION ANGLES AND MOMENTA: " << std::endl;
    for(unsigned int i = 0; i < fDeflectionAngle.size(); i++)
    {
      std::cout << "ANGLE: " << fDeflectionAngle[i] << " DEGREES, " << fMomentum[i] << " GeV" << std::endl;
    }
    std::cout << std::endl;
  };
};

struct TOF{
  std::string fChannel;
  double fTOF;
};

class Spill{

  public:
    unsigned long long fTimeS;	
    unsigned long long fTimeNS;	
    double fFullTimestampNS;

    Spill();
    Spill(unsigned long long cTimeS, unsigned long long cTimeNS);
    bool   isComplete();
    bool   areXBPFConsistent(){return fConsistentXBPF;};
    double getFractionComplete();
    double getSpillVar(std::string const &s_BranchTitle);
    void   addIfBetterMatch(BasicDoubleVar const &basicDoubleVar);
    void   addIfBetterMatch(std::vector<VecBoolVar> const &vec_VecBoolVar);
    void   addIfBetterMatch(std::string const &s_Hardware, VecLongVar const &vec_Frac, VecLongVar const &vec_Coarse, VecLongVar const &vec_Seconds, VecLongVar const &vec_TimestampCount);
    void   constructTimeOfFlight();
    void   findXBPFEvents();
    void   findXTOFEvents();
    void   getSpillVarList(std::vector<std::string> &vec_SpillVarList);
    void   matchXBPFToGT();
    void   matchXTDCToGeneralTrigger();
    void   matchXBPFToGTNSAccurate();
    void   printSpill();
    void   printSpillConcise();
    void   printTDCMatching();
    std::map<std::string,BasicDoubleVar>  getMapBasicDoubleVars(){return fMapBasicDoubleVars;};
    std::map<std::string,AcquisitionXBTF> getMapTDCAcquisitions(){return fMapTDCAcquisitions;};
    std::map<std::pair<unsigned int,std::string>,std::vector<unsigned int>> getMapGenTrigToTDCEvents(){return fMapGenTrigEventsToTDCEvents;};
    std::map<unsigned int,std::vector<TOF>> getMapGenTrigToEventTOFs(){return fMapGenTrigEventToTOFs;};
    std::map<unsigned int,std::vector<XBPFCoincidence>>  getMapGenTrigToXBPFCoincidence(){return fMapGenTrigEventToXBPFCoincidence;}; 
    std::vector<XBPFCoincidence> getVecXBPFCoincidence(){return fVecXBPFCoincidence;};

  private:                            	
    unsigned int fNVarsPossible = 35; 
    bool fConsistentXBPF;
    std::map<std::string,long long> fMapTDCTolerance;

    std::vector<std::string> fListSpillVars = {
      "XBH4.BEND.022.692:I_MEAS",	
      "XBH4.BEND.022.699:I_MEAS",	
      "XBH4.EXPT.NP04.001:COUNTS",
      "XBH4.EXPT.NP04.002:COUNTS",
      "XBH4.EXPT.NP04.003:COUNTS",
      "XBH4.EXPT.NP04.004:COUNTS",
      "XBH4.EXPT.NP04.009:COUNTS",
      "XBH4.EXPT.NP04.010:COUNTS",
      "XBH4.EXPT.NP04.011:COUNTS",
      "XBH4.EXPT.NP04.012:COUNTS",
      "XBH4.XCET.022.713:COUNTS_TRIG",
      "XBH4.XCET.022.713:PRESSURE",
      "XBH4.XCET.022.713:SIMPLE_COUNTS",
      "XBH4.XCET.022.716:COUNTS_TRIG",
      "XBH4.XCET.022.716:PRESSURE",
      "XBH4.XCET.022.716:SIMPLE_COUNTS",
      "XBH4.XCSH.022.694:POS_JAW1_MEAS",
      "XBH4.XCSH.022.694:POS_JAW2_MEAS",
      "XBH4.XSCI.022.680:COUNTS", 
      "SPS.T2:INTENSITY"
    };

    //SPILL VARS
    std::map<std::string,BasicDoubleVar> fMapBasicDoubleVars;
    std::map<std::string,BasicIntVar>    fMapBasicIntVars;

    //XBPFS
    std::map<std::string,AcquisitionXBPF> fMapXBPFAcquisitions;

    //TDCs
    std::map<std::string,AcquisitionXBTF> fMapTDCAcquisitions;

    //XBPF EVENTS.
    std::vector<XBPFCoincidence> fVecXBPFCoincidence;
    std::map<unsigned int,std::vector<XBPFCoincidence>> fMapGenTrigEventToXBPFCoincidence;

    //GENERAL TRIGGER & OTHER TDC ASSOCIATION.
    std::map<std::pair<unsigned int,std::string>,std::vector<unsigned int>> fMapGenTrigEventsToTDCEvents;
    std::map<unsigned int,std::vector<TOF>> fMapGenTrigEventToTOFs;

    bool compareTimestamps(BasicDoubleVar const &basicDoubleVar);                      
    bool compareTimestamps(VecBoolVar     const &vecBoolVar);
    bool compareTimestamps(std::string    const &s_Hardware, VecLongVar const &vecLongVar);
};

#endif
