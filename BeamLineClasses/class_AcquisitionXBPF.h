#ifndef CLASS_AQCUISITIONXBPF_H
#define CLASS_AQCUISITIONXBPF_H

#include <iostream>
#include <array>
#include <bitset>
#include <vector>
#include <map>
#include "class_Spill.h"

struct VecBoolVar;

class AcquisitionXBPF{

  public:
    bool   fIsEmpty;
    double fTimestamp;
    unsigned int fNEventRec;
    unsigned int fNNonZeroEvents;
    double fFirstTriggerTimestamp;
    double fLastTriggerTimestamp;
    double fCurrent;
    struct EventRecordRaw{
           unsigned long fTTS_LSB;
           unsigned long fTTS_MSB;
           unsigned long fETS_LSB;
           unsigned long fETS_MSB;
           std::array<unsigned int, 6> fArrFibres;
    };
    struct EventRecordHR{
           long long fTriggerTimestamp;
           long long fEventTimestamp;
           std::array<bool, 192> fArrFibres;
           std::vector<unsigned int> fFibresList;
           unsigned int fNFibresHit;
           unsigned int fSpan;
    };
    
    AcquisitionXBPF();
    bool         isEmpty();
    double       getFirstTriggerTimestamp();
    double       getLastTriggerTimestamp();
    double       getTimestamp();
    double       getCurrent();
    unsigned int getNEventRec();
    unsigned int getNNonZeroEvents(){return fNNonZeroEvents;};
    std::vector<EventRecordRaw> getDataRaw();
    std::vector<EventRecordHR>  getDataHR();
    void addAcqDetails(double cTimestamp, unsigned int cNEventRec);
    void addEventRecord(EventRecordRaw cEventRecord);
    void addEventRecord(VecBoolVar      cVecBoolVar);
    void clearAcquisition();
    void printRaw();
    void printHR();
    void setEmpty(bool cIsEmpty);
    void setCurrent(double cCurrent);
    void setNNonZeroEvents(unsigned int cNNonZeroEvents){fNNonZeroEvents = cNNonZeroEvents; return;}

  private:
    EventRecordHR fEventRecordHR;
    std::vector<EventRecordRaw> fVecDataRec;
    std::vector<EventRecordHR> fVecDataRecHR;

    EventRecordHR decodeFibres(EventRecordRaw const &cVecDataRec);
    EventRecordHR decodeFibres(VecBoolVar     const &cVecDataRec);

};

#endif
