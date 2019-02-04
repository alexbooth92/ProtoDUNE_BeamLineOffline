#ifndef CLASS_DETECTOR_H
#define CLASS_DETECTOR_H

#include <string>
#include <cmath>
#include "class_AcquisitionXBPF.h"
#include "class_AcquisitionXBTF.h"

class AcquisitionXBPF;
class AcquisitionXBTF;

class Detector{
  public:
    std::string fDetName;
    std::string fType;
    std::string fDetLabel;
    double fFirstTimestamp;
    double fLastTimestamp;
    double fDetFirstTriggerTimestamp;
    double fDetLastTriggerTimestamp;

    Detector();
    Detector(std::string cDetName, std::string cType, std::string cLabel);
    unsigned int getNAcquisitions();
    double getFirstTimestamp();
    double getLastTimestamp();
    double getTimeElapsed();
    double getFirstTriggerTimestamp();
    double getLastTriggerTimestamp();
    std::string getDetName();
    std::string getDetLabel();
    std::string getType();
    void addAcquisition(AcquisitionXBPF cAcq);
    void addAcquisition(AcquisitionXBTF cAcq);
    void getAcquisitions(std::vector<AcquisitionXBPF> &acq);
    void getAcquisitions(std::vector<AcquisitionXBTF> &acq);
    void printDescription();
    void printRaw();
    void printHR();

  private:
    std::vector<AcquisitionXBPF> fVecAcquisitionXBPF;
    std::vector<AcquisitionXBTF> fVecAcquisitionXBTF;

};

#endif
