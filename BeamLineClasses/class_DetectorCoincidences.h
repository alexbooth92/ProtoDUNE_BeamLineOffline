#ifndef CLASS_DETECTORCOINCIDENCES_H
#define CLASS_DETECTORCOINCIDENCES_H

#include <vector>
#include <string>
#include <map>
#include <array>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "class_Detector.h"
#include "class_AcquisitionXBTF.h"
#include "class_AcquisitionXBPF.h"

class Detector;
class AcquisitionXBTF;
class AcquisitionXBPF;

class TFCoincidenceRecord{

  public:
    struct TFCoincidence{
        unsigned int fDetUS;
        unsigned int fAcqUS;
        unsigned int fEventUS;
        unsigned int fDetDS;
        unsigned int fAcqDS;
        unsigned int fEventDS;
        double       fTF;
        std::pair<unsigned int,unsigned int> fMultiplicity;
        bool fUniqueness;
        void setMultiplicity(std::pair<unsigned int,unsigned int> &cMultiplicity){fMultiplicity = cMultiplicity; return;};
        void setUniqueness(bool const &cIsUnique){fUniqueness = cIsUnique; return;};
    };

    TFCoincidenceRecord(){return;};
    std::vector<TFCoincidence> getTFCoincidencesDegenerate();
    std::vector<TFCoincidence> getTFCoincidencesUnique();
    std::pair<unsigned int,unsigned int> getMultiplicity(TFCoincidence const &cTFCoincidence);
    unsigned int getMultiplicity(std::string const &cUSOrDS, unsigned int const &cDet, unsigned int const &cAcq, unsigned int const &cEvent);
    void addCoincidence(unsigned int cDetUS, unsigned int cAcqUS, unsigned int cEventUS, unsigned int cDetDS, unsigned int cAcqDS, unsigned int cEventDS, double cTF);
    void dumpDegenerate (std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TString const &cFilePathName);
    void dumpUnique     (std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TString const &cFilePathName);
    void printDegenerate(std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex);
    void printUnique    (std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex);

  private:
    std::vector<TFCoincidence> fVecTFCoincidesDegenerate;
    std::vector<TFCoincidence> fVecTFCoincidesUnique;
    std::map<std::array<unsigned int,3>,std::vector<TFCoincidence>> fMapUS;
    std::map<std::array<unsigned int,3>,std::vector<TFCoincidence>> fMapDS;
};

class PROFCoincidenceRecord{

  public:
    class PROFCoincidence{

      public:
        unsigned int fAcqPROF1;
        unsigned int fEventPROF1;
        unsigned int fAcqPROF2;
        unsigned int fEventPROF2;
        unsigned int fAcqPROF3;
        unsigned int fEventPROF3;
        std::vector<double> fCosTheta;
        std::vector<double> fTheta;
        std::vector<double> fMomentum;
        std::array<long unsigned int,3> fMultiplicity;
        bool fUniqueness;

        PROFCoincidence(unsigned int cAcqPROF1, unsigned int cEventPROF1, unsigned int cAcqPROF2, unsigned int cEventPROF2, unsigned int cAcqPROF3, unsigned int cEventPROF3)
        {
          fAcqPROF1   = cAcqPROF1;
          fEventPROF1 = cEventPROF1;
          fAcqPROF2   = cAcqPROF2;
          fEventPROF2 = cEventPROF2;
          fAcqPROF3   = cAcqPROF3;
          fEventPROF3 = cEventPROF3;

          return;
        };

        void addDeflectionInfo(std::vector<double> cCosTheta, std::vector<double> cTheta, std::vector<double> cMomentum)
        {
          fCosTheta = cCosTheta;
          fTheta    = cTheta;
          fMomentum = cMomentum;

          return;
        };

        std::vector<double> getCosTheta()
        {
          return fCosTheta;
        };
        std::vector<double> getTheta()
        {
          return fTheta;
        };
        std::vector<double> getMomentum()
        {
          return fMomentum;
        };
        unsigned int getNFibreCombinations()
        {
          return fMomentum.size();
        };
        void setMultiplicity(std::array<long unsigned int,3> &cMultiplicity){fMultiplicity = cMultiplicity; return;};
        void setUniqueness(bool const&cIsUnique){fUniqueness = cIsUnique; return;}
    };
    PROFCoincidenceRecord(){};
    std::vector<PROFCoincidence> getPROFCoincidencesDegenerate();
    std::vector<PROFCoincidence> getPROFCoincidencesUnique();
    std::array<long unsigned int,3> getMultiplicity(PROFCoincidence const &cPROFCoincidence);
    unsigned int getMultiplicity(unsigned int const &cDet, unsigned int const &cAcq, unsigned int const &cEvent);
    void addCoincidence(unsigned int cAcqPROF1, unsigned int cEventPROF1, unsigned int cAcqPROF2, unsigned int cEventPROF2, unsigned int cAcqPROF3, unsigned int cEventPROF3,
                        std::vector<double> cCosTheta, std::vector<double> cTheta, std::vector<double> cMomentum);
    void dumpDegenerate (std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName);
    void dumpUnique     (std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName);
    void printDegenerate(std::map<std::string,Detector> &cMapDetectors);
    void printUnique    (std::map<std::string,Detector> &cMapDetectors);

  private:
    std::vector<PROFCoincidence> fVecPROFCoincidesDegenerate;
    std::vector<PROFCoincidence> fVecPROFCoincidesUnique;
    std::map<std::array<unsigned int,2>,std::vector<PROFCoincidence>> fMapPROF1;
    std::map<std::array<unsigned int,2>,std::vector<PROFCoincidence>> fMapPROF2;
    std::map<std::array<unsigned int,2>,std::vector<PROFCoincidence>> fMapPROF3;
};

class CombinedCoincidenceRecord{
  public:
    class TFPROFCoincidence{
      public:
        //UNTIL WE KNOW THAT THE PATH LENGTH IS DEFINITELY FIXED, LEAVE THE FLEXIBILITY IN THERE.
        TFCoincidenceRecord::TFCoincidence fTFCo;
        std::vector<PROFCoincidenceRecord::PROFCoincidence> fPROFCos;
        std::vector<std::vector<double>> fMass;

        TFPROFCoincidence(TFCoincidenceRecord::TFCoincidence cTFCo, std::vector<PROFCoincidenceRecord::PROFCoincidence> cPROFCos){fTFCo = cTFCo; fPROFCos = cPROFCos; addMasses(); return;};
        unsigned int getMultiplicity()
        {
          return fPROFCos.size();
        }
        std::vector<PROFCoincidenceRecord::PROFCoincidence> getPROFCoincidences()
        {
          return fPROFCos;
        }
        TFCoincidenceRecord::TFCoincidence getTFCoincidence()
        {
          return fTFCo;
        }
        void addPROFCoincidence(PROFCoincidenceRecord::PROFCoincidence cPROFCoincidence)
        {
          fPROFCos.push_back(cPROFCoincidence);
          std::vector<double> tempMass;
          for(unsigned int i = 0; i < cPROFCoincidence.getNFibreCombinations(); i++)
          {
            tempMass.push_back(calculateMass(cPROFCoincidence.getMomentum()[i], 29.00, fTFCo.fTF));  
          }
          fMass.push_back(tempMass);
        }
        void addMasses()
        {
          for(unsigned int i = 0; i < fPROFCos.size(); i++)
          {
            std::vector<double> tempMass;
            for(unsigned int j = 0; j < fPROFCos[i].getNFibreCombinations(); j++)
            {
              tempMass.push_back(calculateMass(fPROFCos[i].getMomentum()[j], 29.00, fTFCo.fTF));  
            }
            fMass.push_back(tempMass);
          }
        }

      private:
        double calculateMass(double const &cMomentum, double const &cPathLength, double const &cTimeOfFlight)
        {
          double a    = 299792458*cTimeOfFlight/cPathLength;
          double mass = (a*a - 1 >= 0 ? (cMomentum/299792458)*std::sqrt(a*a-1) : -1.);

          return mass;
        }
    };
    CombinedCoincidenceRecord(){return;};
    std::vector<TFPROFCoincidence> getTFPROFCoincidencesDegenerate();
    std::vector<TFPROFCoincidence> getTFPROFCoincidencesUnique();
    unsigned int getMultiplicity(PROFCoincidenceRecord::PROFCoincidence const &cPROFCoincidence); //RETURNS THE MULTIPLICITY OF A PARTICULAR PROF COMBINATION ACROSS ALL OF THE TF COINCIDENCES.
    void addTFPROFCoincidence(TFPROFCoincidence cTFPROFCoincidence);
    void addMultiplicityInfo(std::array<unsigned int,6> &cPROFCoincidenceIncides);
    void dumpDegenerate (std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TString const &cFilePathName, TFCoincidenceRecord &cTFCoincidenceRecord);
    void dumpUnique     (std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TString const &cFilePathName, TFCoincidenceRecord &cTFCoincidenceRecord);
    void printDegenerate(std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TFCoincidenceRecord &cTFCoincidenceRecord);
    void printUnique    (std::map<std::string,Detector> &cMapDetectors, std::map<std::string,unsigned int> &cXBTFNameToIndex, TFCoincidenceRecord &cTFCoincidenceRecord);

  private:
    std::vector<TFPROFCoincidence>  fTFPROFCoincidencesDegenerate;
    std::vector<TFPROFCoincidence>  fTFPROFCoincidencesUnique;
    std::map<std::array<unsigned int,6>,unsigned int> fMapPROFCount;

};

#endif
