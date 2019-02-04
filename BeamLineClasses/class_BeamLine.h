#ifndef CLASS_BEAMLINE_H
#define CLASS_BEAMLINE_H

#include <map>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "class_Detector.h"
#include "class_DetectorCoincidences.h"

class TFCoincidenceRecord;
class PROFCoincidenceRecord;

class BeamLine{

  public:
    BeamLine();
    std::vector<double> getCoordinates(std::string const &cDet, std::vector<unsigned int> const &cFibre);
    std::vector<TFCoincidenceRecord::TFCoincidence> getTFCoincidencesDegenerate(); 
    std::vector<TFCoincidenceRecord::TFCoincidence> getTFCoincidencesUnique(); 
    std::vector<PROFCoincidenceRecord::PROFCoincidence> getPROFCoincidencesDegenerate(); 
    std::vector<PROFCoincidenceRecord::PROFCoincidence> getPROFCoincidencesUnique(); 
    std::vector<CombinedCoincidenceRecord::TFPROFCoincidence> getTFPROFCoincidencesDegenerate(); 
    std::vector<CombinedCoincidenceRecord::TFPROFCoincidence> getTFPROFCoincidencesUnique(); 
    void dumpBIData                       (std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName);
    void dumpTFCoincidencesDegenerate     (std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName);
    void dumpTFCoincidencesUnique         (std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName);
    void dumpPROFCoincidencesDegenerate   (std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName);
    void dumpPROFCoincidencesUnique       (std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName);
    void dumpTFPROFCoincidencesDegenerate (std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName);
    void dumpTFPROFCoincidencesUnique     (std::map<std::string,Detector> &cMapDetectors, TString const &cFilePathName);
    void findTFCoincidences               (std::map<std::string,Detector> &cMapDetectors); 
    void findPROFCoincidences             (std::map<std::string,Detector> &cMapDetectors); 
    void findTFPROFCoincidences           (std::map<std::string,Detector> &cMapDetectors, bool const &cUseUniqueCoincidencesOnly); 
    void printTFCoincidencesDegenerate    (std::map<std::string,Detector> &cMapDetectors);
    void printTFCoincidencesUnique        (std::map<std::string,Detector> &cMapDetectors);
    void printPROFCoincidencesDegenerate  (std::map<std::string,Detector> &cMapDetectors);
    void printPROFCoincidencesUnique      (std::map<std::string,Detector> &cMapDetectors);
    void printTFPROFCoincidencesDegenerate(std::map<std::string,Detector> &cMapDetectors); 
    void printTFPROFCoincidencesUnique    (std::map<std::string,Detector> &cMapDetectors); 
    void considerMomenta(std::vector<unsigned int> const &cPROF1Fibs, std::vector<unsigned int> const &cPROF2Fibs, std::vector<unsigned int> const &cPROF3Fibs, double const &cCurrent,
                         std::vector<double> &cCosTheta, std::vector<double> &cTheta, std::vector<double> &cMomentum);

  private:
    TFCoincidenceRecord       fTFCoincidenceRecord;
    PROFCoincidenceRecord     fPROFCoincdenceRecord;
    CombinedCoincidenceRecord fCombinedCoincidenceRecord;
    std::map<std::string,unsigned int> fXBTFDetNameToIndex;
    std::map<unsigned int,std::string> fXBTFIndexToDetName;

    double getCoordinate(unsigned int const &cDet,   unsigned int const &cFibre);
    double getCosTheta  (unsigned int const &cPROF1, unsigned int const &cPROF2, unsigned int const &cPROF4);
    double getMomentum  (double const &cCosTheta, double const &cCurrent);

    unsigned int eventLimit = 500;
    //NANOSECONDS.
    double    toleranceTF     = 500.;
    double    toleranceTFPROF = 500.;
    long long tolerancePROF   = 100.;

    double    TFPROFTimingOffset = 1000.;
    
    //GEOMETRY.
    //RADIANS, METERS FOR OFFSETS, MM FOR FIBRE COORDS (CORRECTED IN THE GET FUNCTION).
    double mag_P1 = 5.82044830e-3;
    double mag_P2 = 0.;
    double mag_P3 = -4.68880000e-6;
    double mag_P4 = 3.24573967e2;
    double Pi     = 3.14159265;
    double Theta0 = 120.030e-3;
    double Offset_XBPF022697 = 1.961;
    //double Offset_XBPF022697 = 2.004;
    std::vector<double> FibreToCoord_XBPF022697 = {
      -95.5, -94.5, -93.5, -92.5, -91.5, -90.5, -89.5, -88.5, -87.5, -86.5, -85.5, -84.5, 
      -83.5, -82.5, -81.5, -80.5, -79.5, -78.5, -77.5, -76.5, -75.5, -74.5, -73.5, -72.5,
      -71.5, -70.5, -69.5, -68.5, -67.5, -66.5, -65.5, -64.5, -63.5, -62.5, -61.5, -60.5,
      -59.5, -58.5, -57.5, -56.5, -55.5, -54.5, -53.5, -52.5, -51.5, -50.5, -49.5, -48.5,
      -47.5, -46.5, -45.5, -44.5, -43.5, -42.5, -41.5, -40.5, -39.5, -38.5, -37.5, -36.5,
      -35.5, -34.5, -33.5, -32.5, -31.5, -30.5, -29.5, -28.5, -27.5, -26.5, -25.5, -24.5,
      -23.5, -22.5, -21.5, -20.5, -19.5, -18.5, -17.5, -16.5, -15.5, -14.5, -13.5, -12.5,
      -11.5, -10.5,  -9.5,  -8.5,  -7.5,  -6.5,  -5.5,  -4.5,  -3.5,  -2.5,  -1.5,  -0.5,
        0.5,   1.5,   2.5,   3.5,   4.5,   5.5,   6.5,   7.5,   8.5,   9.5,  10.5,  11.5,
       12.5,  13.5,  14.5,  15.5,  16.5,  17.5,  18.5,  19.5,  20.5,  21.5,  22.5,  23.5,
       24.5,  25.5,  26.5,  27.5,  28.5,  29.5,  30.5,  31.5,  32.5,  33.5,  34.5,  35.5, 
       36.5,  37.5,  38.5,  39.5,  40.5,  41.5,  42.5,  43.5,  44.5,  45.5,  46.5,  47.5,
       48.5,  49.5,  50.5,  51.5,  52.5,  53.5,  54.5,  55.5,  56.5,  57.5,  58.5,  59.5,
       60.5,  61.5,  62.5,  63.5,  64.5,  65.5,  66.5,  67.5,  68.5,  69.5,  70.5,  71.5,
       72.5,  73.5,  74.5,  75.5,  76.5,  77.5,  78.5,  79.5,  80.5,  81.5,  82.5,  83.5,
       84.5,  85.5,  86.5,  87.5,  88.5,  89.5,  90.5,  91.5,  92.5,  93.5,  94.5, 95.5};

    //double Offset_XBPF022701 = 1.726;
    double Offset_XBPF022701 = 1.71358;
    //double Offset_XBPF022701 = 1.718;
    std::vector<double> FibreToCoord_XBPF022701 = {
      -95.5, -94.5, -93.5, -92.5, -91.5, -90.5, -89.5, -88.5, -87.5, -86.5, -85.5, -84.5, 
      -83.5, -82.5, -81.5, -80.5, -79.5, -78.5, -77.5, -76.5, -75.5, -74.5, -73.5, -72.5,
      -71.5, -70.5, -69.5, -68.5, -67.5, -66.5, -65.5, -64.5, -63.5, -62.5, -61.5, -60.5,
      -59.5, -58.5, -57.5, -56.5, -55.5, -54.5, -53.5, -52.5, -51.5, -50.5, -49.5, -48.5,
      -47.5, -46.5, -45.5, -44.5, -43.5, -42.5, -41.5, -40.5, -39.5, -38.5, -37.5, -36.5,
      -35.5, -34.5, -33.5, -32.5, -31.5, -30.5, -29.5, -28.5, -27.5, -26.5, -25.5, -24.5,
      -23.5, -22.5, -21.5, -20.5, -19.5, -18.5, -17.5, -16.5, -15.5, -14.5, -13.5, -12.5,
      -11.5, -10.5,  -9.5,  -8.5,  -7.5,  -6.5,  -5.5,  -4.5,  -3.5,  -2.5,  -1.5,  -0.5,
        0.5,   1.5,   2.5,   3.5,   4.5,   5.5,   6.5,   7.5,   8.5,   9.5,  10.5,  11.5,
       12.5,  13.5,  14.5,  15.5,  16.5,  17.5,  18.5,  19.5,  20.5,  21.5,  22.5,  23.5,
       24.5,  25.5,  26.5,  27.5,  28.5,  29.5,  30.5,  31.5,  32.5,  33.5,  34.5,  35.5, 
       36.5,  37.5,  38.5,  39.5,  40.5,  41.5,  42.5,  43.5,  44.5,  45.5,  46.5,  47.5,
       48.5,  49.5,  50.5,  51.5,  52.5,  53.5,  54.5,  55.5,  56.5,  57.5,  58.5,  59.5,
       60.5,  61.5,  62.5,  63.5,  64.5,  65.5,  66.5,  67.5,  68.5,  69.5,  70.5,  71.5,
       72.5,  73.5,  74.5,  75.5,  76.5,  77.5,  78.5,  79.5,  80.5,  81.5,  82.5,  83.5,
       84.5,  85.5,  86.5,  87.5,  88.5,  89.5,  90.5,  91.5,  92.5,  93.5,  94.5, 95.5};

    //double Offset_XBPF022702 = 2.189;
    double Offset_XBPF022702 = 2.17325;
    //double Offset_XBPF022702 = 2.728;
    std::vector<double> FibreToCoord_XBPF022702 = {
      -95.5, -94.5, -93.5, -92.5, -91.5, -90.5, -89.5, -88.5, -87.5, -86.5, -85.5, -84.5, 
      -83.5, -82.5, -81.5, -80.5, -79.5, -78.5, -77.5, -76.5, -75.5, -74.5, -73.5, -72.5,
      -71.5, -70.5, -69.5, -68.5, -67.5, -66.5, -65.5, -64.5, -63.5, -62.5, -61.5, -60.5,
      -59.5, -58.5, -57.5, -56.5, -55.5, -54.5, -53.5, -52.5, -51.5, -50.5, -49.5, -48.5,
      -47.5, -46.5, -45.5, -44.5, -43.5, -42.5, -41.5, -40.5, -39.5, -38.5, -37.5, -36.5,
      -35.5, -34.5, -33.5, -32.5, -31.5, -30.5, -29.5, -28.5, -27.5, -26.5, -25.5, -24.5,
      -23.5, -22.5, -21.5, -20.5, -19.5, -18.5, -17.5, -16.5, -15.5, -14.5, -13.5, -12.5,
      -11.5, -10.5,  -9.5,  -8.5,  -7.5,  -6.5,  -5.5,  -4.5,  -3.5,  -2.5,  -1.5,  -0.5,
        0.5,   1.5,   2.5,   3.5,   4.5,   5.5,   6.5,   7.5,   8.5,   9.5,  10.5,  11.5,
       12.5,  13.5,  14.5,  15.5,  16.5,  17.5,  18.5,  19.5,  20.5,  21.5,  22.5,  23.5,
       24.5,  25.5,  26.5,  27.5,  28.5,  29.5,  30.5,  31.5,  32.5,  33.5,  34.5,  35.5, 
       36.5,  37.5,  38.5,  39.5,  40.5,  41.5,  42.5,  43.5,  44.5,  45.5,  46.5,  47.5,
       48.5,  49.5,  50.5,  51.5,  52.5,  53.5,  54.5,  55.5,  56.5,  57.5,  58.5,  59.5,
       60.5,  61.5,  62.5,  63.5,  64.5,  65.5,  66.5,  67.5,  68.5,  69.5,  70.5,  71.5,
       72.5,  73.5,  74.5,  75.5,  76.5,  77.5,  78.5,  79.5,  80.5,  81.5,  82.5,  83.5,
       84.5,  85.5,  86.5,  87.5,  88.5,  89.5,  90.5,  91.5,  92.5,  93.5,  94.5, 95.5};
};

#endif
