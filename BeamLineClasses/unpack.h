#ifndef UNPACK_H
#define UNPACK_H

#include "class_Detector.h"
#include <string>
#include <map>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>

void unpackDIPROOT(std::map<std::string,Detector> &cMapDetector, std::string const &cDataPath, bool const &cZeroSuppress);
void SetHistTitles(TH1I* &h, TString const &title, TString const &xTitle, TString const &yTitle);
void SetHistTitles(TH1D* &h, TString const &title, TString const &xTitle, TString const &yTitle);
void SetHistTitles(TH2D* &h, TString const &title, TString const &xTitle, TString const &yTitle);

#endif
