#ifndef ANALYSE_H
#define ANALYSE_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "class_Spill.h"
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TKey.h>

std::map<std::pair<unsigned long long,unsigned long long>, Spill> matchVarsBySpill(std::string const &s_InputFile);
void printSpills    (std::map<std::pair<unsigned long long,unsigned long long>, Spill> &map_Spills);
void dumpToEventTree(std::map<std::pair<unsigned long long,unsigned long long>, Spill> &map_Spills, std::string const &s_TimeInfo, std::string const &s_OutDir);

#endif
