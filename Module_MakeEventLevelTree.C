#include "BeamLineClasses/class_BeamLine.h"
#include "BeamLineClasses/class_Spill.h"
#include "BeamLineClasses/unpack.h"
#include "BeamLineClasses/plotting.h"
#include "BeamLineClasses/Analyse.h"
#include <iostream>
#include <string>
#include <vector>
#include <inttypes.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TROOT.h>
#include <TInterpreter.h>


int main(int argc, char *argv[])
{
  if(argc<4)
  {
    std::cout << "USAGE: ./Module_MakeEventLevelTree.exe <INPUT FILE PATH, [Marcel's files of form SpillBySpill_*_raw.root]> <OUTPUT FILE APPENDAGE, [YYYY-MM-DD_HH:MM_YYY-MM-DD_HH:DD]> <OUTDIR, [WITH FINAL FORWARD SLASH]>" << std::endl; 
    std::exit(0);
  }

  std::string s_InputFile  = argv[1];
  std::string s_TimeInfo   = argv[2];
  std::string s_OutDir     = argv[3];

  gInterpreter->GenerateDictionary("vector<ULong64_t>","vector");

  std::map<std::pair<unsigned long long,unsigned long long>, Spill> map_Spills = matchVarsBySpill(s_InputFile);
  std::cout << "ALL MATCHING DONE! DUMPING TO EVENT TREE" << std::endl;
  dumpToEventTree(map_Spills, s_TimeInfo, s_OutDir);
  printSpills(map_Spills);

  return 0;
}
