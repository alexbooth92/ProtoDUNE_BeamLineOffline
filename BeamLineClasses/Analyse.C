#include "Analyse.h"
#include "class_Spill.h"

unsigned long long convertTIAToUTC = 37;

std::pair<unsigned long long,unsigned long long> decompseXBPFTimeStamp(unsigned long long const &timeS, unsigned long long const &timeNS)
{
  std::pair<unsigned long long,unsigned long long> decomposedPair;

  unsigned long long MSB = (timeS - convertTIAToUTC)*1000000000;
  //unsigned long long LSB = timeNS*8;
  unsigned long long LSB = timeNS;
  MSB+=LSB;

  std::string s_AcqStamp = std::to_string(MSB);
  decomposedPair = {std::stoi(s_AcqStamp.substr(0,10)), std::stoi(s_AcqStamp.substr(10))};

  return decomposedPair;
}

std::map<std::pair<unsigned long long,unsigned long long>, Spill> unpackToDataTypes(std::string const &s_InputFile, 
                                                                                    std::vector<BasicDoubleVar> &vec_BasicDoubleVar, 
                                                                                    std::map<std::string,std::vector<BasicIntVar>> &map_BasicIntVar,
                                                                                    std::map<std::string,std::vector<VecBoolVar>>  &map_VecBoolVar, 
                                                                                    std::map<std::string,std::vector<VecLongVar>>  &map_VecLongVar)
{
  std::map<std::pair<unsigned long long,unsigned long long>, Spill> map_Spills;

  TFile *f_In = new TFile((TString)s_InputFile, "READ");
  TIter iter(f_In->GetListOfKeys());
  TKey *key;

  std::map<std::string,TTree*> map_NameToTrees;

  while((key=(TKey*)iter()))
  {
    TTree *t = (TTree*)key->ReadObj();
    map_NameToTrees[(std::string)t->GetName()] = t;  
  }

  std::vector<unsigned long long> testa, testb;
  for(auto tree : map_NameToTrees)
  {
    TIter iter_Branch(tree.second->GetListOfBranches());
    TKey *key_Branch;
    unsigned int  dataType            = 0;
    unsigned long long timeS          = 0;
    unsigned long long timeNS         = 0;
    int                basicIntVar    = 0;
    double             basicDoubleVar = 0.;
    std::vector<unsigned long long> *vecLongVar = 0;
    std::vector<bool>               *vecBoolVar = 0;

    tree.second->SetBranchAddress("timeS",  &timeS);
    tree.second->SetBranchAddress("timeNS", &timeNS);

    std::string s_BranchTitle = "";
    while((key_Branch=(TKey*)iter_Branch()))
    {
      s_BranchTitle = (std::string)key_Branch->GetTitle();
      if((s_BranchTitle.find("timeS") == std::string::npos) && (s_BranchTitle.find("timeNS") == std::string::npos))
      {
        if(s_BranchTitle.find("/D") != std::string::npos)
        {
          tree.second->SetBranchAddress((TString)s_BranchTitle.substr(0,s_BranchTitle.length()-2), &basicDoubleVar); 
          dataType = 1;
        }
        else if(s_BranchTitle.find("/I") != std::string::npos)
        {
          tree.second->SetBranchAddress((TString)(s_BranchTitle.substr(0,s_BranchTitle.length()-2)), &basicIntVar);
          dataType = 4;
        }
        else if(s_BranchTitle.find("profile") != std::string::npos)
        {
          tree.second->SetBranchAddress((TString)s_BranchTitle, &vecBoolVar); 
          dataType = 2;
        }
        else
        {
          tree.second->SetBranchAddress((TString)s_BranchTitle, &vecLongVar); 
          dataType = 3;
        }
      }
    }

    for(unsigned int i = 0; i < tree.second->GetEntries(); i++)
    //for(unsigned int i = 0; i < 5; i++)
    {
      tree.second->GetEntry(i);
      //std::cout << "THE TREE TITLE IS: " << tree.first << " THE TIME S IS: " << timeS << std::endl;
      if(dataType==1)
      {
        vec_BasicDoubleVar.push_back({tree.first, timeS, timeNS, basicDoubleVar});
      }
      else if(dataType==2)
      {
        std::pair<unsigned long long,unsigned long long> decomposedParts = decompseXBPFTimeStamp(timeS, timeNS);
        map_VecBoolVar[tree.first].push_back({tree.first, decomposedParts.first, decomposedParts.second, *vecBoolVar});
        vecBoolVar->clear();
      }
      else if(dataType==4)
      {
        map_BasicIntVar[tree.first].push_back({tree.first, timeS, timeNS, basicIntVar});
      }
      else if(dataType==3)
      {
        if(vecLongVar->size()!=0)
        {
          std::vector<unsigned long long> vecLongVarUTC = *vecLongVar;
          if(tree.first.find("SECONDS")!=std::string::npos)
          {
            std::for_each(vecLongVarUTC.begin(), vecLongVarUTC.end(), [](unsigned long long &seconds){seconds-=convertTIAToUTC;});
          }
          if(tree.first=="XBH4.GENERALTRIGGER:SECONDS")
          {
            //TIMESTAMP THE SPILL WITH THE FIRST GENERAL TRIGGER TIME IN UTC.
            Spill newSpill(vecLongVarUTC.at(0), 0);
            if(!(map_Spills.count({vecLongVarUTC.at(0), 0})))
            {
              map_Spills[{vecLongVarUTC.at(0), 0}] = newSpill;
              map_VecLongVar[tree.first].push_back({tree.first, timeS, timeNS, vecLongVarUTC});
            }
          }
          else
          {
            map_VecLongVar[tree.first].push_back({tree.first, timeS, timeNS, vecLongVarUTC});
          }
          vecLongVar->clear();
          vecLongVarUTC.clear();
        }
      }
    }
  }

  return map_Spills;
}


std::map<std::string,std::vector<std::vector<VecBoolVar>>> sortVecBoolVarToSpillsGroups(std::map<std::string,std::vector<VecBoolVar>> const &map_VecBoolVar)
{
  std::map<std::string,std::vector<std::vector<VecBoolVar>>> map_VecBoolVarSorted;

  for(auto it : map_VecBoolVar)
  {
    bool isNewSpill = false;
    std::vector<VecBoolVar> vec_VecBoolVar;
    std::vector<VecBoolVar> sortedByTime = it.second;
    std::sort(sortedByTime.begin(), sortedByTime.end(), [](const VecBoolVar &left, const VecBoolVar &right){return ((left.fTimeS*1000000000+left.fTimeNS)<(right.fTimeS*1000000000+right.fTimeNS));});
    for(unsigned int i = 0; i < it.second.size(); i++)
    {
      if(i==0 || isNewSpill==true)
      {
        vec_VecBoolVar.push_back(sortedByTime[i]);
        isNewSpill = false;
      }
      else
      {
        unsigned long long deltaSecs = sortedByTime[i].fTimeS - sortedByTime[i-1].fTimeS;
        if(deltaSecs < 6)
        {
          vec_VecBoolVar.push_back(sortedByTime[i]);
        }
        else
        {
          if(vec_VecBoolVar.size()>1)
          {
            map_VecBoolVarSorted[it.first].push_back(vec_VecBoolVar);
          }
          vec_VecBoolVar.clear();
          isNewSpill = true;
          i--;
        }

        if(i == sortedByTime.size()-1)
        {
          map_VecBoolVarSorted[it.first].push_back(vec_VecBoolVar);
          vec_VecBoolVar.clear();
          break;
        }
      }
    }
  }

  return map_VecBoolVarSorted;
}


std::map<std::string,std::map<std::string,std::vector<VecLongVar>>> combineVecLongVars(std::map<std::string,std::vector<VecLongVar>>  &map_VecLongVar,
                                                                                       std::map<std::string,std::vector<BasicIntVar>> &map_BasicIntVar)
{
  std::map<std::string,std::map<std::string,std::vector<VecLongVar>>> map_VecLongVarCombined;
  for(auto it : map_VecLongVar)
  {
    std::string s_Hardware = it.first.substr(0, it.first.find(":"));
    std::string s_Variable = it.first.substr(it.first.find(":") + 1);
    for(auto jt : map_BasicIntVar)
    {
      if(jt.first.find(s_Hardware)!=std::string::npos)
      {
        for(unsigned int i = 0; i < jt.second.size(); i++)
        {
          for(unsigned int j = 0; j < it.second.size(); j++)
          {
            if(jt.second[i].fTimeS==it.second[j].fTimeS)
            {
              map_VecLongVarCombined[s_Hardware][s_Variable].push_back(it.second[j]);
              if(s_Variable=="FRAC")
              {
                map_VecLongVarCombined[s_Hardware]["TIMESTAMPCOUNT"].push_back({jt.second[i].fBranchTitle, jt.second[i].fTimeS, jt.second[i].fTimeNS, {(unsigned long long)jt.second[i].fVar}}); 
              }
              break;
            }
          }
        }
      }
    }
  }

  return map_VecLongVarCombined;
}


void combineDataToSpill(std::map<std::pair<unsigned long long,unsigned long long>, Spill>   &map_Spills, std::vector<BasicDoubleVar> const &vec_BasicDoubleVar, 
                        std::map<std::string,std::vector<std::vector<VecBoolVar>>>          &map_VecBoolVarSorted,
                        std::map<std::string,std::map<std::string,std::vector<VecLongVar>>> &map_VecLongVarCombined)
{
  for(auto it : map_Spills)
  {
    for(unsigned int i = 0; i < vec_BasicDoubleVar.size(); i++)
    {
      map_Spills[it.first].addIfBetterMatch(vec_BasicDoubleVar[i]);
    }
    for(auto jt : map_VecBoolVarSorted)
    {
      for(unsigned int i = 0; i < jt.second.size(); i++)
      {
        map_Spills[it.first].addIfBetterMatch(map_VecBoolVarSorted[jt.first][i]);
      }
    }
    for(auto jt : map_VecLongVarCombined)
    {
      for(unsigned int i = 0; i < jt.second["TIMESTAMPCOUNT"].size(); i++)
      {
        map_Spills[it.first].addIfBetterMatch(jt.first, jt.second["FRAC"][i], jt.second["COARSE"][i], jt.second["SECONDS"][i], jt.second["TIMESTAMPCOUNT"][i]);
      }
    }
  }

  return;
}


void checkXBPFToGTConsistency(std::map<std::pair<unsigned long long,unsigned long long>, Spill> &map_Spills)
{
  for(auto it : map_Spills)
  {
    map_Spills[it.first].matchXBPFToGT();
  }

  return;
}


void matchXTDCEventWise(std::map<std::pair<unsigned long long,unsigned long long>, Spill> &map_Spills)
{
  for(auto it : map_Spills)
  {
    map_Spills[it.first].matchXTDCToGeneralTrigger();
    map_Spills[it.first].constructTimeOfFlight();
    //map_Spills[it.first].printTDCMatching();
  }

  return;
}


unsigned int rankEvent(std::vector<unsigned int> const &XBH4XBPF022697_NFibresHit,std::vector<unsigned int> const &XBH4XBPF022701_NFibresHit,std::vector<unsigned int> const &XBH4XBPF022702_NFibresHit, 
                       std::vector<unsigned int> const &XBH4XBPF022697_Span,      std::vector<unsigned int> const &XBH4XBPF022701_Span,     std::vector<unsigned int> const &XBH4XBPF022702_Span,
                       unsigned int const &nMatchedTOFs, unsigned int const &XBH4XTDC022713_nMatched, unsigned int const &XBH4XTDC022716_nMatched)
{
  unsigned int rank = 10;
  bool uniqueEventMatchedInXBPFSpectrometer = false;
  bool singleHitInAllXBPFSpectrometer       = false;
  bool uniqueButAdjacentFibres              = false;
  bool uniqueTOFChannel                     = false;
  bool uniqueXCETMatched                    = false;

  if(XBH4XBPF022697_NFibresHit.size()==1 && XBH4XBPF022701_NFibresHit.size()==1 && XBH4XBPF022702_NFibresHit.size()==1)
  { 
    uniqueEventMatchedInXBPFSpectrometer = true; 
    if(XBH4XBPF022697_NFibresHit[0]==1 && XBH4XBPF022701_NFibresHit[0]==1 && XBH4XBPF022702_NFibresHit[0]==1)
    {
      singleHitInAllXBPFSpectrometer = true;
    }
    else if(uniqueEventMatchedInXBPFSpectrometer && ((XBH4XBPF022697_NFibresHit[0]<=2 && XBH4XBPF022697_Span[0]<=1) &&
                                                     (XBH4XBPF022701_NFibresHit[0]<=2 && XBH4XBPF022701_Span[0]<=1) &&
                                                     (XBH4XBPF022702_NFibresHit[0]<=2 && XBH4XBPF022702_Span[0]<=1)))
    {
      uniqueButAdjacentFibres = true;
    }
  }
  if(nMatchedTOFs==1){ uniqueTOFChannel = true; }
  if(XBH4XTDC022713_nMatched<=1 && XBH4XTDC022716_nMatched<=1){ uniqueXCETMatched = true; }

  if(uniqueEventMatchedInXBPFSpectrometer && singleHitInAllXBPFSpectrometer && uniqueTOFChannel && uniqueXCETMatched)
  {
    rank = 1;
  }
  else if(uniqueEventMatchedInXBPFSpectrometer && uniqueButAdjacentFibres && uniqueTOFChannel && uniqueXCETMatched)
  {
    rank = 2;
  }
  else if(uniqueEventMatchedInXBPFSpectrometer && singleHitInAllXBPFSpectrometer)
  {
    rank = 3;
  }

  return rank;
}


std::map<std::pair<unsigned long long,unsigned long long>, Spill> matchVarsBySpill(std::string const &s_InputFile)
{
  std::vector<BasicDoubleVar> vec_BasicDoubleVar;
  std::map<std::string,std::vector<BasicIntVar>> map_BasicIntVar;
  std::map<std::string,std::vector<VecBoolVar>>  map_VecBoolVar;
  std::map<std::string,std::vector<VecLongVar>>  map_VecLongVar;

  std::cout << "UNPACKING TO DATA TYPES" << std::endl;
  std::map<std::pair<unsigned long long,unsigned long long>, Spill> map_Spills = unpackToDataTypes(s_InputFile, vec_BasicDoubleVar, map_BasicIntVar, map_VecBoolVar, map_VecLongVar);

  std::cout << "SORTING VECBOOLVARS TO INTO SPILL GROUPS" << std::endl;
  std::map<std::string,std::vector<std::vector<VecBoolVar>>> map_VecBoolVarSorted = sortVecBoolVarToSpillsGroups(map_VecBoolVar);

  std::cout << "COMBINING VECLONGVAR VARIABLES" << std::endl;
  std::map<std::string,std::map<std::string,std::vector<VecLongVar>>> map_VecLongVarCombined = combineVecLongVars(map_VecLongVar, map_BasicIntVar);

  std::cout << "COMBINING DATA TYPES WITHIN SPILL CLASS OBJECTS" << std::endl;
  combineDataToSpill(map_Spills, vec_BasicDoubleVar, map_VecBoolVarSorted, map_VecLongVarCombined);

  std::cout << "CHECKING XBPFs ARE CONSISTENT WITH GENERAL TRIGGERS" << std::endl;
  checkXBPFToGTConsistency(map_Spills);

  std::cout << "MATCHING EVENT-WISE XCETs AND XTOFs TO GENERAL TRIGGERS" << std::endl;
  matchXTDCEventWise(map_Spills);

  return map_Spills;
}


void printSpills(std::map<std::pair<unsigned long long,unsigned long long>, Spill> &map_Spills)
{
  for(auto it : map_Spills)
  {
    //it.second.printSpill();
    it.second.printSpillConcise();
  }

  return;
}


std::string getPID(double const &tof, bool const &XBH4XTDC022713_On, bool const &XBH4XTDC022716_On, bool const &C2IsHighPressure, 
                   double const &refMomentum)
{
  std::string PID = "nopid";

  bool statusCL = 0; bool statusCH = 0;
  if(C2IsHighPressure){
    statusCL = XBH4XTDC022713_On;
    statusCH = XBH4XTDC022716_On;
  }
  else{    
    statusCH = XBH4XTDC022713_On;
    statusCL = XBH4XTDC022716_On;
  }

  if(refMomentum==1){
    if     (tof < 105 && statusCH) {PID =     "e";}
    else if(tof < 110 && !statusCH){PID = "mu/pi";}
    else if(tof > 110 && tof < 160 && !statusCH){PID = "P";}
  }
  else if(refMomentum==2){
    if     (tof < 105 && statusCH) {PID = "e";}
    else if(tof < 103 && !statusCH){PID = "mu/pi";}
    else if(tof > 103 && tof < 160 && !statusCH){PID = "P";}
  }
  else if(refMomentum==3){
    if     (statusCL==true  && statusCH==true) {PID = "e";    }
    else if(statusCL==false && statusCH==true) {PID = "mu/pi";}
    else if(statusCL==false && statusCH==false){PID = "K/P";  }
  }
  else if(refMomentum==6 || refMomentum==7){
    if     (statusCL==true  && statusCH==true) {PID = "e/mu/pi";}
    else if(statusCL==false && statusCH==true) {PID = "K";      }
    else if(statusCL==false && statusCH==false){PID = "P";      }
  }

  return PID;
}


void assignValidTOFAndPID(double &validTOF, std::string &validChannel, std::string &PID, 
                          std::vector<double> const &tof, std::vector<std::string> const &tofChannel,
                          std::vector<unsigned long long> const &XBH4_XTDC_022_713_TimestampNS, 
                          std::vector<unsigned long long> const &XBH4_XTDC_022_716_TimestampNS, 
                          double const &refMomentum, bool const C2IsHighPressure)
{
  validTOF     = -99.;
  validChannel = "nouniquietof";
  PID          = "nopid";

  std::vector<double> vAA;
  std::vector<double> vAB;
  std::vector<double> vBA;
  std::vector<double> vBB;

  for(unsigned int i = 0; i < tof.size(); i++)
  {
    if     (tofChannel.at(i)=="AA"){vAA.push_back(tof.at(i));}
    else if(tofChannel.at(i)=="AB"){vAB.push_back(tof.at(i));}
    else if(tofChannel.at(i)=="BA"){vBA.push_back(tof.at(i));}
    else if(tofChannel.at(i)=="BB"){vBB.push_back(tof.at(i));}
  }

  if(vAA.size()+vAB.size()+vBA.size()+vBB.size()>0)
  {
    int statusCL = 0; int statusCH = 0;
    if(C2IsHighPressure){
      if(XBH4_XTDC_022_713_TimestampNS.size()>0) statusCL=1;
      if(XBH4_XTDC_022_716_TimestampNS.size()>0) statusCH=1;
    }
    else{
      if(XBH4_XTDC_022_713_TimestampNS.size()>0) statusCH=1;
      if(XBH4_XTDC_022_716_TimestampNS.size()>0) statusCL=1;
    }
    double mini = 1e15;
    double maxi = 0;
    std::string validChannel_Mini = "";
    std::string validChannel_Maxi = "";
    for(int i=0; i < vAA.size(); i++){
      if(mini > vAA[i]){
        mini = vAA[i];
        validChannel_Mini = "AA";
      }
      if(maxi < vAA[i]){
        maxi = vAA[i];
        validChannel_Maxi = "AA";
      }
    }
    for(int i=0; i < vAB.size(); i++){
      if(mini > vAB[i]){
        mini = vAB[i];
        validChannel_Mini = "AB";
      }
      if(maxi < vAB[i]){
        maxi = vAB[i];
        validChannel_Maxi = "AB";
      }
    }
    for(int i=0; i < vBA.size(); i++){
      if(mini > vBA[i]){
        mini = vBA[i];
        validChannel_Mini = "BA";
      }
      if(maxi < vBA[i]){
        maxi = vBA[i];
        validChannel_Maxi = "BA";
      }
    }
    for(int i=0; i < vBB.size(); i++){
      if(mini > vBB[i]){
        mini = vBB[i];
        validChannel_Mini = "BB";
      }
      if(maxi < vBB[i]){
        maxi = vBB[i];
        validChannel_Maxi = "BB";
      }
    }

    if(refMomentum==1 || refMomentum==2){
      PID = getPID(mini, XBH4_XTDC_022_713_TimestampNS.size() > 0 ? true : false, 
                         XBH4_XTDC_022_716_TimestampNS.size() > 0 ? true : false, C2IsHighPressure, refMomentum);
      if     (PID=="e")    {validTOF = mini; validChannel = validChannel_Mini;}
      else if(PID=="mu/pi"){validTOF = mini; validChannel = validChannel_Mini;}
      else if(getPID(maxi, XBH4_XTDC_022_713_TimestampNS.size() > 0 ? true : false, 
                           XBH4_XTDC_022_716_TimestampNS.size() > 0 ? true : false, C2IsHighPressure, refMomentum)=="P"){
        validChannel = validChannel_Maxi; 
        validTOF     = maxi;    
        PID          = "P"; 
      }
    }
    else if(maxi-mini < 10){
      validChannel = validChannel_Mini;
      validTOF     = mini;
      PID          = getPID(mini, XBH4_XTDC_022_713_TimestampNS.size() > 0 ? true : false, 
                                  XBH4_XTDC_022_716_TimestampNS.size() > 0 ? true : false, C2IsHighPressure, refMomentum);
    }
  }

  return;
}


void dumpToEventTree(std::map<std::pair<unsigned long long,unsigned long long>, Spill> &map_Spills, std::string const &s_TimeInfo, std::string const &s_OutDir)
{
  TFile *f_Out = new TFile((TString)s_OutDir+"EventTree_"+(TString)s_TimeInfo+".root", "RECREATE");
  TTree *t_Out = new TTree("EventTree","EventTree");

  unsigned int eventRank;
  unsigned long long spillTimeStamp; //UTC UNIX TIMESTAMP OF FIRST GENERAL TRIGGER IN THE SPILL THE EVENT IS ASSOCIATED WITH. UNIQUE SPILL IDENTIFIER.
  double fractionComplete;
  double referenceMomentum;
  double XBH4BEND022692I_MEAS;	
  double XBH4BEND022699I_MEAS;	
  double XBH4EXPTNP04001COUNTS;
  double XBH4EXPTNP04002COUNTS;
  double XBH4EXPTNP04003COUNTS;
  double XBH4EXPTNP04004COUNTS;
  double XBH4EXPTNP04009COUNTS;
  double XBH4EXPTNP04010COUNTS;
  double XBH4EXPTNP04011COUNTS;
  double XBH4EXPTNP04012COUNTS;
  double XBH4XCET022713COUNTS_TRIG;
  double XBH4XCET022713PRESSURE;
  double XBH4XCET022713SIMPLE_COUNTS;
  double XBH4XCET022716COUNTS_TRIG;
  double XBH4XCET022716PRESSURE;
  double XBH4XCET022716SIMPLE_COUNTS;
  double XBH4XCSH022694POS_JAW1_MEAS;
  double XBH4XCSH022694POS_JAW2_MEAS;
  double XBH4XSCI022680COUNTS; 
  double SPST2INTENSITY;
  
  unsigned int XBH4GENERALTRIGGER_NCounts; 
  unsigned int XBH4XTDC022713_NCounts;
  unsigned int XBH4XTDC022716_NCounts;
  unsigned int XBTF022687A_NCounts; 
  unsigned int XBTF022687B_NCounts; 
  unsigned int XBTF022716A_NCounts; 
  unsigned int XBTF022716B_NCounts;

  //TIMESTAMPNS IS THE NANOSCOND PIECE OF THE TIMESTAMP WITHOUT THE FRAC CLOCK. TIMESTAMPFRACACCURACY IS THE NS AND FRACTIONAL NS CONTRIBUTION TO THE TIMESTAMP FROM FRAC.
  //TO GET THE TIMESTAMP TO THE MAXIMUM ACCURACY YOU NEED TO SIMPLY ADD TIMESTAMPNS AND TIMESTAMPFRACACCURACY, NOT DONE HERE BECAUSE YOU
  //CANT FIT THAT MANY DIGITS INTO A LONG LONG.
  unsigned long long XBH4GENERALTRIGGER_TimestampNS;          double XBH4GENERALTRIGGER_TimestampFracAccuracy; 
  std::vector<unsigned long long> XBH4XTDC022713_TimestampNS; std::vector<double> XBH4XTDC022713_TimestampFracAccuracy;
  std::vector<unsigned long long> XBH4XTDC022716_TimestampNS; std::vector<double> XBH4XTDC022716_TimestampFracAccuracy;
  std::vector<unsigned long long> XBTF022687A_TimestampNS;    std::vector<double> XBTF022687A_TimestampFracAccuracy;
  std::vector<unsigned long long> XBTF022687B_TimestampNS;    std::vector<double> XBTF022687B_TimestampFracAccuracy;
  std::vector<unsigned long long> XBTF022716A_TimestampNS;    std::vector<double> XBTF022716A_TimestampFracAccuracy;
  std::vector<unsigned long long> XBTF022716B_TimestampNS;    std::vector<double> XBTF022716B_TimestampFracAccuracy;

  std::vector<std::string> tofChannel;
  std::vector<double> tof;
  bool XBH4XTDC022713_On;
  bool XBH4XTDC022716_On;

  std::vector<unsigned long long> XBH4XBPF022697_TimestampNS; std::vector<unsigned long long> XBH4XBPF022701_TimestampNS; 
  std::vector<unsigned long long> XBH4XBPF022707_TimestampNS; std::vector<unsigned long long> XBH4XBPF022716_TimestampNS;
  std::vector<unsigned long long> XBH4XBPF022698_TimestampNS; std::vector<unsigned long long> XBH4XBPF022702_TimestampNS;    
  std::vector<unsigned long long> XBH4XBPF022708_TimestampNS; std::vector<unsigned long long> XBH4XBPF022717_TimestampNS;

  std::vector<unsigned int> XBH4XBPF022697_NFibresHit; std::vector<unsigned int> XBH4XBPF022701_NFibresHit; 
  std::vector<unsigned int> XBH4XBPF022707_NFibresHit; std::vector<unsigned int> XBH4XBPF022716_NFibresHit;
  std::vector<unsigned int> XBH4XBPF022698_NFibresHit; std::vector<unsigned int> XBH4XBPF022702_NFibresHit;     
  std::vector<unsigned int> XBH4XBPF022708_NFibresHit; std::vector<unsigned int> XBH4XBPF022717_NFibresHit;

  std::vector<unsigned int> XBH4XBPF022697_Span; std::vector<unsigned int> XBH4XBPF022701_Span;           
  std::vector<unsigned int> XBH4XBPF022707_Span; std::vector<unsigned int> XBH4XBPF022716_Span;
  std::vector<unsigned int> XBH4XBPF022698_Span; std::vector<unsigned int> XBH4XBPF022702_Span;           
  std::vector<unsigned int> XBH4XBPF022708_Span; std::vector<unsigned int> XBH4XBPF022717_Span;

  std::vector<std::vector<double>> XBH4XBPF022697_HitCoordinates; std::vector<std::vector<double>> XBH4XBPF022701_HitCoordinates; 
  std::vector<std::vector<double>> XBH4XBPF022707_HitCoordinates; std::vector<std::vector<double>> XBH4XBPF022716_HitCoordinates;
  std::vector<std::vector<double>> XBH4XBPF022698_HitCoordinates; std::vector<std::vector<double>> XBH4XBPF022702_HitCoordinates; 
  std::vector<std::vector<double>> XBH4XBPF022708_HitCoordinates; std::vector<std::vector<double>> XBH4XBPF022717_HitCoordinates;

  std::vector<unsigned int>        nMomenta;
  std::vector<std::vector<double>> deflectionAngle;
  std::vector<std::vector<double>> reconstructedMomentum;

  double      validTOF;
  std::string validChannel;
  std::string PID;

  t_Out->Branch("eventRank",                       &eventRank,                   "eventRank/I"                      );
  t_Out->Branch("spillTimeStamp",                  &spillTimeStamp,              "spillTimeStamp/l"                 );
  t_Out->Branch("fractionComplete",                &fractionComplete,            "fractionComplete/D"               ); 
  t_Out->Branch("referenceMomentum",               &referenceMomentum,           "referenceMomentum/D"              );
  t_Out->Branch("XBH4.BEND.022.692_I_MEAS",        &XBH4BEND022692I_MEAS,        "XBH4.BEND.022.692_I_MEAS/D"       );	
  t_Out->Branch("XBH4.BEND.022.699_I_MEAS",        &XBH4BEND022699I_MEAS,        "XBH4.BEND.022.699_I_MEAS/D"       );	
  t_Out->Branch("XBH4.EXPT.NP04.001_COUNTS",       &XBH4EXPTNP04001COUNTS,       "XBH4.EXPT.NP04.001_COUNTS/D"      );
  t_Out->Branch("XBH4.EXPT.NP04.002_COUNTS",       &XBH4EXPTNP04002COUNTS,       "XBH4.EXPT.NP04.002_COUNTS/D"      );
  t_Out->Branch("XBH4.EXPT.NP04.003_COUNTS",       &XBH4EXPTNP04003COUNTS,       "XBH4.EXPT.NP04.003_COUNTS/D"      );
  t_Out->Branch("XBH4.EXPT.NP04.004_COUNTS",       &XBH4EXPTNP04004COUNTS,       "XBH4.EXPT.NP04.004_COUNTS/D"      );
  t_Out->Branch("XBH4.EXPT.NP04.009_COUNTS",       &XBH4EXPTNP04009COUNTS,       "XBH4.EXPT.NP04.009_COUNTS/D"      );
  t_Out->Branch("XBH4.EXPT.NP04.010_COUNTS",       &XBH4EXPTNP04010COUNTS,       "XBH4.EXPT.NP04.010_COUNTS/D"      );
  t_Out->Branch("XBH4.EXPT.NP04.011_COUNTS",       &XBH4EXPTNP04011COUNTS,       "XBH4.EXPT.NP04.011_COUNTS/D"      );
  t_Out->Branch("XBH4.EXPT.NP04.012_COUNTS",       &XBH4EXPTNP04012COUNTS,       "XBH4.EXPT.NP04.012_COUNTS/D"      );
  t_Out->Branch("XBH4.XCET.022.713_COUNTS_TRIG",   &XBH4XCET022713COUNTS_TRIG,   "XBH4.XCET.022.713_COUNTS_TRIG/D"  );
  t_Out->Branch("XBH4.XCET.022.713_PRESSURE",      &XBH4XCET022713PRESSURE,      "XBH4.XCET.022.713_PRESSURE/D"     );
  t_Out->Branch("XBH4.XCET.022.713_SIMPLE_COUNTS", &XBH4XCET022713SIMPLE_COUNTS, "XBH4.XCET.022.713_SIMPLE_COUNTS/D");
  t_Out->Branch("XBH4.XCET.022.716_COUNTS_TRIG",   &XBH4XCET022716COUNTS_TRIG,   "XBH4.XCET.022.716_COUNTS_TRIG/D"  );
  t_Out->Branch("XBH4.XCET.022.716_PRESSURE",      &XBH4XCET022716PRESSURE,      "XBH4.XCET.022.716_PRESSURE/D"     );
  t_Out->Branch("XBH4.XCET.022.716_SIMPLE_COUNTS", &XBH4XCET022716SIMPLE_COUNTS, "XBH4.XCET.022.716_SIMPLE_COUNTS/D");
  t_Out->Branch("XBH4.XCSH.022.694_POS_JAW1_MEAS", &XBH4XCSH022694POS_JAW1_MEAS, "XBH4.XCSH.022.694_POS_JAW1_MEAS/D");
  t_Out->Branch("XBH4.XCSH.022.694_POS_JAW2_MEAS", &XBH4XCSH022694POS_JAW2_MEAS, "XBH4.XCSH.022.694_POS_JAW2_MEAS/D");
  t_Out->Branch("XBH4.XSCI.022.680_COUNTS",        &XBH4XSCI022680COUNTS,        "XBH4.XSCI.022.680_COUNTS/D"       ); 
  t_Out->Branch("SPS.T2_INTENSITY",                &SPST2INTENSITY,              "XBH4.XTDC.022.713/D"              );

  t_Out->Branch("XBH4GENERALTRIGGER_NCounts", &XBH4GENERALTRIGGER_NCounts, "XBH4GENERALTRIGGER_NCounts/I"); 
  t_Out->Branch("XBH4XTDC022713_NCounts",     &XBH4XTDC022713_NCounts,     "XBH4XTDC022713_NCounts/I"    );
  t_Out->Branch("XBH4XTDC022716_NCounts",     &XBH4XTDC022716_NCounts,     "XBH4XTDC022716_NCounts/I"    );
  t_Out->Branch("XBTF022687A_NCounts",        &XBTF022687A_NCounts,        "XBTF022687A_NCounts/I"       ); 
  t_Out->Branch("XBTF022687B_NCounts",        &XBTF022687B_NCounts,        "XBTF022687B_NCounts/I"       ); 
  t_Out->Branch("XBTF022716A_NCounts",        &XBTF022716A_NCounts,        "XBTF022716A_NCounts/I"       ); 
  t_Out->Branch("XBTF022716B_NCounts",        &XBTF022716B_NCounts,        "XBTF022716B_NCounts/I"       );

  t_Out->Branch("XBH4GENERALTRIGGER_TimestampNS", &XBH4GENERALTRIGGER_TimestampNS, "XBH4GENERALTRIGGER_TimestampNS/l"); 
  t_Out->Branch("XBH4GENERALTRIGGER_TimestampFracAccuracy", &XBH4GENERALTRIGGER_TimestampFracAccuracy, "XBH4GENERALTRIGGER_TimestampFracAccuracy/D"); 
  t_Out->Branch("XBH4.XTDC.022.713_TimestampNS",  &XBH4XTDC022713_TimestampNS);
  t_Out->Branch("XBH4.XTDC.022.713_TimestampFracAccuracy",  &XBH4XTDC022713_TimestampFracAccuracy);
  t_Out->Branch("XBH4.XTDC.022.716_TimestampNS",  &XBH4XTDC022716_TimestampNS);
  t_Out->Branch("XBH4.XTDC.022.716_TimestampFracAccuracy",  &XBH4XTDC022716_TimestampFracAccuracy);
  t_Out->Branch("XBTF.022.687A_TimestampNS",      &XBTF022687A_TimestampNS   );
  t_Out->Branch("XBTF.022.687A_TimestampFracAccuracy",      &XBTF022687A_TimestampFracAccuracy   );
  t_Out->Branch("XBTF.022.687B_TimestampNS",      &XBTF022687B_TimestampNS   );
  t_Out->Branch("XBTF.022.687B_TimestampFracAccuracy",      &XBTF022687B_TimestampFracAccuracy   );
  t_Out->Branch("XBTF.022.716A_TimestampNS",      &XBTF022716A_TimestampNS   );
  t_Out->Branch("XBTF.022.716A_TimestampFracAccuracy",      &XBTF022716A_TimestampFracAccuracy   );
  t_Out->Branch("XBTF.022.716B_TimestampNS",      &XBTF022716B_TimestampNS   );
  t_Out->Branch("XBTF.022.716B_TimestampFracAccuracy",      &XBTF022716B_TimestampFracAccuracy   );

  t_Out->Branch("TOFChannel", &tofChannel);
  t_Out->Branch("TOF",        &tof);
  t_Out->Branch("XBH4.XTDC.022.713_On", &XBH4XTDC022713_On, "XBH4.XTDC.022.713_On/B"); 
  t_Out->Branch("XBH4.XTDC.022.716_On", &XBH4XTDC022716_On, "XBH4.XTDC.022.716_On/B"); 

  t_Out->Branch("XBH4.XBPF.022.697_TimestampNS",    &XBH4XBPF022697_TimestampNS   );  
  t_Out->Branch("XBH4.XBPF.022.697_NFibresHit",     &XBH4XBPF022697_NFibresHit    );   
  t_Out->Branch("XBH4.XBPF.022.697_Span",           &XBH4XBPF022697_Span          );         
  t_Out->Branch("XBH4.XBPF.022.697_HitCoordinates", &XBH4XBPF022697_HitCoordinates); 
  t_Out->Branch("XBH4.XBPF.022.698_TimestampNS",    &XBH4XBPF022698_TimestampNS   );  
  t_Out->Branch("XBH4.XBPF.022.698_NFibresHit",     &XBH4XBPF022698_NFibresHit    );   
  t_Out->Branch("XBH4.XBPF.022.698_Span",           &XBH4XBPF022698_Span          );         
  t_Out->Branch("XBH4.XBPF.022.698_HitCoordinates", &XBH4XBPF022698_HitCoordinates); 
  t_Out->Branch("XBH4.XBPF.022.701_TimestampNS",    &XBH4XBPF022701_TimestampNS   );  
  t_Out->Branch("XBH4.XBPF.022.701_NFibresHit",     &XBH4XBPF022701_NFibresHit    );   
  t_Out->Branch("XBH4.XBPF.022.701_Span",           &XBH4XBPF022701_Span          );         
  t_Out->Branch("XBH4.XBPF.022.701_HitCoordinates", &XBH4XBPF022701_HitCoordinates); 
  t_Out->Branch("XBH4.XBPF.022.702_TimestampNS",    &XBH4XBPF022702_TimestampNS   );  
  t_Out->Branch("XBH4.XBPF.022.702_NFibresHit",     &XBH4XBPF022702_NFibresHit    );   
  t_Out->Branch("XBH4.XBPF.022.702_Span",           &XBH4XBPF022702_Span          );         
  t_Out->Branch("XBH4.XBPF.022.702_HitCoordinates", &XBH4XBPF022702_HitCoordinates); 
  t_Out->Branch("XBH4.XBPF.022.707_TimestampNS",    &XBH4XBPF022707_TimestampNS   );  
  t_Out->Branch("XBH4.XBPF.022.707_NFibresHit",     &XBH4XBPF022707_NFibresHit    );   
  t_Out->Branch("XBH4.XBPF.022.707_Span",           &XBH4XBPF022707_Span          );         
  t_Out->Branch("XBH4.XBPF.022.707_HitCoordinates", &XBH4XBPF022707_HitCoordinates); 
  t_Out->Branch("XBH4.XBPF.022.708_TimestampNS",    &XBH4XBPF022708_TimestampNS   );  
  t_Out->Branch("XBH4.XBPF.022.708_NFibresHit",     &XBH4XBPF022708_NFibresHit    );   
  t_Out->Branch("XBH4.XBPF.022.708_Span",           &XBH4XBPF022708_Span          );         
  t_Out->Branch("XBH4.XBPF.022.708_HitCoordinates", &XBH4XBPF022708_HitCoordinates); 
  t_Out->Branch("XBH4.XBPF.022.716_TimestampNS",    &XBH4XBPF022716_TimestampNS   );  
  t_Out->Branch("XBH4.XBPF.022.716_NFibresHit",     &XBH4XBPF022716_NFibresHit    );   
  t_Out->Branch("XBH4.XBPF.022.716_Span",           &XBH4XBPF022716_Span          );         
  t_Out->Branch("XBH4.XBPF.022.716_HitCoordinates", &XBH4XBPF022716_HitCoordinates); 
  t_Out->Branch("XBH4.XBPF.022.717_TimestampNS",    &XBH4XBPF022717_TimestampNS   );  
  t_Out->Branch("XBH4.XBPF.022.717_NFibresHit",     &XBH4XBPF022717_NFibresHit    );   
  t_Out->Branch("XBH4.XBPF.022.717_Span",           &XBH4XBPF022717_Span          );         
  t_Out->Branch("XBH4.XBPF.022.717_HitCoordinates", &XBH4XBPF022717_HitCoordinates); 

  t_Out->Branch("nPossibleMomenta",      &nMomenta             );
  t_Out->Branch("deflectionAngle",       &deflectionAngle      );
  t_Out->Branch("reconstructedMomentum", &reconstructedMomentum);

  t_Out->Branch("validTOF",        &validTOF, "validTOF/D");
  t_Out->Branch("validTOFChannel", &validChannel          );
  t_Out->Branch("PID",             &PID                   );

  for(auto it : map_Spills)
  {
    std::map<std::string,BasicDoubleVar>  mapBasicDoubleVars = it.second.getMapBasicDoubleVars();
    std::map<std::string,AcquisitionXBTF> mapTDCAcquisitions = it.second.getMapTDCAcquisitions();
    std::map<std::pair<unsigned int,std::string>,std::vector<unsigned int>> mapGenTrigEventsToTDCEvents = it.second.getMapGenTrigToTDCEvents();
    std::map<unsigned int,std::vector<TOF>> mapGenTrigEventToTOFs = it.second.getMapGenTrigToEventTOFs();
    std::map<unsigned int,std::vector<XBPFCoincidence>> mapVecXBPFCoincidence = it.second.getMapGenTrigToXBPFCoincidence();

    std::vector<AcquisitionXBTF::EventRecordHR> GenTrigRecords          = mapTDCAcquisitions["XBH4.GENERALTRIGGER"].getDataHR();
    std::vector<AcquisitionXBTF::EventRecordHR> XBH4XTDC022713_Records  = mapTDCAcquisitions["XBH4.XTDC.022.713"]  .getDataHR();
    std::vector<AcquisitionXBTF::EventRecordHR> XBH4XTDC022716_Records  = mapTDCAcquisitions["XBH4.XTDC.022.716"]  .getDataHR();
    std::vector<AcquisitionXBTF::EventRecordHR> XBTF022687A_Records     = mapTDCAcquisitions["XBTF.022.687A"]      .getDataHR(); 
    std::vector<AcquisitionXBTF::EventRecordHR> XBTF022687B_Records     = mapTDCAcquisitions["XBTF.022.687B"]      .getDataHR(); 
    std::vector<AcquisitionXBTF::EventRecordHR> XBTF022716A_Records     = mapTDCAcquisitions["XBTF.022.716A"]      .getDataHR(); 
    std::vector<AcquisitionXBTF::EventRecordHR> XBTF022716B_Records     = mapTDCAcquisitions["XBTF.022.716B"]      .getDataHR();

    for(unsigned int i = 0; i < mapTDCAcquisitions["XBH4.GENERALTRIGGER"].getTimestampCount(); i++)
    {
      XBH4GENERALTRIGGER_NCounts = 0; 
      XBH4XTDC022713_NCounts     = 0;
      XBH4XTDC022716_NCounts     = 0;
      XBTF022687A_NCounts        = 0; 
      XBTF022687B_NCounts        = 0; 
      XBTF022716A_NCounts        = 0; 
      XBTF022716B_NCounts        = 0;

      XBH4XTDC022713_TimestampNS.clear(); XBH4XTDC022713_TimestampFracAccuracy.clear();
      XBH4XTDC022716_TimestampNS.clear(); XBH4XTDC022716_TimestampFracAccuracy.clear();
      XBTF022687A_TimestampNS.clear();    XBTF022687A_TimestampFracAccuracy.clear();
      XBTF022687B_TimestampNS.clear();    XBTF022687B_TimestampFracAccuracy.clear();
      XBTF022716A_TimestampNS.clear();    XBTF022716A_TimestampFracAccuracy.clear();
      XBTF022716B_TimestampNS.clear();    XBTF022716B_TimestampFracAccuracy.clear();
      //DEFAULT FALSE, IF THERE IS A SIGNAL FROM THEM, THIS IS REFLECTED IN THE LOOPS BELOW.
      XBH4XTDC022713_On = false;
      XBH4XTDC022716_On = false;

      tofChannel.clear();
      tof.clear();
      XBH4XBPF022697_HitCoordinates.clear(); XBH4XBPF022701_HitCoordinates.clear(); XBH4XBPF022707_HitCoordinates.clear(); XBH4XBPF022716_HitCoordinates.clear();
      XBH4XBPF022698_HitCoordinates.clear(); XBH4XBPF022702_HitCoordinates.clear(); XBH4XBPF022708_HitCoordinates.clear(); XBH4XBPF022717_HitCoordinates.clear();
      deflectionAngle.clear();
      reconstructedMomentum.clear();

      XBH4XBPF022697_TimestampNS.clear(); XBH4XBPF022701_TimestampNS.clear(); XBH4XBPF022707_TimestampNS.clear(); XBH4XBPF022716_TimestampNS.clear();
      XBH4XBPF022697_NFibresHit.clear();  XBH4XBPF022701_NFibresHit.clear();  XBH4XBPF022707_NFibresHit.clear();  XBH4XBPF022716_NFibresHit.clear();
      XBH4XBPF022697_Span.clear();        XBH4XBPF022701_Span.clear();        XBH4XBPF022707_Span.clear();        XBH4XBPF022716_Span.clear();
      XBH4XBPF022698_TimestampNS.clear(); XBH4XBPF022702_TimestampNS.clear(); XBH4XBPF022708_TimestampNS.clear(); XBH4XBPF022717_TimestampNS.clear();
      XBH4XBPF022698_NFibresHit.clear();  XBH4XBPF022702_NFibresHit.clear();  XBH4XBPF022708_NFibresHit.clear();  XBH4XBPF022717_NFibresHit.clear();
      XBH4XBPF022698_Span.clear();        XBH4XBPF022702_Span.clear();        XBH4XBPF022708_Span.clear();        XBH4XBPF022717_Span.clear();
      nMomenta.clear();

      spillTimeStamp              = it.first.first;     
      fractionComplete            = it.second.getFractionComplete();
      referenceMomentum           = it.second.getReferenceMomemtum();
      XBH4BEND022692I_MEAS        = mapBasicDoubleVars["XBH4.BEND.022.692:I_MEAS"       ].fVar;	
      XBH4BEND022699I_MEAS        = mapBasicDoubleVars["XBH4.BEND.022.699:I_MEAS"       ].fVar;	
      XBH4EXPTNP04001COUNTS       = mapBasicDoubleVars["XBH4.EXPT.NP04.001:COUNTS"      ].fVar;
      XBH4EXPTNP04002COUNTS       = mapBasicDoubleVars["XBH4.EXPT.NP04.002:COUNTS"      ].fVar;
      XBH4EXPTNP04003COUNTS       = mapBasicDoubleVars["XBH4.EXPT.NP04.003:COUNTS"      ].fVar;
      XBH4EXPTNP04004COUNTS       = mapBasicDoubleVars["XBH4.EXPT.NP04.004:COUNTS"      ].fVar;
      XBH4EXPTNP04009COUNTS       = mapBasicDoubleVars["XBH4.EXPT.NP04.009:COUNTS"      ].fVar;
      XBH4EXPTNP04010COUNTS       = mapBasicDoubleVars["XBH4.EXPT.NP04.010:COUNTS"      ].fVar;
      XBH4EXPTNP04011COUNTS       = mapBasicDoubleVars["XBH4.EXPT.NP04.011:COUNTS"      ].fVar;
      XBH4EXPTNP04012COUNTS       = mapBasicDoubleVars["XBH4.EXPT.NP04.012:COUNTS"      ].fVar;
      XBH4XCET022713COUNTS_TRIG   = mapBasicDoubleVars["XBH4.XCET.022.713:COUNTS_TRIG"  ].fVar;
      XBH4XCET022713PRESSURE      = mapBasicDoubleVars["XBH4.XCET.022.713:PRESSURE"     ].fVar;
      XBH4XCET022713SIMPLE_COUNTS = mapBasicDoubleVars["XBH4.XCET.022.713:SIMPLE_COUNTS"].fVar;
      XBH4XCET022716COUNTS_TRIG   = mapBasicDoubleVars["XBH4.XCET.022.716:COUNTS_TRIG"  ].fVar;
      XBH4XCET022716PRESSURE      = mapBasicDoubleVars["XBH4.XCET.022.716:PRESSURE"     ].fVar;
      XBH4XCET022716SIMPLE_COUNTS = mapBasicDoubleVars["XBH4.XCET.022.716:SIMPLE_COUNTS"].fVar;
      XBH4XCSH022694POS_JAW1_MEAS = mapBasicDoubleVars["XBH4.XCSH.022.694:POS_JAW1_MEAS"].fVar;
      XBH4XCSH022694POS_JAW2_MEAS = mapBasicDoubleVars["XBH4.XCSH.022.694:POS_JAW2_MEAS"].fVar;
      XBH4XSCI022680COUNTS        = mapBasicDoubleVars["XBH4.XSCI.022.680:COUNTS"       ].fVar; 
      SPST2INTENSITY              = mapBasicDoubleVars["SPS.T2:INTENSITY"               ].fVar;

      XBH4GENERALTRIGGER_NCounts = (mapTDCAcquisitions.count( "XBH4.GENERALTRIGGER") ? mapTDCAcquisitions["XBH4.GENERALTRIGGER"].getTimestampCount() : 0);  
      XBH4XTDC022713_NCounts     = (mapTDCAcquisitions.count( "XBH4.XTDC.022.713"  ) ? mapTDCAcquisitions["XBH4.XTDC.022.713"]  .getTimestampCount() : 0); 
      XBH4XTDC022716_NCounts     = (mapTDCAcquisitions.count( "XBH4.XTDC.022.716"  ) ? mapTDCAcquisitions["XBH4.XTDC.022.716"]  .getTimestampCount() : 0); 
      XBTF022687A_NCounts        = (mapTDCAcquisitions.count( "XBTF.022.687A"      ) ? mapTDCAcquisitions["XBTF.022.687A"]      .getTimestampCount() : 0);  
      XBTF022687B_NCounts        = (mapTDCAcquisitions.count( "XBTF.022.687B"      ) ? mapTDCAcquisitions["XBTF.022.687B"]      .getTimestampCount() : 0);  
      XBTF022716A_NCounts        = (mapTDCAcquisitions.count( "XBTF.022.716A"      ) ? mapTDCAcquisitions["XBTF.022.716A"]      .getTimestampCount() : 0);  
      XBTF022716B_NCounts        = (mapTDCAcquisitions.count( "XBTF.022.716B"      ) ? mapTDCAcquisitions["XBTF.022.716B"]      .getTimestampCount() : 0); 

      XBH4GENERALTRIGGER_TimestampNS = GenTrigRecords[i].fNSAccurate; XBH4GENERALTRIGGER_TimestampFracAccuracy = GenTrigRecords[i].fFracPieceHR;         
      if(mapGenTrigEventsToTDCEvents.count({i,"XBH4.XTDC.022.713"}))
      {
        for(unsigned int j = 0; j < mapGenTrigEventsToTDCEvents[{i,"XBH4.XTDC.022.713"}].size(); j++){
          XBH4XTDC022713_TimestampNS          .push_back(XBH4XTDC022713_Records[mapGenTrigEventsToTDCEvents[{i,"XBH4.XTDC.022.713"}][j]].fNSAccurate );
          XBH4XTDC022713_TimestampFracAccuracy.push_back(XBH4XTDC022713_Records[mapGenTrigEventsToTDCEvents[{i,"XBH4.XTDC.022.713"}][j]].fFracPieceHR);
        }
      }
      if(mapGenTrigEventsToTDCEvents.count({i,"XBH4.XTDC.022.716"}))
      {
        for(unsigned int j = 0; j < mapGenTrigEventsToTDCEvents[{i,"XBH4.XTDC.022.716"}].size(); j++){
          XBH4XTDC022716_TimestampNS          .push_back(XBH4XTDC022716_Records[mapGenTrigEventsToTDCEvents[{i,"XBH4.XTDC.022.716"}][j]].fNSAccurate );
          XBH4XTDC022716_TimestampFracAccuracy.push_back(XBH4XTDC022716_Records[mapGenTrigEventsToTDCEvents[{i,"XBH4.XTDC.022.716"}][j]].fFracPieceHR);
        }
      }

      if(XBH4XTDC022713_TimestampNS.size()>0){XBH4XTDC022713_On = true;}
      if(XBH4XTDC022716_TimestampNS.size()>0){XBH4XTDC022716_On = true;}

      if(mapGenTrigEventsToTDCEvents.count({i,"XBTF.022.687A"}))
      {
        for(unsigned int j = 0; j < mapGenTrigEventsToTDCEvents[{i,"XBTF.022.687A"}].size(); j++){
          XBTF022687A_TimestampNS          .push_back(XBTF022687A_Records[mapGenTrigEventsToTDCEvents[{i,"XBTF.022.687A"}][j]].fNSAccurate );
          XBTF022687A_TimestampFracAccuracy.push_back(XBTF022687A_Records[mapGenTrigEventsToTDCEvents[{i,"XBTF.022.687A"}][j]].fFracPieceHR);
        }
      }
      if(mapGenTrigEventsToTDCEvents.count({i,"XBTF.022.687B"}))
      {
        for(unsigned int j = 0; j < mapGenTrigEventsToTDCEvents[{i,"XBTF.022.687B"}].size(); j++){
          XBTF022687B_TimestampNS          .push_back(XBTF022687B_Records[mapGenTrigEventsToTDCEvents[{i,"XBTF.022.687B"}][j]].fNSAccurate );
          XBTF022687B_TimestampFracAccuracy.push_back(XBTF022687B_Records[mapGenTrigEventsToTDCEvents[{i,"XBTF.022.687B"}][j]].fFracPieceHR);
        }
      }
      if(mapGenTrigEventsToTDCEvents.count({i,"XBTF.022.716A"}))
      {
        for(unsigned int j = 0; j < mapGenTrigEventsToTDCEvents[{i,"XBTF.022.716A"}].size(); j++){
          XBTF022716A_TimestampNS          .push_back(XBTF022716A_Records[mapGenTrigEventsToTDCEvents[{i,"XBTF.022.716A"}][j]].fNSAccurate );
          XBTF022716A_TimestampFracAccuracy.push_back(XBTF022716A_Records[mapGenTrigEventsToTDCEvents[{i,"XBTF.022.716A"}][j]].fFracPieceHR);
        }
      }
      if(mapGenTrigEventsToTDCEvents.count({i,"XBTF.022.716B"}))
      {
        for(unsigned int j = 0; j < mapGenTrigEventsToTDCEvents[{i,"XBTF.022.716B"}].size(); j++){
          XBTF022716B_TimestampNS          .push_back(XBTF022716B_Records[mapGenTrigEventsToTDCEvents[{i,"XBTF.022.716B"}][j]].fNSAccurate );
          XBTF022716B_TimestampFracAccuracy.push_back(XBTF022716B_Records[mapGenTrigEventsToTDCEvents[{i,"XBTF.022.716B"}][j]].fFracPieceHR);
        }
      }

      for(unsigned int j = 0; j < mapGenTrigEventToTOFs[i].size(); j++)
      {
        tofChannel.push_back(mapGenTrigEventToTOFs[i][j].fChannel);
        tof.push_back(mapGenTrigEventToTOFs[i][j].fTOF);
      }

      std::vector<XBPFCoincidence> vecXBPFCoincidence = mapVecXBPFCoincidence[i];
      for(unsigned int j = 0; j < vecXBPFCoincidence.size(); j++)
      {
        XBH4XBPF022697_TimestampNS   .push_back((unsigned long long)vecXBPFCoincidence[j].f697.fTriggerTimestamp); 
        XBH4XBPF022697_NFibresHit    .push_back((unsigned int)vecXBPFCoincidence[j].f697.fFibresHit); 
        XBH4XBPF022697_Span          .push_back((unsigned int)vecXBPFCoincidence[j].f697.fSpan);             
        XBH4XBPF022697_HitCoordinates.push_back(vecXBPFCoincidence[j].f697.fPhysicalCoordinates); 
        XBH4XBPF022698_TimestampNS   .push_back((unsigned long long)vecXBPFCoincidence[j].f698.fTriggerTimestamp); 
        XBH4XBPF022698_NFibresHit    .push_back((unsigned int)vecXBPFCoincidence[j].f698.fFibresHit); 
        XBH4XBPF022698_Span          .push_back((unsigned int)vecXBPFCoincidence[j].f698.fSpan);             
        XBH4XBPF022698_HitCoordinates.push_back(vecXBPFCoincidence[j].f698.fPhysicalCoordinates); 
        XBH4XBPF022701_TimestampNS   .push_back((unsigned long long)vecXBPFCoincidence[j].f701.fTriggerTimestamp); 
        XBH4XBPF022701_NFibresHit    .push_back((unsigned int)vecXBPFCoincidence[j].f701.fFibresHit); 
        XBH4XBPF022701_Span          .push_back((unsigned int)vecXBPFCoincidence[j].f701.fSpan);             
        XBH4XBPF022701_HitCoordinates.push_back(vecXBPFCoincidence[j].f701.fPhysicalCoordinates); 
        XBH4XBPF022702_TimestampNS   .push_back((unsigned long long)vecXBPFCoincidence[j].f702.fTriggerTimestamp); 
        XBH4XBPF022702_NFibresHit    .push_back((unsigned int)vecXBPFCoincidence[j].f702.fFibresHit); 
        XBH4XBPF022702_Span          .push_back((unsigned int)vecXBPFCoincidence[j].f702.fSpan);             
        XBH4XBPF022702_HitCoordinates.push_back(vecXBPFCoincidence[j].f702.fPhysicalCoordinates); 
        XBH4XBPF022707_TimestampNS   .push_back((unsigned long long)vecXBPFCoincidence[j].f707.fTriggerTimestamp); 
        XBH4XBPF022707_NFibresHit    .push_back((unsigned int)vecXBPFCoincidence[j].f707.fFibresHit); 
        XBH4XBPF022707_Span          .push_back((unsigned int)vecXBPFCoincidence[j].f707.fSpan);             
        XBH4XBPF022707_HitCoordinates.push_back(vecXBPFCoincidence[j].f707.fPhysicalCoordinates); 
        XBH4XBPF022708_TimestampNS   .push_back((unsigned long long)vecXBPFCoincidence[j].f708.fTriggerTimestamp); 
        XBH4XBPF022708_NFibresHit    .push_back((unsigned int)vecXBPFCoincidence[j].f708.fFibresHit); 
        XBH4XBPF022708_Span          .push_back((unsigned int)vecXBPFCoincidence[j].f708.fSpan);             
        XBH4XBPF022708_HitCoordinates.push_back(vecXBPFCoincidence[j].f708.fPhysicalCoordinates); 
        XBH4XBPF022716_TimestampNS   .push_back((unsigned long long)vecXBPFCoincidence[j].f716.fTriggerTimestamp); 
        XBH4XBPF022716_NFibresHit    .push_back((unsigned int)vecXBPFCoincidence[j].f716.fFibresHit); 
        XBH4XBPF022716_Span          .push_back((unsigned int)vecXBPFCoincidence[j].f716.fSpan);             
        XBH4XBPF022716_HitCoordinates.push_back(vecXBPFCoincidence[j].f716.fPhysicalCoordinates); 
        XBH4XBPF022717_TimestampNS   .push_back((unsigned long long)vecXBPFCoincidence[j].f717.fTriggerTimestamp); 
        XBH4XBPF022717_NFibresHit    .push_back((unsigned int)vecXBPFCoincidence[j].f717.fFibresHit); 
        XBH4XBPF022717_Span          .push_back((unsigned int)vecXBPFCoincidence[j].f717.fSpan);             
        XBH4XBPF022717_HitCoordinates.push_back(vecXBPFCoincidence[j].f717.fPhysicalCoordinates); 

        nMomenta             .push_back(vecXBPFCoincidence[j].fMomentum.size());
        deflectionAngle      .push_back(vecXBPFCoincidence[j].fDeflectionAngle);
        reconstructedMomentum.push_back(vecXBPFCoincidence[j].fMomentum);
      }

      eventRank = rankEvent(XBH4XBPF022697_NFibresHit, XBH4XBPF022701_NFibresHit, XBH4XBPF022702_NFibresHit, 
          XBH4XBPF022697_Span,       XBH4XBPF022701_Span,       XBH4XBPF022702_Span,      
          mapGenTrigEventToTOFs[i].size(), XBH4XTDC022713_TimestampNS.size(), XBH4XTDC022716_TimestampNS.size());

      //DO PID FOR NON-DEGENERATE TIME OF FLIGHTS.
      if(mapGenTrigEventToTOFs[i].size()==1 && XBH4XTDC022713_TimestampNS.size()<2 && XBH4XTDC022716_TimestampNS.size()<2)
      {
        validTOF     = tof.at(0); 
        validChannel = tofChannel.at(0);
        PID          = getPID(validTOF, XBH4XTDC022713_On, XBH4XTDC022716_On, (XBH4XCET022716PRESSURE >= XBH4XCET022713PRESSURE ? true : false), 
                              referenceMomentum);
      }
      //TO SELECT A VALID TOF AND DO A PID ANALYSIS, NEED TO HAVE NON-DEGENERATE XCET SIGNALS.
      else if(XBH4XTDC022713_TimestampNS.size()<2 && XBH4XTDC022716_TimestampNS.size()<2)
      {
        assignValidTOFAndPID(validTOF, validChannel, PID, tof, tofChannel, XBH4XTDC022713_TimestampNS, XBH4XTDC022716_TimestampNS, referenceMomentum, 
                            (XBH4XCET022716PRESSURE >= XBH4XCET022713PRESSURE ? true : false));
      }
      else
      {
        validTOF     = -99.;
        validChannel = "nouniquietof";
        PID          = "nopid";
      }

      t_Out->Fill();
    }
  }

  t_Out->Write();
  delete f_Out;

  return;
}
