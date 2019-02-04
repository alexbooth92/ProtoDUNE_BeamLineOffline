#include "plotting.h"

void plot_dTAcq(std::map<std::string,Detector> &map_Detector, std::map<std::string,std::vector<TH1I*>> &map_NameToVecHist_I)
{
  for(auto it_Detector : map_Detector)
  {
    if(it_Detector.second.getType()=="XBPF")
    {
      TH1I *h_dTAcq = new TH1I("h_dTAcq_"+(TString)it_Detector.second.getDetName(),"h_dTAcq_"+(TString)it_Detector.second.getDetName(), 71, -0.5, 70.5);

      SetHistTitles(h_dTAcq, (TString)it_Detector.second.getDetName(), "Time Between DIP Data Acquisitions, (s)", "Acquisition Pairs");

      std::vector<AcquisitionXBPF> vec_DetAcq;
      it_Detector.second.getAcquisitions(vec_DetAcq);
      if(vec_DetAcq.size()!=0)
      {
        for(unsigned int j = 0; j < vec_DetAcq.size()-1; j++)
        {
          double dTAcq = vec_DetAcq.at(j+1).getTimestamp() - vec_DetAcq.at(j).getTimestamp();
          h_dTAcq->Fill(dTAcq);
        }
      }
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_dTAcq);
    }
    else if(it_Detector.second.getType()=="XBTF")
    {
      TH1I *h_dTAcq = new TH1I("h_dTAcq_"+(TString)it_Detector.second.getDetName(),"h_dTAcq_"+(TString)it_Detector.second.getDetName(), 71, -0.5, 70.5);

      SetHistTitles(h_dTAcq, (TString)it_Detector.second.getDetName(), "Time Between DIP Data Acquisitions, (s)", "Acquisition Pairs");

      std::vector<AcquisitionXBTF> vec_DetAcq;
      it_Detector.second.getAcquisitions(vec_DetAcq);
      if(vec_DetAcq.size()!=0)
      {
        for(unsigned int j = 0; j < vec_DetAcq.size()-1; j++)
        {
          double dTAcq = vec_DetAcq.at(j+1).getTimestamp() - vec_DetAcq.at(j).getTimestamp();
          h_dTAcq->Fill(dTAcq);
        }
      }
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_dTAcq);
    }
  }

  return;
}

void plot_AcqLevel(std::map<std::string,Detector> &map_Detector, std::map<std::string,std::vector<TH1I*>> &map_NameToVecHist_I)
{
  for(auto it_Detector : map_Detector)
  {
    if(it_Detector.second.getType()=="XBPF")
    {
      TH1I *h_NEventsPerAcq = new TH1I("h_NEventsPerAcq_"+(TString)it_Detector.second.getDetName(),"h_NEventsPerAcq_"+(TString)it_Detector.second.getDetName(), 3001, -0.5, 3000.5);
      TH1I *h_FracBuffer    = new TH1I("h_FracBuffer_"   +(TString)it_Detector.second.getDetName(),"h_FracBuffer_"   +(TString)it_Detector.second.getDetName(), 50,  -0.05,   1.05);

      SetHistTitles(h_NEventsPerAcq, (TString)it_Detector.second.getDetName(), "Number of Events / Acqusition",         "Acquisitions");
      SetHistTitles(h_FracBuffer,    (TString)it_Detector.second.getDetName(), "Fraction of EventsData Array Filled",   "Acquisitions");

      std::vector<AcquisitionXBPF> vec_DetAcq;
      it_Detector.second.getAcquisitions(vec_DetAcq);
      for(unsigned int j = 0; j < vec_DetAcq.size(); j++)
      {
        h_NEventsPerAcq->Fill(vec_DetAcq.at(j).getNEventRec());
        h_FracBuffer->Fill(vec_DetAcq.at(j).getNEventRec()/(double)3000);
      }
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_NEventsPerAcq);
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_FracBuffer);
    }
    else if(it_Detector.second.getType()=="XBTF")
    {
      TH1I *h_NEventsPerAcq = new TH1I("h_NEventsPerAcq_"+(TString)it_Detector.second.getDetName(),"h_NEventsPerAcq_"+(TString)it_Detector.second.getDetName(), 30001, -0.5, 30000.5);
      TH1I *h_FracBuffer    = new TH1I("h_FracBuffer_"   +(TString)it_Detector.second.getDetName(),"h_FracBuffer_"   +(TString)it_Detector.second.getDetName(), 50,    -0.05,   1.05);

      SetHistTitles(h_NEventsPerAcq, (TString)it_Detector.second.getDetName(), "Number of Events / Acqusition",         "Acquisitions");
      SetHistTitles(h_FracBuffer,    (TString)it_Detector.second.getDetName(), "Fraction of EventsData Array Filled",   "Acquisitions");

      std::vector<AcquisitionXBTF> vec_DetAcq;
      it_Detector.second.getAcquisitions(vec_DetAcq);
      for(unsigned int j = 0; j < vec_DetAcq.size(); j++)
      {
        h_NEventsPerAcq->Fill(vec_DetAcq.at(j).getTimestampCount());
        h_FracBuffer->Fill(vec_DetAcq.at(j).getTimestampCount()/(double)30000);
      }
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_NEventsPerAcq);
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_FracBuffer);
    }
  }

  return;
}

void plot_EventLevel(std::map<std::string,Detector> &map_Detector, std::map<std::string,std::vector<TH1I*>> &map_NameToVecHist_I, std::map<std::string,std::vector<TH1D*>> &map_NameToVecHist_D)
{
  for(auto it_Detector : map_Detector)
  {
    if(it_Detector.second.getType()=="XBPF")
    {
      double tLow        = it_Detector.second.getFirstTriggerTimestamp();   
      double tHigh       = it_Detector.second.getLastTriggerTimestamp();   
      double timeElapsed = tHigh - tLow;
      double sensitivity = 1.e-3;
      long   nBins       = std::ceil(timeElapsed/sensitivity);
      double upperbin    = tLow + (double)nBins*(double)sensitivity;
      if(tLow==2.e9)
      {
        nBins = 1;
        tLow  = it_Detector.second.getFirstTimestamp();
        upperbin = tLow + sensitivity;
      }

      TH1I *h_NFibresPerEvent        = new TH1I("h_NFibresPerEvent_"       +(TString)it_Detector.second.getDetName(),"h_NFibresPerEvent_"       
                                                                           +(TString)it_Detector.second.getDetName(),  11, -0.5,  10.5);
      TH1I *h_Span                   = new TH1I("h_Span_"                  +(TString)it_Detector.second.getDetName(),"h_Span_"                  
                                                                           +(TString)it_Detector.second.getDetName(),   7, -0.5,   6.5);
      TH1I *h_TriggerTS              = new TH1I("h_TriggerTS_"             +(TString)it_Detector.second.getDetName(),"h_TriggerTS_"             
                                                                           +(TString)it_Detector.second.getDetName(),nBins, tLow, upperbin);
      TH1I *h_EventTS                = new TH1I("h_EventTS_"               +(TString)it_Detector.second.getDetName(),"h_EventTS_"               
                                                                           +(TString)it_Detector.second.getDetName(),nBins, tLow, upperbin); 
      TH1I *h_EventTimeSinceFirst    = new TH1I("h_EventTimeSinceFirst_"   +(TString)it_Detector.second.getDetName(),"h_EventTimeSinceFirst_"   
                                                                           +(TString)it_Detector.second.getDetName(),100,    0,   5e9); 
      TH1I *h_TriggerTimeSinceFirst  = new TH1I("h_TriggerTimeSinceFirst_" +(TString)it_Detector.second.getDetName(),"h_TriggerTimeSinceFirst_" 
                                                                           +(TString)it_Detector.second.getDetName(),100,    0,   5e9); 
      TH1D *h_NHitsPerFibre          = new TH1D("h_NHitsPerFibre_"         +(TString)it_Detector.second.getDetName(),"h_NHitsPerFibre_"
                                                                           +(TString)it_Detector.second.getDetName(), 192, -0.5, 191.5);

      SetHistTitles(h_NFibresPerEvent      , (TString)it_Detector.second.getDetName(), "Number of Fibres Hit / Event",                   "Events"  ); 
      SetHistTitles(h_Span                 , (TString)it_Detector.second.getDetName(), "Hit Fibre Span",                                 "Events"  ); 
      SetHistTitles(h_TriggerTS            , (TString)it_Detector.second.getDetName(), "Trigger Timestamp, Unix",                        "Triggers"); 
      SetHistTitles(h_EventTS              , (TString)it_Detector.second.getDetName(), "Event Timestamp, Unix",                          "Events"  ); 
      SetHistTitles(h_EventTimeSinceFirst  , (TString)it_Detector.second.getDetName(), "Event Time Elapsed Since First Event, (ns)",     "Events"  ); 
      SetHistTitles(h_TriggerTimeSinceFirst, (TString)it_Detector.second.getDetName(), "Trigger Time Elapsed Since First Trigger, (ns)", "Triggers"); 
      SetHistTitles(h_NHitsPerFibre        , (TString)it_Detector.second.getDetName(), "Number of Hits / Fibre,(Fibre No.)",             "Hits"    ); 

      std::vector<AcquisitionXBPF> vec_DetAcq;
      it_Detector.second.getAcquisitions(vec_DetAcq);
      unsigned int eventCount = 0;
      for(unsigned int j = 0; j < vec_DetAcq.size(); j++)
      {
        bool b_IsFirstRecord = true;
        AcquisitionXBPF::EventRecordHR firstRecord;
        unsigned int nEvents = vec_DetAcq.at(j).getDataHR().size();
        unsigned int nNZeroEvents = vec_DetAcq.at(j).getNNonZeroEvents();
        for(unsigned int k = 0; k < nNZeroEvents; k++)
        //for(unsigned int k = 0; k < nEvents; k++)
        {
          if(b_IsFirstRecord)
          {
            firstRecord = vec_DetAcq.at(j).getDataHR().at(0);
            b_IsFirstRecord = false;
          }
          AcquisitionXBPF::EventRecordHR currentRecord = vec_DetAcq.at(j).getDataHR().at(k);

          h_NFibresPerEvent->Fill(currentRecord.fNFibresHit < 10 ? currentRecord.fNFibresHit : 10);
          h_Span     ->Fill((currentRecord.fSpan < 6) ? currentRecord.fSpan : 6);
          h_TriggerTS->Fill(currentRecord.fTriggerTimestamp/1.e9);
          h_EventTS  ->Fill(currentRecord.fEventTimestamp/1.e9);
          h_TriggerTimeSinceFirst->Fill(currentRecord.fTriggerTimestamp - firstRecord.fTriggerTimestamp);
          h_EventTimeSinceFirst  ->Fill(currentRecord.fEventTimestamp   - firstRecord.fEventTimestamp);
          for(unsigned int l = 0; l < currentRecord.fFibresList.size(); l++)
          {
            h_NHitsPerFibre->Fill(currentRecord.fFibresList[l]);
          }
        }
        eventCount+=nEvents;
      }
      /*
      if(eventCount!=0)
      {
        h_NHitsPerFibre->Scale(1./(double)eventCount);
      }*/

      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_NFibresPerEvent);
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_Span);
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_TriggerTS);
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_EventTS);
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_TriggerTimeSinceFirst);
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_EventTimeSinceFirst);
      map_NameToVecHist_D[it_Detector.second.getDetName()].push_back(h_NHitsPerFibre);
    }
    else if(it_Detector.second.getType()=="XBTF")
    {
      double tLow  = it_Detector.second.getFirstTriggerTimestamp();   
      double tHigh = it_Detector.second.getLastTriggerTimestamp();   
      double timeElapsed = tHigh - tLow;
      double sensitivity = 1.e-3;
      int    nBins       = std::ceil(timeElapsed/sensitivity);
      double upperbin    = tLow + (double)nBins*(double)sensitivity;
      if(tLow==2.e9)
      {
        nBins = 1;
        tLow  = it_Detector.second.getFirstTimestamp();
        upperbin = tLow + sensitivity;
      }

      TH1I *h_EventTS_SecondsPart     = new TH1I("h_EventTS_SecondsPart_"    +(TString)it_Detector.second.getDetName(), "h_EventTS_SecondsPart_"             
                                                                             +(TString)it_Detector.second.getDetName(), 100, tLow, upperbin);
      TH1I *h_EventTS_NanosecondsPart = new TH1I("h_EventTS_NanosecondsPart_"+(TString)it_Detector.second.getDetName(), "h_EventTS_NanosecondsPart_"             
                                                                             +(TString)it_Detector.second.getDetName(), 100, 0, 1e9);
      TH1I *h_EventTS                 = new TH1I("h_EventTS_"                +(TString)it_Detector.second.getDetName(), "h_EventTS_"               
                                                                             +(TString)it_Detector.second.getDetName(), 100, tLow, upperbin); 
      TH1I *h_EventTimeSinceFirst     = new TH1I("h_EventTimeSinceFirst_"    +(TString)it_Detector.second.getDetName(), "h_EventTimeSinceFirst_"   
                                                                             +(TString)it_Detector.second.getDetName(), 100, 0, 5e9); 

      SetHistTitles(h_EventTS_SecondsPart,     (TString)it_Detector.second.getDetName(), "Event Timestamp Seconds Part, Unix",         "Events");
      SetHistTitles(h_EventTS_NanosecondsPart, (TString)it_Detector.second.getDetName(), "Event Timestamp Nanoseconds Part, (ns)",     "Events");
      SetHistTitles(h_EventTS,                 (TString)it_Detector.second.getDetName(), "Event Timestamp, (s)",                       "Events");
      SetHistTitles(h_EventTimeSinceFirst,     (TString)it_Detector.second.getDetName(), "Event Time Elapsed Since First Event, (ns)", "Events");

      std::vector<AcquisitionXBTF> vec_DetAcq;
      it_Detector.second.getAcquisitions(vec_DetAcq);
      for(unsigned int j = 0; j < vec_DetAcq.size(); j++)
      {
        bool b_IsFirstRecord = true;
        AcquisitionXBTF::EventRecordHR firstRecord;
        for(unsigned int k = 0; k < vec_DetAcq.at(j).getDataHR().size(); k++)
        {
          if(b_IsFirstRecord)
          {
            firstRecord = vec_DetAcq.at(j).getDataHR().at(0);
            b_IsFirstRecord = false;
          }
          AcquisitionXBTF::EventRecordHR currentRecord = vec_DetAcq.at(j).getDataHR().at(k);

          h_EventTS_SecondsPart    ->Fill(currentRecord.fSeconds);
          h_EventTS_NanosecondsPart->Fill(currentRecord.fSubSeconds);
          h_EventTS                ->Fill(currentRecord.fSeconds_FullTime);
          h_EventTimeSinceFirst    ->Fill((((double)currentRecord.fSeconds-(double)firstRecord.fSeconds)+(currentRecord.fSubSeconds - firstRecord.fSubSeconds)/1.e9)*1.e9);
        }
      }

      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_EventTS_SecondsPart);
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_EventTS_NanosecondsPart);
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_EventTS);
      map_NameToVecHist_I[it_Detector.second.getDetName()].push_back(h_EventTimeSinceFirst);
    }
  }

  return;
}

void plotPROF(std::vector<PROFCoincidenceRecord::PROFCoincidence> &vec_PROFCO, std::vector<TH1D*> &vec_PROFHist1D, std::vector<TH2D*> &vec_PROFHist2D, TString const &sType) 
{
  TH1D *h_CosTheta   = new TH1D("h_CosTheta_"  +sType, "h_CosTheta_"  +sType, 100, -1.1,  1.1);
  TH1D *h_Momentum   = new TH1D("h_Momentum_"  +sType, "h_Momentum_"  +sType, 100,    0,   12);
  TH1D *h_Theta      = new TH1D("h_Theta_"     +sType, "h_Theta_"     +sType, 100,    0,   12);
  TH1D *h_FibComb    = new TH1D("h_FibComb_"   +sType, "h_FibComb_"   +sType,  21, -0.5, 20.5);
  TH1D *h_MultiPROF1 = new TH1D("h_MultiPROF1_"+sType, "h_MultiPROF1_"+sType,   7, -0.5,  6.5);
  TH1D *h_MultiPROF2 = new TH1D("h_MultiPROF2_"+sType, "h_MultiPROF2_"+sType,   7, -0.5,  6.5);
  TH1D *h_MultiPROF3 = new TH1D("h_MultiPROF3_"+sType, "h_MultiPROF3_"+sType,   7, -0.5,  6.5);

  if(sType=="UNIQ")
  {
    SetHistTitles(h_CosTheta  , "Cosine of Deflection Angle, Unique Events",            "CosTheta",          "Events");
    SetHistTitles(h_Momentum  , "Momentum, Unique Events",                              "Momentum, (GeV)",   "Events");
    SetHistTitles(h_Theta     , "Deflection Angle, Unique Events",                      "Deflection, (Deg)", "Events");
    SetHistTitles(h_FibComb   , "Number of Possible Fibre Combinations / Unique Event", "Combinations",      "Events");
    SetHistTitles(h_MultiPROF1, "Event Coincidence Degeneracy, XBPF022697", "Combinations","Events");
    SetHistTitles(h_MultiPROF2, "Event Coincidence Degeneracy, XBPF022701", "Combinations","Events");
    SetHistTitles(h_MultiPROF3, "Event Coincidence Degeneracy, XBPF022702", "Combinations","Events");
  }
  else
  {
    SetHistTitles(h_CosTheta  , "Cosine of Deflection Angle, Degenerate Events",            "CosTheta",          "Events");
    SetHistTitles(h_Momentum  , "Momentum, Degenerate Events",                              "Momentum, (GeV)",   "Events");
    SetHistTitles(h_Theta     , "Deflection Angle, Degenerate Events",                      "Deflection, (Deg)", "Events");
    SetHistTitles(h_FibComb   , "Number of Possible Fibre Combinations / Degenerate Event", "Combinations",      "Events");
    SetHistTitles(h_MultiPROF1, "Event Coincidence Degeneracy, XBPF022697", "Combinations","Events");
    SetHistTitles(h_MultiPROF2, "Event Coincidence Degeneracy, XBPF022701", "Combinations","Events");
    SetHistTitles(h_MultiPROF3, "Event Coincidence Degeneracy, XBPF022702", "Combinations","Events");
  }

  for(unsigned int i = 0; i < vec_PROFCO.size(); i++)
  {
    for(unsigned int j = 0; j < vec_PROFCO[i].getMomentum().size(); j++)
    {
      if(sType=="UNIQ")
      {
        if(vec_PROFCO[i].getMomentum().size()==1)
        {
          h_CosTheta->Fill(vec_PROFCO[i].getCosTheta()[j]);  
          h_Momentum->Fill(vec_PROFCO[i].getMomentum()[j]);  
          h_Theta   ->Fill(vec_PROFCO[i].getTheta()[j]*180./3.141593);
          h_FibComb   ->Fill((vec_PROFCO[i].getNFibreCombinations()<=19 ? vec_PROFCO[i].getNFibreCombinations() : 20));
          h_MultiPROF1->Fill((vec_PROFCO[i].fMultiplicity[0]<=5 ? vec_PROFCO[i].fMultiplicity[0] : 6));
          h_MultiPROF2->Fill((vec_PROFCO[i].fMultiplicity[1]<=5 ? vec_PROFCO[i].fMultiplicity[1] : 6));
          h_MultiPROF3->Fill((vec_PROFCO[i].fMultiplicity[2]<=5 ? vec_PROFCO[i].fMultiplicity[2] : 6));
        }
      }
      else
      {
        h_CosTheta->Fill(vec_PROFCO[i].getCosTheta()[j]);  
        h_Momentum->Fill(vec_PROFCO[i].getMomentum()[j]);  
        h_Theta   ->Fill(vec_PROFCO[i].getTheta()[j]*180./3.141593);
        h_FibComb   ->Fill((vec_PROFCO[i].getNFibreCombinations()<=19 ? vec_PROFCO[i].getNFibreCombinations() : 20));
        h_MultiPROF1->Fill((vec_PROFCO[i].fMultiplicity[0]<=5 ? vec_PROFCO[i].fMultiplicity[0] : 6));
        h_MultiPROF2->Fill((vec_PROFCO[i].fMultiplicity[1]<=5 ? vec_PROFCO[i].fMultiplicity[1] : 6));
        h_MultiPROF3->Fill((vec_PROFCO[i].fMultiplicity[2]<=5 ? vec_PROFCO[i].fMultiplicity[2] : 6));
      }
    }
  }

  vec_PROFHist1D.push_back(h_CosTheta  );
  vec_PROFHist1D.push_back(h_Momentum  );
  vec_PROFHist1D.push_back(h_Theta     );
  vec_PROFHist1D.push_back(h_FibComb   );
  vec_PROFHist1D.push_back(h_MultiPROF1);
  vec_PROFHist1D.push_back(h_MultiPROF2);
  vec_PROFHist1D.push_back(h_MultiPROF3);

  return;
}

void plotTF(std::vector<TFCoincidenceRecord::TFCoincidence> &vec_TFCO, std::vector<TH1D*> &vec_TFHist1D, std::vector<TH2D*> &vec_TFHist2D, TString const &sType) 
{
  TH1D *h_TF      = new TH1D("h_TF_"     +sType,"h_TF_"      +sType, 200,    0,  500);
  TH1D *h_TF_ROI  = new TH1D("h_TF_ROI_" +sType,"h_TF_ROI_"  +sType, 134,  140,  140+134*0.3);
  TH1D *h_MultiUS = new TH1D("h_MultiUS_"+sType, "h_MultiUS_"+sType,   7, -0.5,  6.5); 
  TH1D *h_MultiDS = new TH1D("h_MultiDS_"+sType, "h_MultiDS_"+sType,   7, -0.5,  6.5); 

  if(sType=="UNIQ")
  {
    SetHistTitles(h_TF,      "Time of Flight from XBTFs, Unique Events", "Time, (ns)",   "Events");
    SetHistTitles(h_TF_ROI,  "Time of Flight from XBTFs, Unique Events", "Time, (ns)",   "Events");
    SetHistTitles(h_MultiUS, "Event Coincidence Degeneracy, US",         "Combinations", "Events");
    SetHistTitles(h_MultiDS, "Event Coincidence Degeneracy, DS",         "Combinations", "Events");
  }
  else
  {
    SetHistTitles(h_TF,      "Time of Flight from XBTFs, Degenerate Events", "Time, (ns)",   "Events");
    SetHistTitles(h_TF_ROI,  "Time of Flight from XBTFs, Degenerate Events", "Time, (ns)",   "Events");
    SetHistTitles(h_MultiUS, "Event Coincidence Degeneracy, US",             "Combinations", "Events");
    SetHistTitles(h_MultiDS, "Event Coincidence Degeneracy, DS",             "Combinations", "Events");
  }

  for(unsigned int i = 0; i < vec_TFCO.size(); i++)
  {
    h_TF->Fill(vec_TFCO[i].fTF);
    h_TF_ROI->Fill(vec_TFCO[i].fTF);
    h_MultiUS->Fill((vec_TFCO[i].fMultiplicity.first <=5 ? vec_TFCO[i].fMultiplicity.first  : 6));
    h_MultiDS->Fill((vec_TFCO[i].fMultiplicity.second<=5 ? vec_TFCO[i].fMultiplicity.second : 6));
  }

  vec_TFHist1D.push_back(h_TF);
  vec_TFHist1D.push_back(h_TF_ROI);
  vec_TFHist1D.push_back(h_MultiUS);
  vec_TFHist1D.push_back(h_MultiDS);

  return;
}

void plotTFPROF(std::vector<CombinedCoincidenceRecord::TFPROFCoincidence> &vec_TFPROFCO, std::vector<TH1D*> &vec_TFPROFHist1D, std::vector<TH2D*> &vec_TFPROFHist2D, TString const &sType) 
{
  TH1D *h_Multi = new TH1D("h_Multi_"+sType, "h_Multi_"+sType,   7, -0.5,  6.5); 
  TH1D *h_Mass  = new TH1D("h_Mass_"+sType,  "h_Mass_" +sType, 100,    0, 1000);
  TH2D *h_TFMom = new TH2D("h_TFMom_"+sType, "h_TFMom_"+sType, 100,    0, 12, 134, 140, 140+134*0.3);

  if(sType=="UNIQ")
  {
    SetHistTitles(h_Mass,  "Mass, Unique Events",                        "Mass, (GeV)",     "Time, (ns)");
    SetHistTitles(h_TFMom, "Time of Flight vs. Momentum, Unique Events", "Momentum, (GeV)", "Time, (ns)");
  }
  else
  {
    SetHistTitles(h_Mass,  "Mass, Degenerate Events",                        "Mass, (GeV)",     "Time, (ns)");
    SetHistTitles(h_TFMom, "Time of Flight vs. Momentum, Degenerate Events", "Momentum, (GeV)", "Time, (ns)");
  }

  for(unsigned int i = 0; i < vec_TFPROFCO.size(); i++)
  {
    TFCoincidenceRecord::TFCoincidence currentTFCO = vec_TFPROFCO[i].fTFCo;
    std::vector<PROFCoincidenceRecord::PROFCoincidence> currentPROFCOs = vec_TFPROFCO[i].fPROFCos;
    std::vector<std::vector<double>> currentMasses = vec_TFPROFCO[i].fMass;
    h_Multi->Fill((vec_TFPROFCO[i].getMultiplicity()<=5 ?  vec_TFPROFCO[i].getMultiplicity() : 6));
    for(unsigned int j = 0; j < currentMasses.size(); j++)
    {
      for(unsigned int k = 0; k < currentMasses[j].size(); k++)
      {
        h_Mass ->Fill(currentMasses[j][k]);
        h_TFMom->Fill(currentPROFCOs[j].getMomentum()[k], currentTFCO.fTF);
      }
    }
  }

  vec_TFPROFHist1D.push_back(h_Mass);
  vec_TFPROFHist2D.push_back(h_TFMom);

  return;
}
