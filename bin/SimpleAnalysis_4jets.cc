#include <string>
#include <iostream>
#include <vector>
#include <chrono>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h"

#include "Analysis/Core/interface/Analysis.h"
#include <boost/filesystem.hpp>

using namespace analysis;
using namespace analysis::tools;

void printUsage(const char* bin_name) {
   std::cout << "Usage: " << bin_name << " [rootFileList_name] [num_jets]" << std::endl;
}

// =============================================================================================   
int main(int argc, char* argv[])
{
   bool isMC = true;
   bool isbbbb = true;
   std::string inputList = "rootFileList.txt";
   std::string outputRoot = "histograms.root";
   std::string json = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
   unsigned int num_jets = 4;
   int bWP = 1;
   float btagcut[3] = { 0.46, 0.84, 0.92 };
   // Cuts                                         // <<<<===== CMSDAS
   float ptmin[4]   = { 100.0, 100.0, 40.0, 30. };
   float etamax[4]  = {   2.2,   2.2,  2.2, 2.4 };
   float btagmin[4] = { btagcut[bWP], btagcut[bWP], btagcut[bWP], btagcut[0] };
   float nonbtag    = 0.46;
   float dRmin      = 1.;
   float detamax    = 1.55;
   
   if (argc >= 2) {
      inputList = argv[1];
   }
   if (argc == 3) {
      int arg;
      try {
         arg = std::stoi(argv[2]);
      } catch (...) {
         printUsage(argv[0]);
         return 1;
      }
      if (arg <= 0) {
         std::cout << "num_jets must be greater than 0." << std::endl;
         return 1;
      } else {
         num_jets = arg;
      }
   } else if (argc > 3) {
      std::cout << "Too many arguments." << std::endl;
      printUsage(argv[0]);
      return 1;
   }
   
   outputRoot = boost::filesystem::basename(inputList) + "_hist.root";
   
   
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   // Input files list
   Analysis analysis(inputList);
   
   //analysis.addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsPuppiReapplyJEC");
   analysis.addTree<Jet> ("Jets","MssmHbb/Events/selectedUpdatedPatJetsPuppi");
   
   std::vector<std::string> triggerObjects;
   triggerObjects.push_back("hltL1sDoubleJetC100");
   triggerObjects.push_back("hltDoubleJetsC100");
   triggerObjects.push_back("hltBTagCaloCSVp014DoubleWithMatching");
   triggerObjects.push_back("hltDoublePFJetsC100");
   triggerObjects.push_back("hltDoublePFJetsC100MaxDeta1p6");

   for ( auto & obj : triggerObjects )
      analysis.addTree<TriggerObject> (obj, Form("MssmHbb/Events/selectedPatTrigger/%s", obj.c_str()));
   
   analysis.triggerResults("MssmHbb/Events/TriggerResults");
   std::string hltPath = "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v";
   
   
   if( !isMC ) analysis.processJsonFile(json);
   
   TFile hout(outputRoot.c_str(),"recreate");
   
   std::map<std::string, TH1F*> h1;
   h1["total"]    = new TH1F("total", "", 1, 0., 1.);
   h1["totalSelected"]    = new TH1F("totalSelected", "", 1, 0., 1.);
   h1["n"]        = new TH1F("n" , "" , 30, 0, 30);
   h1["n_csv"]    = new TH1F("n_csv" , "" , 30, 0, 30);
   h1["n_ptmin20"]= new TH1F("n_ptmin20" , "" , 30, 0, 30);
   h1["n_ptmin20_csv"] = new TH1F("n_ptmin20_csv" , "" , 30, 0, 30);
   h1["n_ptmin30_csv"] = new TH1F("n_ptmin30_csv" , "" , 30, 0, 30);
   for ( unsigned int i = 0 ; i < num_jets ; ++i )
   {
      if (i < 2)
      {
         h1[Form("pt_%i",i)]         = new TH1F(Form("pt_%i",i) , "" , 100, 0, 1000);
         h1[Form("eta_%i",i)]        = new TH1F(Form("eta_%i",i) , "" , 100, -5, 5);
         h1[Form("phi_%i",i)]        = new TH1F(Form("phi_%i",i) , "" , 100, -4, 4);
         h1[Form("btag_%i",i)]       = new TH1F(Form("btag_%i",i) , "" , 100, 0, 1);
         
         h1[Form("pt_%i_csv",i)]     = new TH1F(Form("pt_%i_csv",i) , "" , 100, 0, 1000);
         h1[Form("eta_%i_csv",i)]    = new TH1F(Form("eta_%i_csv",i) , "" , 100, -5, 5);
         h1[Form("phi_%i_csv",i)]    = new TH1F(Form("phi_%i_csv",i) , "" , 100, -4, 4);
         h1[Form("btag_%i_csv",i)]   = new TH1F(Form("btag_%i_csv",i) , "" , 100, 0, 1);
      }
      else
      {
         h1[Form("pt_%i",i)]         = new TH1F(Form("pt_%i",i) , "" ,  50, 0, 200);
         h1[Form("eta_%i",i)]        = new TH1F(Form("eta_%i",i) , "" , 100, -5, 5);
         h1[Form("phi_%i",i)]        = new TH1F(Form("phi_%i",i) , "" , 100, -4, 4);
         h1[Form("btag_%i",i)]       = new TH1F(Form("btag_%i",i) , "" , 100, 0, 1);
         
         h1[Form("pt_%i_csv",i)]     = new TH1F(Form("pt_%i_csv",i) , "" , 50, 0, 200);
         h1[Form("eta_%i_csv",i)]    = new TH1F(Form("eta_%i_csv",i) , "" , 100, -5, 5);
         h1[Form("phi_%i_csv",i)]    = new TH1F(Form("phi_%i_csv",i) , "" , 100, -4, 4);
         h1[Form("btag_%i_csv",i)]   = new TH1F(Form("btag_%i_csv",i) , "" , 100, 0, 1);
      }
   }
   h1["m12"]     = new TH1F("m12"     , "" , 50, 0, 1000);
   h1["m12_csv"] = new TH1F("m12_csv" , "" , 50, 0, 1000);
   
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events." << std::endl;
   
   auto start_time = std::chrono::steady_clock::now();
   
   // Cut flow
   // 0: triggered events
   // 1: 3+ idloose jets
   // 2: kinematics
   // 3: matched to online
   // 4: delta R
   // 5: delta eta
   // 6: btag (bbnb)
   int nsel[10] = { };
   int nmatch[10] = { };
   
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {
      int njets = 0;
      int njets_csv_pt20 = 0;
      int njets_csv_pt30 = 0;
      bool goodEvent = true;
      
      if ( i > 0 && i%100000==0 ) std::cout << i << " events processed!" << std::endl;
      
      analysis.event(i);
      h1["total"]->Fill(0.5);
      if (! isMC )
      {
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
      }
      
      int triggerFired = analysis.triggerResult(hltPath);
      if ( !triggerFired ) continue;
      
      ++nsel[0];
      
      // Jets - std::shared_ptr< Collection<Jet> >
      auto slimmedJets = analysis.collection<Jet>("Jets");
      std::vector<Jet *> selectedJets;
      for ( int j = 0 ; j < slimmedJets->size() ; ++j )
      {
         Jet* jet = &slimmedJets->at(j);
         if (jet->idLoose())
            selectedJets.push_back(jet);
      }
      if ( selectedJets.size() < num_jets ) continue;
      ++nsel[1];
      
      // Fill histograms
      for ( int j = 0 ; j < (int)selectedJets.size() ; ++j )
      {
         if ( selectedJets[j]->pt() < 20. ) continue;
         ++njets;
      }
      h1["n"] -> Fill(selectedJets.size());
      h1["n_ptmin20"] -> Fill(njets);
      h1["m12"] -> Fill((selectedJets[0]->p4() + selectedJets[1]->p4()).M());
      for ( unsigned int j = 0; j < num_jets; ++j )
      {
         Jet* jet = selectedJets[j];
         h1[Form("pt_%i",j)]   -> Fill(jet->pt());
         h1[Form("eta_%i",j)]  -> Fill(jet->eta());
         h1[Form("phi_%i",j)]  -> Fill(jet->phi());
         h1[Form("btag_%i",j)] -> Fill(jet->btag());
      }
      
      // Kinematic selection
      for ( unsigned int j = 0; j < num_jets; ++j )
      {
         Jet* jet = selectedJets[j];
         if ( jet->pt() < ptmin[j] || fabs(jet->eta()) > etamax[j] )
         {
            goodEvent = false;
            break;
         }
      }
      if ( ! goodEvent ) continue;
      ++nsel[2];
      
      // Delta R selection
      for ( unsigned int j1 = 0; j1 < num_jets - 1; ++j1 )
      {
         const Jet& jet1 = *selectedJets[j1];
         for ( unsigned int j2 = j1 + 1; j2 < num_jets; ++j2 )
         {
            const Jet& jet2 = *selectedJets[j2];
            if ( jet1.deltaR(jet2) < dRmin )
            {
               goodEvent = false;
               break;
            }
         }
      }
      if ( ! goodEvent ) continue;
      ++nsel[3];
      
      // Delta eta selection - 2 leading jets
      if ( fabs(selectedJets[0]->eta() - selectedJets[1]->eta()) > detamax ) continue;
      ++nsel[4];
      
      for (unsigned int j = 0; j < num_jets; ++j)
      {
         Jet* jet = selectedJets[j];
         if (isbbbb)
         {
            if ( j < num_jets && jet->btag() < btagmin[j] )
            //if ( j < num_jets && (jet->btag("btag_deepb") + jet->btag("btag_deepbb")) < btagmin[j] )
            {
               goodEvent = false;
               break;
            }
         }
         else
         {
            if ( j < num_jets && jet->btag() > nonbtag )
            //if ( j < num_jets && (jet->btag("btag_deepb") + jet->btag("btag_deepbb")) < btagmin[j] )
            {
               goodEvent = false;
               break;
            }
         }
      }
      if ( ! goodEvent ) continue;
      
      ++nsel[5];
        
      // match offline to online
      analysis.match<Jet, TriggerObject>("Jets", triggerObjects, 0.5);
      // Are the leading to jets matched?
      bool matched[10] = {true,true,true,true,true,true,true,true,true,true};
      for ( unsigned int j = 0; j < 2; ++j )
      {
         Jet* jet = selectedJets[j];
         for ( size_t io = 0; io < triggerObjects.size() ; ++io )
         {       
            if ( ! jet->matched(triggerObjects[io]) ) matched[io] = false;
         }
      }
      
      for ( size_t io = 0; io < triggerObjects.size() ; ++io )
      {
         if (matched[io])
            ++nmatch[io];
         else
         {
            goodEvent = false;
            break;
         }
      }
      if ( ! goodEvent ) continue;

      ++nsel[6];
      
      // Fill histograms of passed bbnb btagging selection
      for ( int j = 0 ; j < (int)selectedJets.size() ; ++j )
      {
         if (selectedJets[j]->pt() >= 20.)
            ++njets_csv_pt20;
         if (selectedJets[j]->pt() >= 30.)
            ++njets_csv_pt30;
      }
      h1["totalSelected"]->Fill(0.5);
      h1["n_csv"] -> Fill(selectedJets.size());
      h1["n_ptmin20_csv"] -> Fill(njets_csv_pt20);
      h1["n_ptmin30_csv"] -> Fill(njets_csv_pt30);
      for ( unsigned int j = 0; j < num_jets; ++j )
      {
         Jet* jet = selectedJets[j];
         h1[Form("pt_%i_csv",j)]   -> Fill(jet->pt());
         h1[Form("eta_%i_csv",j)]  -> Fill(jet->eta());
         h1[Form("phi_%i_csv",j)]  -> Fill(jet->phi());
         h1[Form("btag_%i_csv",j)] -> Fill(jet->btag());
      }

      //if ( !isbbbb ) h1["m12_csv"] -> Fill((selectedJets[0]->p4() + selectedJets[1]->p4()).M());
      h1["m12_csv"] -> Fill((selectedJets[0]->p4() + selectedJets[1]->p4()).M());
   }
   
   auto end_time = std::chrono::steady_clock::now();
   std::cout << "Analysis took " << 
      std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
      << " ms." << std::endl;
   float eff = (double)nsel[6] / analysis.size();
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
// PRINT OUTS
   
   // Cut flow
   // 0: triggered events
   // 1: 3+ idloose jets
   // 2: matched to online
   // 3: kinematics
   // 4: delta R
   // 5: delta eta
   // 6: btag (bbnb)
   
   double fracAbs[10];
   double fracRel[10];
   std::string cuts[10];
   cuts[0] = "Triggered";
   cuts[1] = "Triple idloose-jet";
   cuts[2] = "Triple jet kinematics";
   cuts[3] = "Delta R(i;j)";
   cuts[4] = "Delta eta(j1;j2)";
   if (isbbbb)
      cuts[5] = "btagged (bbbb)";
   else
      cuts[5] = "btagged (bbbnb)";
   cuts[6] = "Matched to online j1;j2";
   
   printf ("%-23s  %10s  %10s  %10s \n", std::string("Cut flow").c_str(), std::string("# events").c_str(), std::string("absolute").c_str(), std::string("relative").c_str() ); 
   for ( int i = 0; i < 7; ++i )
   {
      fracAbs[i] = double(nsel[i])/nsel[0];
      if ( i>0 )
         fracRel[i] = double(nsel[i])/nsel[i-1];
      else
         fracRel[i] = 1.;
      printf ("%-23s  %10d  %10.3f  %10.3f \n", cuts[i].c_str(), nsel[i], fracAbs[i], fracRel[i] ); 
   }
   /*
   // CSV output
   printf ("%-23s , %10s , %10s , %10s \n", std::string("Cut flow").c_str(), std::string("# events").c_str(), std::string("absolute").c_str(), std::string("relative").c_str() ); 
   for ( int i = 0; i < 7; ++i )
      printf ("%-23s , %10d , %10.3f , %10.3f \n", cuts[i].c_str(), nsel[i], fracAbs[i], fracRel[i] );
   */

   // Trigger objects counts   
   std::cout << std::endl;
   printf ("%-40s  %10s \n", std::string("Trigger object").c_str(), std::string("# events").c_str() ); 
   for ( size_t io = 0; io < triggerObjects.size() ; ++io )
   {
      printf ("%-40s  %10d \n", triggerObjects[io].c_str(), nmatch[io] ); 
   }
   
   // Efficiency
   std::cout << std::endl;
   std::cout << "Efficiency: " << eff << std::endl;
}

