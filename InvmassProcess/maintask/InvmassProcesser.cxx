#include "../config/config.h"
#include "../include/toollib.h"
#include "../include/toollib.cpp"
#include "../include/InvmassProcesser.h"

#include "TLorentzVector.h"

//init steeings



void InvmassProcesser(TString dataset, TString target)
{  
  TH1F* h_invmass[4];
  for(int i= 0; i < 4; i++)
  {
    h_invmass[i] = new TH1F(Form("h_invmass_%d",i),"",massbinnum,massmin,massmax);
  }
  TH1F* h_invmass_after[4];
  for(int i= 0; i < 4; i++)
  {
    h_invmass_after[i] = new TH1F(Form("h_invmass_after_%d",i),"",massbinnum,massmin,massmax);
  }
  std::vector<std::vector<TH1F*>> h_invmass_matter_before;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_invmass_matter_before_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_invmass_matter_before_temp.push_back(new TH1F(Form("h_invmass_matter_before_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_invmass_matter_before.push_back(h_invmass_matter_before_temp);
  }
  std::vector<std::vector<TH1F*>> h_invmass_matter_after;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_invmass_matter_after_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_invmass_matter_after_temp.push_back(new TH1F(Form("h_invmass_matter_after_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_invmass_matter_after.push_back(h_invmass_matter_after_temp);
  }
  std::vector<std::vector<TH1F*>> h_invmass_anti_before;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_invmass_anti_before_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_invmass_anti_before_temp.push_back(new TH1F(Form("h_invmass_anti_before_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_invmass_anti_before.push_back(h_invmass_anti_before_temp);
  }
  std::vector<std::vector<TH1F*>> h_invmass_anti_after;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_invmass_anti_after_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_invmass_anti_after_temp.push_back(new TH1F(Form("h_invmass_anti_after_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_invmass_anti_after.push_back(h_invmass_anti_after_temp);
  }
  std::vector<std::vector<TH1F*>> h_bkg_matter_before;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_bkg_matter_before_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_bkg_matter_before_temp.push_back(new TH1F(Form("h_bkg_matter_before_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_bkg_matter_before.push_back(h_bkg_matter_before_temp);
  }
  std::vector<std::vector<TH1F*>> h_bkg_matter_after;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_bkg_matter_after_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_bkg_matter_after_temp.push_back(new TH1F(Form("h_bkg_matter_after_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_bkg_matter_after.push_back(h_bkg_matter_after_temp);
  }
  std::vector<std::vector<TH1F*>> h_bkg_anti_before;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_bkg_anti_before_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_bkg_anti_before_temp.push_back(new TH1F(Form("h_bkg_anti_before_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_bkg_anti_before.push_back(h_bkg_anti_before_temp);
  }
  std::vector<std::vector<TH1F*>> h_bkg_anti_after;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_bkg_anti_after_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_bkg_anti_after_temp.push_back(new TH1F(Form("h_bkg_anti_after_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_bkg_anti_after.push_back(h_bkg_anti_after_temp);
  }
  std::vector<std::vector<TH1F*>> h_invmass_total_before;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_invmass_total_before_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_invmass_total_before_temp.push_back(new TH1F(Form("h_invmass_total_before_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_invmass_total_before.push_back(h_invmass_total_before_temp);
  }
  std::vector<std::vector<TH1F*>> h_invmass_total_after;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_invmass_total_after_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_invmass_total_after_temp.push_back(new TH1F(Form("h_invmass_total_after_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_invmass_total_after.push_back(h_invmass_total_after_temp);
  }
  std::vector<std::vector<TH1F*>> h_bkg_total_before;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_bkg_total_before_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_bkg_total_before_temp.push_back(new TH1F(Form("h_bkg_total_before_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_bkg_total_before.push_back(h_bkg_total_before_temp);
  }
  std::vector<std::vector<TH1F*>> h_bkg_total_after;
  for(int i = 0; i < config::nCenbins; i++)
  {
    std::vector<TH1F*> h_bkg_total_after_temp;
    for(int j = 0; j < config::nPtbins; j++)
    {
      h_bkg_total_after_temp.push_back(new TH1F(Form("h_bkg_total_after_%d_%d",i,j),"",massbinnum,massmin,massmax));
    }
    h_bkg_total_after.push_back(h_bkg_total_after_temp);
  }
  gStyle->SetOptStat(0000);
  gStyle->SetImageScaling(50.);
  ifstream initxt;
  initxt.open("../input/pathtoinputfile.txt");
  if (!initxt.is_open()) { cout << "::: Error ::: Unable to open input text file" << endl;return;}
  std::vector<TString> rootFiles;
  std::string  line = "";
  while (std::getline(initxt, line)) 
  {
    if (!line.empty()) // 检查是否为空行
    {  
      rootFiles.push_back(TString(line));
    }
  }
  initxt.close();
  cout << "number of root files: " << rootFiles.size() << endl;
  int i =1;
  for(const auto& pathtofile : rootFiles)
  {
    cout << "processing file " << i << " " << pathtofile << endl;
    i++;
    TFile* inputfile = tool::FileReader(pathtofile.Data());
    if (inputfile == NULL) { cout << "::: Error ::: inputfile not found" << endl; return; }
    TList   *keyList = inputfile->GetListOfKeys();    // 获取文件中的对象列表
    TIter   next(keyList);                     // 创建一个迭代器来遍历列表
    TKey    *key;                              // 定义一个指针，用于存储当前遍历到的对象的键
    while((key = (TKey*)next())) 
    {             
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TDirectoryFile")) continue;
      TString dirName = (TString ) key->GetName();
      TTree *TreeCandidate_1=NULL;
      TDirectory *dir = (TDirectory*)inputfile->Get(Form("%s:/%s",inputfile->GetName(),dirName.Data()));
      dir->GetObject(TreeName_1.Data(),TreeCandidate_1);
      if(!TreeCandidate_1) continue;
      else
      {
        fTreeChain = TreeCandidate_1;
        tool::SetTreeBranchAddress(fTreeChain,PiKDCascadeDecayTree_BranchAddress_int,PiKDCascadeDecayTree_BranchAddress_float,PiKDCascadeDecayTree_BranchAddress_double,PiKDCascadeDecayTree_Branches);
        Long64_t nEntries = fTreeChain->GetEntries();
        for(Long64_t i = 0; i < nEntries; ++i)
        {
          fTreeChain->GetEntry(i);
          TLorentzVector vPion, vKaon, vDeuteron, vV0, vCascade;
          vPion.SetXYZM(T1PionPx,T1PionPy,T1PionPz,PionMass);
          vKaon.SetXYZM(T1KaonPx,T1KaonPy,T1KaonPz,KaonMass);
          vDeuteron.SetXYZM(T1DeuteronPx,T1DeuteronPy,T1DeuteronPz,DeuteronMass);
          vV0 = vPion+vDeuteron;
          vCascade = vV0+vKaon;
          Double_t MassV0 = vV0.M();
          Double_t Invmass = vCascade.M();
          Float_t Pt_cascade = vCascade.Pt();
          h_invmass[T1DecayChannel-1]->Fill(Invmass);
          for(int i = 0; i < config::nCenbins; i++)
          {
            for(int j = 0; j < config::nPtbins; j++)
            {
              if((T1EventPercentile > config::Cen_bins[i]) && (T1EventPercentile <= config::Cen_bins[i+1]) && (Pt_cascade > config::Pt_bins[j]) && (Pt_cascade <= config::Pt_bins[j+1]))
              {
                // cout <<"i: "<< i << " j: " << j << endl;
                // cout<<"decaychannel" << T1DecayChannel << endl;
                if(T1DecayChannel == 1) 
                {
                  h_invmass_matter_before[i][j]->Fill(Invmass);
                  h_invmass_total_before[i][j]->Fill(Invmass);
                }
                else if(T1DecayChannel == 2)
                {
                  h_invmass_anti_before[i][j]->Fill(Invmass);
                  h_invmass_total_before[i][j]->Fill(Invmass);
                } 
                else if(T1DecayChannel == 3) 
                {
                  h_bkg_matter_before[i][j]->Fill(Invmass);
                  h_bkg_total_before[i][j]->Fill(Invmass);
                }
                else if(T1DecayChannel == 4)
                {
                  h_bkg_anti_before[i][j]->Fill(Invmass);
                  h_bkg_total_before[i][j]->Fill(Invmass);
                } 
              }
            }
          }
          if(T1TVCPA < config::cascade_cpa_min) continue;
          if(T1SVCPA < config::v0_cpa_min) continue;
          if(T1TVDCA > config::cascade_dca_max) continue;
          if(T1SVDCA > config::v0_dca_max) continue;
          if(TMath::Abs(T1KaonTPCnSigma) > config::abs_kaon_nsigmatpc_max) continue;
          if(TMath::Abs(T1PionTPCnSigma) > config::abs_pion_nsigmatpc_max) continue;
          if(TMath::Abs(T1DeuteronTPCnSigma) > config::abs_deuteron_nsigmatpc_max) continue;
          if((T1TVtoPV < config::tvtopv_min) || (T1TVtoPV > config::tvtopv_max)) continue;
          if((T1SVtoPV < config::svtopv_min) || (T1SVtoPV > config::svtopv_max)) continue;
          if((T1SVtoTV < config::svtotv_min) || (T1SVtoTV > config::svtotv_max)) continue;
          if(T1KaonTPCNcls < config::kaon_tpcncls_min) continue;
          if(T1PionTPCNcls < config::pion_tpcncls_min) continue;
          if(T1DeuteronTPCNcls < config::deuteron_tpcncls_min) continue;
          if(T1KaonDCAtoPVXY < config::kaon_dcatopv_min) continue;
          if(T1DeuteronDCAtoPVXY < config::deuteron_dcatopv_min) continue;
          if(T1PionDCAtoPVXY < config::pion_dcatopv_min) continue;
          h_invmass_after[T1DecayChannel-1]->Fill(Invmass);
          for(int i = 0; i < config::nCenbins; i++)
          {
            for(int j = 0; j < config::nPtbins; j++)
            {
              if((T1EventPercentile > config::Cen_bins[i]) && (T1EventPercentile <= config::Cen_bins[i+1]) && (Pt_cascade > config::Pt_bins[j]) && (Pt_cascade <= config::Pt_bins[j+1]))
              {
                if(T1DecayChannel == 1)
                {
                  h_invmass_matter_after[i][j]->Fill(Invmass);
                  h_invmass_total_after[i][j]->Fill(Invmass);
                } 
                else if(T1DecayChannel == 2)
                {
                  h_invmass_anti_after[i][j]->Fill(Invmass);
                  h_invmass_total_after[i][j]->Fill(Invmass);
                } 
                else if(T1DecayChannel == 3)
                {
                  h_bkg_matter_after[i][j]->Fill(Invmass);
                  h_bkg_total_after[i][j]->Fill(Invmass);
                } 
                else if(T1DecayChannel == 4)
                {
                  h_bkg_anti_after[i][j]->Fill(Invmass);
                  h_bkg_total_after[i][j]->Fill(Invmass);
                } 
              }
            }
          }
        }  
      }
    }
  }
  
  DrawSignalOnly(h_invmass[0],h_invmass[1],dataset.Data(),target.Data(),"beforecut");
  xpos += 50;
  DrawSignalOnly(h_invmass_after[0],h_invmass_after[1],dataset.Data(),target.Data(),"aftercut");
  xpos += 50;
  ypos += 100;
  DrawSignalvsBKG(h_invmass[0],h_invmass[1],h_invmass[2],h_invmass[3],dataset.Data(),target.Data(),"beforecut");
  xpos += 50;
  DrawSignalvsBKG(h_invmass_after[0],h_invmass_after[1],h_invmass_after[2],h_invmass_after[3],dataset.Data(),target.Data(),"aftercut");
  xpos += 50;
  ypos += 100;
  Draw_Split_Invmass_Spectrum(h_invmass_total_before,h_bkg_total_before,config::nCenbins,config::nPtbins,dataset.Data(),target.Data(),"tota_be");
  xpos += 50;
  Draw_Split_Invmass_Spectrum(h_invmass_total_after,h_bkg_total_after,config::nCenbins,config::nPtbins,dataset.Data(),target.Data(),"tota_af");
  xpos += 50;
  Draw_Split_Invmass_Spectrum(h_invmass_matter_before,h_bkg_matter_before,config::nCenbins,config::nPtbins,dataset.Data(),target.Data(),"matt_be");
  xpos += 50;
  Draw_Split_Invmass_Spectrum(h_invmass_matter_after,h_bkg_matter_after,config::nCenbins,config::nPtbins,dataset.Data(),target.Data(),"matt_af");
  xpos += 50;
  Draw_Split_Invmass_Spectrum(h_invmass_anti_before,h_bkg_anti_before,config::nCenbins,config::nPtbins,dataset.Data(),target.Data(),"anti_be");
  xpos += 50;
  Draw_Split_Invmass_Spectrum(h_invmass_anti_after,h_bkg_anti_after,config::nCenbins,config::nPtbins,dataset.Data(),target.Data(),"anti_af");
  cout<<"====================="<<Form("End of InvmassProcess for : %s_%s",dataset.Data(),target.Data())<<"====================="<<endl;
}


void DrawSignalOnly(TH1* h_signal,TH1* h_signal_anti,TString dataset,TString target,TString beoraf)
{
  tool::SetTitleTH1(h_signal,"counts",0.1,0.7,"InvM_{#pi+K+D} (Gev/c^{2})",0.07,0.8);
  tool::SetTitleTH1(h_signal_anti,"counts",0.1,0.7,"InvM_{#pi+K+D} (Gev/c^{2})",0.07,0.8);
  tool::SetAxisTH1(h_signal,0.0,0.0,massmin,massmax,0.04,0.04);
  tool::SetAxisTH1(h_signal_anti,0.0,0.0,massmin,massmax,0.04,0.04);
  tool::SetMarkerTH1(h_signal,"Matter",33,1.5,4,4);
  tool::SetMarkerTH1(h_signal_anti,"AntiMatter",33,1.5,2,2);
  TCanvas *canvas = tool::GetCanvas(Form("canvas_SignalOnly_%s_%s_%s",dataset.Data(),target.Data(),beoraf.Data()),xpos,ypos,cansizeX,cansizeY,0,1,0.02,0.02,0.14,0.01);
  canvas->cd();
  TPad* padleft = tool::GetPad(Form("padleft%s",dataset.Data()),0.0,0.0,0.48,1.0,0.07,0.15,0.14,0.10);
  padleft->Draw();
  padleft->cd();
  padleft->SetTicks();
  padleft->SetGridy();
  double yMax = h_signal->GetMaximum();
  h_signal->SetMaximum(yMax * 1.1); 
  h_signal->Draw("PE0");
  TLegend* legendleft = new TLegend(0.45,0.7,0.675,0.9);
  legendleft->SetBorderSize(0);
  legendleft->SetFillColor(0);
  legendleft->SetTextFont(42);
  legendleft->SetTextSize(0.05);
  legendleft->AddEntry((TObject*)0,Form("matter_total_%s",beoraf.Data()),"");
  legendleft->AddEntry(h_signal,Form("signal:%s",dataset.Data()),"P");
  legendleft->Draw();  

  canvas->cd();
  TPad* padright = tool::GetPad(Form("padright%s",dataset.Data()),0.52,0.0,0.98,1.0,0.07,0.15,0.14,0.10);
  padright->Draw();
  padright->cd();
  padright->SetTicks();
  padright->SetGridy();
  yMax = h_signal_anti->GetMaximum();
  h_signal_anti->SetMaximum(yMax * 1.1); 
  h_signal_anti->Draw("PE0");
  TLegend* legendright = new TLegend(0.45,0.7,0.675,0.9);
  legendright->SetBorderSize(0);
  legendright->SetFillColor(0);
  legendright->SetTextFont(42);
  legendright->SetTextSize(0.05);
  legendright->AddEntry((TObject*)0,Form("anti_total_%s",beoraf.Data()),"");
  legendright->AddEntry(h_signal_anti,Form("signal_anti:%s",dataset.Data()),"P");
  legendright->Draw();  
 
  TString saving_path = Form("../output/%s_%s/",dataset.Data(),target.Data());
  saving_path = tool::create_folder(saving_path);
  TString saving_name = Form("SignalOnly_total_%s.pdf",beoraf.Data());
  canvas->SaveAs(Form("%s%s",saving_path.Data(),saving_name.Data()));
}

void DrawSignalvsBKG(TH1* h_signal,TH1* h_signal_anti,TH1* h_bkg,TH1* h_bkg_anti,TString dataset,TString target,TString beoraf)
{
  tool::SetTitleTH1(h_signal,"counts",0.1,0.7,"InvM_{#pi+K+D} (Gev/c^{2})",0.07,0.8);
  tool::SetTitleTH1(h_signal_anti,"counts",0.1,0.7,"InvM_{#pi+K+D} (Gev/c^{2})",0.07,0.8);
  tool::SetTitleTH1(h_bkg,"counts",0.1,0.7,"InvM_{#pi+K+D} (Gev/c^{2})",0.07,0.8);
  tool::SetTitleTH1(h_bkg_anti,"counts",0.1,0.7,"InvM_{#pi+K+D} (Gev/c^{2})",0.07,0.8);
  tool::SetAxisTH1(h_signal,0.0,0.0,massmin,massmax,0.04,0.04);
  tool::SetAxisTH1(h_signal_anti,0.0,0.0,massmin,massmax,0.04,0.04);
  tool::SetAxisTH1(h_bkg,0.0,0.0,massmin,massmax,0.04,0.04);
  tool::SetAxisTH1(h_bkg_anti,0.0,0.0,massmin,massmax,0.04,0.04);
  tool::SetMarkerTH1(h_signal,"Matter",27,1.5,2,2);
  tool::SetMarkerTH1(h_signal_anti,"AntiMatter",27,1.5,2,2);
  tool::SetMarkerTH1(h_bkg,"",33,1.5,1,1);
  tool::SetMarkerTH1(h_bkg_anti,"",33,1.5,1,1);
  TCanvas *canvas = tool::GetCanvas(Form("c1_%s_%s_%s",dataset.Data(),target.Data(),beoraf.Data()),xpos,ypos,cansizeX,cansizeY,0,1,0.02,0.02,0.14,0.01);
  canvas->cd();
  TPad* padleft = tool::GetPad(Form("padleft%s",dataset.Data()),0.0,0.0,0.48,1.0,0.07,0.15,0.14,0.10);
  padleft->Draw();
  padleft->cd();
  padleft->SetTicks();
  padleft->SetGridy();
  double yMax = std::max(h_signal->GetMaximum(), h_bkg->GetMaximum());
  h_signal->SetMaximum(yMax * 1.1); 
  h_bkg->SetMaximum(yMax * 1.1);
  h_signal->Draw("PE0");
  h_bkg->Draw("SAMEPE0");
  TLegend* legendleft = new TLegend(0.45,0.7,0.675,0.9);
  legendleft->SetBorderSize(0);
  legendleft->SetFillColor(0);
  legendleft->SetTextFont(42);
  legendleft->SetTextSize(0.05);
  legendleft->AddEntry((TObject*)0,Form("matter_total_%s",beoraf.Data()),"");
  legendleft->AddEntry(h_signal,Form("signal:%s",dataset.Data()),"P");
  legendleft->AddEntry(h_bkg,Form("bkg:%s",dataset.Data()),"P");
  legendleft->Draw();  

  canvas->cd();
  TPad* padright = tool::GetPad(Form("padright%s",dataset.Data()),0.52,0.0,0.98,1.0,0.07,0.15,0.14,0.10);
  padright->Draw();
  padright->cd();
  padright->SetTicks();
  padright->SetGridy();
  yMax = std::max(h_signal_anti->GetMaximum(), h_bkg_anti->GetMaximum());
  h_signal_anti->SetMaximum(yMax * 1.1); 
  h_bkg_anti->SetMaximum(yMax * 1.1);
  h_signal_anti->Draw("PE0");
  h_bkg_anti->Draw("SAMEPE0");
  TLegend* legendright = new TLegend(0.45,0.7,0.675,0.9);
  legendright->SetBorderSize(0);
  legendright->SetFillColor(0);
  legendright->SetTextFont(42);
  legendright->SetTextSize(0.05);
  legendright->AddEntry((TObject*)0,Form("anti_total_%s",beoraf.Data()),"");
  legendright->AddEntry(h_signal_anti,Form("signal_anti:%s",dataset.Data()),"P");
  legendright->AddEntry(h_bkg_anti,Form("bkg_anti:%s",dataset.Data()),"P");
  legendright->Draw();  
 
  TString saving_path = Form("../output/%s_%s/",dataset.Data(),target.Data());
  saving_path = tool::create_folder(saving_path);
  TString saving_name = Form("SignalvsRotationBKG_total_%s.pdf",beoraf.Data());
  canvas->SaveAs(Form("%s%s",saving_path.Data(),saving_name.Data()));
}

void Draw_Split_Invmass_Spectrum(std::vector<std::vector<TH1F*>> h_signal, std::vector<std::vector<TH1F*>> h_bkg, int xbins, int ybins, TString dataset, TString target, TString type)
{
  //SAMPLE type : matt_be, anti_af, tota_be
  TString species;
  TString beoraf;
  int species_int = 0;
  int beoraf_int = 0;
  Int_t marker_type = 0;
  Int_t marker_color = 0;
  Int_t marker_color_bkg = 1;
  if(type(0,4) == "matt")
  {
    species = "Matter";
    species_int = 1;
    marker_color = 2;
  }   
  else if (type(0,4) == "anti")
  {
    species = "Antimatter";
    species_int = 2;
    marker_color = 9;
  }   
  else if (type(0,4) == "tota")
  {
    species = "Total";
    species_int = 3;
    marker_color = 4;
  } 
  if(type(5,2) == "be") 
  {
    beoraf = "before_cut";
    beoraf_int = 1;
    marker_type = 29;
  }  
  else if (type(5,2) == "af")
  {
    beoraf = "after_cut";
    beoraf_int = 2;
    marker_type =47; 
  }  
  for (int i = 0; i < xbins; i++)
  {
    for (int j = 0; j < ybins; j++)
    {
     tool::SetTitleTH1(h_signal[i][j],"counts",0.06,0.8,"InvM_{#pi+K+D} (Gev/c^{2})",0.05,0.8);
     tool::SetTitleTH1(h_bkg[i][j],   "counts",0.06,0.8,"InvM_{#pi+K+D} (Gev/c^{2})",0.05,0.8);
     tool::SetAxisTH1(h_signal[i][j],0.0,0.0,massmin,massmax,0.04,0.04);
     tool::SetAxisTH1(h_bkg[i][j],   0.0,0.0,massmin,massmax,0.04,0.04);
     tool::SetMarkerTH1(h_signal[i][j],Form("%s_%s",species.Data(),beoraf.Data()), marker_type, 1.5, marker_color, marker_color);
     tool::SetMarkerTH1(h_bkg[i][j],   "", marker_type, 1.5, marker_color_bkg, marker_color_bkg);
    }
  }
  TCanvas* canvas = tool::GetCanvas(Form("can_split_spec_%s_%s_%s",dataset.Data(),target.Data(),type.Data()),xpos,ypos,1800,1000,0,1,0.02,0.02,0.14,0.01);
  canvas->Divide(ybins,xbins);
  TLegend *leg[xbins][ybins];
  for (int i = 0; i < xbins; i++)
  {
    for (int j = 0; j < ybins; j++)
    {
      canvas->cd((i*ybins)+j+1);
      gPad->SetGridy();
      double yMax = std::max(h_signal[i][j]->GetMaximum(), h_bkg[i][j]->GetMaximum());
      h_signal[i][j]->SetMaximum(yMax * 1.1); 
      h_bkg[i][j]->SetMaximum(yMax * 1.1);
      h_signal[i][j]->Draw("PE0"); 
      h_bkg[i][j]->Draw("SAMEPE0");
      leg[i][j] = new TLegend(0.5,0.7,0.9,0.9);
      leg[i][j]->SetBorderSize(0);
      leg[i][j]->SetFillColor(0);
      leg[i][j]->SetTextFont(42);
      leg[i][j]->SetTextSize(0.05);
      leg[i][j]->AddEntry((TObject*)0,Form("%s",dataset.Data()),"");
      leg[i][j]->AddEntry((TObject*)0,"bkg_Rot_D #pm 5^{#circ}","");
      leg[i][j]->AddEntry(h_signal[i][j],"SIGNAL","lep");
      leg[i][j]->AddEntry(h_bkg[i][j],   "BKG","lep");
      leg[i][j]->Draw("SAME");
      tool::drawMyTextNDC(0.5,0.6,0.05,Form("Cen: %.2f-%.2f",config::Cen_bins[i],config::Cen_bins[i+1]),kRed);
      tool::drawMyTextNDC(0.5,0.5,0.05,Form("Pt: %.2f-%.2f",config::Pt_bins[j],config::Pt_bins[j+1]),kRed);
    }
  }
  TString saving_path = Form("../output/%s_%s/",dataset.Data(),target.Data());
  saving_path = tool::create_folder(saving_path);
  TString saving_name = Form("Invmass_split_spec_%s_%s_%s",dataset.Data(),target.Data(),type.Data());
  canvas->SaveAs(Form("%s%s.pdf",saving_path.Data(),saving_name.Data()));
}
