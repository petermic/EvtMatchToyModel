#include "TRandom.h"
#include "RejectionSampler.h"

void evtmatch_MBOnMB(){
  TFile* f = TFile::Open("/raid5/data/wangj/home/phoDcorrRun2018/mcsub/rootfiles/hydjet_hydjet__PbPb_MC/savehist.root");
  TH1F* nD = (TH1F*)f->Get("hnG");
  TH1F* njt = (TH1F*)f->Get("hnjt");
  RejectionSampler rD(nD);
  RejectionSampler rjt(njt);
  TRandomMixMax rsignal;
  TRandomMixMax rphi;
  float D_signal_ratio = 0.434291;
  float jet_signal_ratio = 0.887814;
  float dphi_signal_mean = 0.;
  float dphi_signal_width = 0.2;

  unsigned int n_evts = 100000000;
  // M = MinBias (signal + background), X = Mixed event
  TH1F* npairs_mm = new TH1F("npairs_mm","n_pairs (MB D, MB jet)",20,0,20);
  TH1F* npairs_mx = new TH1F("npairs_mx","n_pairs (MB D, Mixed jet)",20,0,20);
  TH1F* npairs_xm = new TH1F("npairs_xm","n_pairs (Mixed D, MB jet)",20,0,20);
  TH1F* npairs_xx = new TH1F("npairs_xx","n_pairs (Mixed D, Mixed jet)",20,0,20);
  TH1F* npairs_signal = new TH1F("npairs_signal","n_pairs (signal)",20,0,20);

  TH1F* dphi_mm = new TH1F("dphi_mm","dphi (MB D, MB jet)",100,0,M_PI);
  TH1F* dphi_mx = new TH1F("dphi_mx","dphi (MB D, Mixed jet)",100,0,M_PI);
  TH1F* dphi_xm = new TH1F("dphi_xm","dphi (Mixed D, MB jet)",100,0,M_PI);
  TH1F* dphi_xx = new TH1F("dphi_xx","dphi (Mixed D, Mixed jet)",100,0,M_PI);
  TH1F* dphi_signal = new TH1F("dphi_signal","dphi (signal)",100,0,M_PI);
  TH1F* dphi_ratio = new TH1F("dphi_ratio","dphi subtracted / signal",100,0,M_PI);

  for(int i=0;i<n_evts;i++){
    // # of D and jet in triggered MB event
    int nD_mb = rD.getIntValue();
    int njt_mb = rjt.getIntValue();
    // # of D and jet in other MB event, for mixing
    int nD_mix = rD.getIntValue();
    int njt_mix = rjt.getIntValue();

    int np_mm = nD_mb * njt_mb;
    int np_mx = nD_mb * njt_mix;
    int np_xm = nD_mix * njt_mb;
    int np_xx = nD_mix * njt_mix;

    for(int j=0;j<np_mm;j++){
      if(rsignal.Uniform(0,1)>D_signal_ratio && rsignal.Uniform(0,1)>jet_signal_ratio){
        dphi_mm->Fill(rphi.Gaus(dphi_signal_mean,dphi_signal_width));
        dphi_signal->Fill(rphi.Gaus(dphi_signal_mean,dphi_signal_width));
      }
      else dphi_mm->Fill(rphi.Uniform(M_PI));
    }
    for(int j=0;j<np_mx;j++){
      dphi_mx->Fill(rphi.Uniform(M_PI));
    }
    for(int j=0;j<np_xm;j++){
      dphi_xm->Fill(rphi.Uniform(M_PI));
    }
    for(int j=0;j<np_xx;j++){
      if(rsignal.Uniform(0,1)>D_signal_ratio && rsignal.Uniform(0,1)>jet_signal_ratio) dphi_xx->Fill(rphi.Gaus(dphi_signal_mean,dphi_signal_width));
      else dphi_xx->Fill(rphi.Uniform(M_PI));
    }

    npairs_mm->Fill(np_mm);
    npairs_mx->Fill(np_mx);
    npairs_xm->Fill(np_xm);
    npairs_xx->Fill(np_xx);
  }

  dphi_mx->Scale(jet_signal_ratio);
  dphi_xm->Scale(D_signal_ratio);
  dphi_xx->Scale((jet_signal_ratio*D_signal_ratio));

  dphi_mm->Sumw2();
  dphi_mx->Sumw2();
  dphi_xm->Sumw2();
  dphi_xx->Sumw2();

  TH1F* dphi_sub = new TH1F("dphi_sub","dphi (subtracted)",100,0,M_PI);
  dphi_sub->Add(dphi_mm,dphi_mx,1.,-1.);
  dphi_sub->Add(dphi_xm,-1.);
  dphi_sub->Add(dphi_xx);
  dphi_sub->Scale(1./(1.+jet_signal_ratio*D_signal_ratio));

  dphi_sub->Sumw2();

  gStyle->SetOptStat(0);
/*
  TCanvas* c = new TCanvas("c","c",800,600);
  c->SetLogy();

  npairs_mm->SetMinimum(-100);
  npairs_mm->SetMaximum(pow(10,7));
  npairs_mm->SetTitle("");

  npairs_mm->SetLineColor(kOrange+2);
  npairs_mx->SetLineColor(kBlue);
  npairs_xm->SetLineColor(kGreen+2);
  npairs_xx->SetLineColor(kViolet);
  npairs_signal->SetLineColor(kRed);

  npairs_mm->Draw("E");
  npairs_mx->Draw("Esame");
  npairs_xm->Draw("Esame");
  npairs_xx->Draw("Esame");
  npairs_signal->Draw("Esame");

  TLegend* legend = new TLegend(0.8,0.8,1.,1.);
  legend->AddEntry(npairs_mm);
  legend->AddEntry(npairs_mx);
  legend->AddEntry(npairs_xm);
  legend->AddEntry(npairs_xx);
  legend->Draw("same");
*/
  TCanvas* d = new TCanvas("d","d",800,600);
  
  dphi_mm->SetMinimum(0.);
  dphi_mm->SetTitle("");

  dphi_mm->SetLineColor(kOrange+2);
  dphi_mx->SetLineColor(kGreen+2);
  dphi_xm->SetLineColor(kBlue);
  dphi_xx->SetLineColor(kViolet);
  dphi_signal->SetLineColor(kRed);
  dphi_sub->SetLineColor(kBlack);

  dphi_mm->Draw("E");
  dphi_mx->Draw("Esame");
  dphi_xm->Draw("Esame");
  dphi_xx->Draw("Esame");
  dphi_signal->Draw("Esame");
  dphi_sub->Draw("Esame");

  TLegend* dlegend = new TLegend(0.7,0.8,1.,1.);
  dlegend->AddEntry(dphi_mm,"dphi (MB D, MB jet)");
  dlegend->AddEntry(dphi_mx);
  dlegend->AddEntry(dphi_xm);
  dlegend->AddEntry(dphi_xx);
  dlegend->AddEntry(dphi_signal);
  dlegend->AddEntry(dphi_sub);
  dlegend->Draw("same");

  TCanvas* e = new TCanvas("e","e",800,600);
  dphi_ratio->Divide(dphi_sub,dphi_signal);
  dphi_ratio->SetMinimum(0.);
  dphi_ratio->SetMaximum(2.);
  dphi_ratio->Draw("E");
}
