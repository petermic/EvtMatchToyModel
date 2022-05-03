#include "RejectionSampler.h"
#include "TRandom.h"

void samplertest(){
  TFile* f = TFile::Open("/raid5/data/wangj/home/phoDcorrRun2018/mcsub/rootfiles/hydjet_hydjet__PbPb_MC/savehist.root");
  TH1F* nD = (TH1F*)f->Get("hnG");
  TH1F* njt = (TH1F*)f->Get("hnjt");
  RejectionSampler rD(nD);
  RejectionSampler rjt(njt);
  TRandomMixMax rpythia;
  TRandomMixMax rphi;
  double nD_coeff = 1./6.;
  double njt_coeff = 1./8.;
  double dphi_signal_mean = 0.;
  double dphi_signal_width = 0.2;

  unsigned int n_evts = 100000000;
  // e = embed, m = mix
  TH1F* npairs_ee = new TH1F("npairs_ee","n_pairs (Pythia+HYDJET D, Pythia+HYDJET jet)",41,-20,20);
  TH1F* npairs_em = new TH1F("npairs_em","n_pairs (Pythia+HYDJET D, HYDJET jet)",41,-20,20);
  TH1F* npairs_me = new TH1F("npairs_me","n_pairs (HYDJET D, Pythia+HYDJET jet)",41,-20,20);
  TH1F* npairs_mm = new TH1F("npairs_mm","n_pairs (HYDJET D, HYDJET jet)",41,-20,20);
  TH1F* npairs_pythia = new TH1F("npairs_pythia","n_pairs (Pythia)",41,-20,20);

  TH1F* dphi_ee = new TH1F("dphi_ee","dphi (Pythia+HYDJET D, Pythia+HYDJET jet)",100,0,M_PI);
  TH1F* dphi_em = new TH1F("dphi_em","dphi (Pythia+HYDJET D, HYDJET jet)",100,0,M_PI);
  TH1F* dphi_me = new TH1F("dphi_me","dphi (HYDJET D, Pythia+HYDJET jet)",100,0,M_PI);
  TH1F* dphi_mm = new TH1F("dphi_mm","dphi (HYDJET D, HYDJET jet)",100,0,M_PI);
  TH1F* dphi_pythia = new TH1F("dphi_pythia","dphi (Pythia)",100,0,M_PI);
  TH1F* dphi_ratio = new TH1F("dphi_ratio","Subtracted / Pythia",100,0,M_PI);

  int nD_hydjet = 0;
  int nD_total = 0;
  int njt_hydjet = 0;
  int njt_total = 0;

  for(int i=0;i<n_evts;i++){
    int nD_mb = rD.getIntValue();
    int njt_mb = rjt.getIntValue();
    int nD_pythia = round(rpythia.Exp(nD_coeff));
    int njt_pythia = round(rpythia.Exp(njt_coeff));
    int nD_mix = rD.getIntValue();
    int njt_mix = rjt.getIntValue();

    nD_hydjet += nD_mb;
    nD_total += (nD_mb + nD_pythia);
    njt_hydjet += njt_mb;
    njt_total += (njt_mb + njt_pythia);
    
    int np_ee = (nD_mb + nD_pythia) * (njt_mb + njt_pythia);
    int np_em = (nD_mb + nD_pythia) * njt_mix;
    int np_me = nD_mix * (njt_mb + njt_pythia);
    int np_mm = nD_mix * njt_mix;
    int np_pythia = nD_pythia * njt_pythia;

    for(int j=0;j<np_ee;j++){
      if(j<np_pythia) dphi_ee->Fill(rphi.Gaus(dphi_signal_mean,dphi_signal_width));
      else dphi_ee->Fill(rphi.Uniform(M_PI));
    }
    for(int j=0;j<np_em;j++){
      dphi_em->Fill(rphi.Uniform(M_PI));
    }
    for(int j=0;j<np_me;j++){
      dphi_me->Fill(rphi.Uniform(M_PI));
    }
    for(int j=0;j<np_mm;j++){
      dphi_mm->Fill(rphi.Uniform(M_PI));
    }
    for(int j=0;j<np_pythia;j++){
      dphi_pythia->Fill(rphi.Gaus(dphi_signal_mean,dphi_signal_width));
    }

    npairs_ee->Fill(np_ee);
    npairs_em->Fill(np_em);
    npairs_me->Fill(np_me);
    npairs_mm->Fill(np_mm);
    npairs_pythia->Fill(np_pythia);
  }

  dphi_ee->Sumw2();
  dphi_em->Sumw2();
  dphi_me->Sumw2();
  dphi_mm->Sumw2();

  TH1F* dphi_sub = new TH1F("dphi_sub","dphi (subtracted)",100,0,M_PI);
  dphi_sub->Add(dphi_ee,dphi_em,1.,-1.);
  dphi_sub->Add(dphi_me,-1.);
  dphi_sub->Add(dphi_mm);

  gStyle->SetOptStat(0);
/*
  TCanvas* c = new TCanvas("c","c",800,600);
  c->SetLogy();

  npairs_ee->SetMinimum(10);
  npairs_ee->SetMaximum(pow(10,7));
  npairs_ee->SetTitle("");

  npairs_ee->SetLineColor(kBlack);
  npairs_em->SetLineColor(kBlue);
  npairs_me->SetLineColor(kGreen+2);
  npairs_mm->SetLineColor(kViolet);
  npairs_pythia->SetLineColor(kRed);

  npairs_ee->Draw("E");
  npairs_em->Draw("Esame");
  npairs_me->Draw("Esame");
  npairs_mm->Draw("Esame");
  npairs_pythia->Draw("Esame");

  TLegend* legend = new TLegend(0.8,0.8,1.,1.);
  legend->AddEntry(npairs_ee);
  legend->AddEntry(npairs_em);
  legend->AddEntry(npairs_me);
  legend->AddEntry(npairs_mm);
  legend->Draw("same");
*/
  TCanvas* d = new TCanvas("d","d",800,600);
  
  dphi_ee->SetMinimum(0);
  dphi_ee->SetTitle("");

  dphi_ee->SetLineColor(kOrange+2);
  dphi_em->SetLineColor(kBlue);
  dphi_me->SetLineColor(kGreen+2);
  dphi_mm->SetLineColor(kViolet);
  dphi_pythia->SetLineColor(kRed);
  dphi_sub->SetLineColor(kBlack);

  dphi_ee->Draw("E");
  dphi_em->Draw("Esame");
  dphi_me->Draw("Esame");
  dphi_mm->Draw("Esame");
  dphi_pythia->Draw("Esame");
  dphi_sub->Draw("Esame");

  TLegend* dlegend = new TLegend(0.6,0.8,1.,1.);
  dlegend->AddEntry(dphi_ee,"dphi (Pythia+HYDJET D, Pythia+HYDJET jet)");
  dlegend->AddEntry(dphi_em);
  dlegend->AddEntry(dphi_me);
  dlegend->AddEntry(dphi_mm);
  dlegend->AddEntry(dphi_pythia);
  dlegend->AddEntry(dphi_sub);
  dlegend->Draw("same");

  TCanvas* e = new TCanvas("e","e",800,600);
  dphi_pythia->Sumw2();
  dphi_ratio->Divide(dphi_sub,dphi_pythia);
  dphi_ratio->SetMinimum(0.);
  dphi_ratio->SetMaximum(2.);
  dphi_ratio->Draw("E");

  std::cout << "D ratio: HYDJET/total = " << (float)njt_hydjet/(float)njt_total << std::endl;
  std::cout << "jet ratio: HYDJET/total = " << (float)nD_hydjet/(float)nD_total << std::endl;
}
