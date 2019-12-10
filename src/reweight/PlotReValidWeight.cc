#include "NSyst.h"

#include "EventSummaryHelper.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace neut;
using namespace neut::rew;

int main(int argc, char const *argv[]) {

  if (argc != 4) {
    std::cout << "[ERROR]: Expected to be passed: File1.root:title "
                 "File2.root:title Output.pdf"
              << std::endl;
    return 1;
  }

  std::string file1name = argv[1];
  std::string title1 = "";
  if (file1name.find_first_of(':') != std::string::npos) {
    title1 = file1name.substr(file1name.find_first_of(':') + 1);
    std::cout << "title1 " << title1 << ", file1name " << file1name << ", "
              << file1name.find_first_of(':') << std::endl;
    file1name = file1name.substr(0, file1name.find_first_of(':'));
  }

  std::string file2name = argv[2];
  std::string title2 = "";
  if (file2name.find_first_of(':') != std::string::npos) {
    title2 = file2name.substr(file2name.find_first_of(':') + 1);
    std::cout << "title2 " << title2 << ", file2name " << file2name << ", "
              << file2name.find_first_of(':') << std::endl;
    file2name = file2name.substr(0, file2name.find_first_of(':'));
  }

  std::vector<std::vector<std::map<int, HistReadBlob> > > AllTheHists1;
  std::vector<std::vector<std::map<int, HistReadBlob> > > AllTheHists2;

  TFile *f1 = TFile::Open(file1name.c_str(), "READ");
  TFile *f2 = TFile::Open(file2name.c_str(), "READ");

  std::vector<NSyst_t> Dials;
  std::vector<std::string> DialNames;

  Dials.push_back(kXSecTwkDial_MaCCQE);
  Dials.push_back(kXSecTwkDial_MaRES);
  Dials.push_back(kXSecTwkDial_CA5RES);
  Dials.push_back(kXSecTwkDial_BgSclRES);
  Dials.push_back(kXSecTwkDial_BgSclLMCPiBarRES);
  // Dials.push_back(kCascTwkDial_FrAbs_pi);
  // Dials.push_back(kCascTwkDial_FrInelLow_pi);
  // Dials.push_back(kCascTwkDial_FrInelHigh_pi);
  // Dials.push_back(kCascTwkDial_FrCExLow_pi);
  // Dials.push_back(kCascTwkDial_FrCExHigh_pi);
  // Dials.push_back(kCascTwkDial_FrPiProd_pi);

  std::map<NSyst_t, std::vector<int> > RelevantModes;

  RelevantModes[kXSecTwkDial_MaCCQE].push_back(1);
  RelevantModes[kXSecTwkDial_MaCCQE].push_back(-1);

  RelevantModes[kXSecTwkDial_MaRES].push_back(11);
  RelevantModes[kXSecTwkDial_MaRES].push_back(-11);
  RelevantModes[kXSecTwkDial_MaRES].push_back(12);
  RelevantModes[kXSecTwkDial_MaRES].push_back(-12);
  RelevantModes[kXSecTwkDial_MaRES].push_back(13);
  RelevantModes[kXSecTwkDial_MaRES].push_back(-13);

  RelevantModes[kXSecTwkDial_CA5RES] = RelevantModes[kXSecTwkDial_MaRES];
  RelevantModes[kXSecTwkDial_BgSclRES] = RelevantModes[kXSecTwkDial_MaRES];

  std::map<int, std::string> ModeNames;
  ModeNames[0] = "All Channels";
  ModeNames[1] = "True CCQE";
  ModeNames[11] = "True RES p#pi^{+}";
  ModeNames[12] = "True RES p#pi^{0}";
  ModeNames[13] = "True RES n#pi^{+}";

  size_t NDials = Dials.size();

  AllTheHists1.resize(NDials);
  AllTheHists2.resize(NDials);

  for (size_t d_it = 0; d_it < NDials; ++d_it) {
    DialNames.push_back(NSyst::AsString(Dials[d_it]));

    TDirectory *dial_dir_f1 = f1->GetDirectory(DialNames[d_it].c_str());
    if (!dial_dir_f1) {
      std::cout << "[ERROR]: Failed to read directory: " << DialNames[d_it]
                << " from " << file1name << std::endl;
      exit(1);
    }
    AllTheHists1[d_it].resize(7);
    for (int i = -100; i < 100; ++i) {
      std::stringstream ss("");
      ss << "mode_" << i;
      TDirectory *mode_dir = dial_dir_f1->GetDirectory(ss.str().c_str());

      if (!mode_dir) {
        continue;
      }
      for (int dv = -3; dv < 4; ++dv) {
        ss.str("");
        ss << "tweak_" << (dv < 0 ? "m" : "") << std::abs(dv);
        AllTheHists1[d_it][dv + 3].insert(
            std::make_pair(i, HistReadBlob(mode_dir, ss.str())));
      }
    }

    TDirectory *dial_dir_f2 = f2->GetDirectory(DialNames[d_it].c_str());
    if (!dial_dir_f2) {

      std::string dir = DialNames[d_it];
      if (Dials[d_it] == kXSecTwkDial_MaRES) {
        dir = "MaNFFRES";
      }
      dial_dir_f2 = f2->GetDirectory(dir.c_str());
      if (!dial_dir_f2) {
        std::cout << "[ERROR]: Failed to read directory: " << DialNames[d_it]
                  << " or " << dir << " from " << file2name << std::endl;
        exit(1);
      }
    }
    AllTheHists2[d_it].resize(7);
    for (int i = -100; i < 100; ++i) {
      std::stringstream ss("");
      ss << "mode_" << i;
      TDirectory *mode_dir = dial_dir_f2->GetDirectory(ss.str().c_str());

      if (!mode_dir) {
        continue;
      }
      for (int dv = -3; dv < 4; ++dv) {
        ss.str("");
        ss << "tweak_" << (dv < 0 ? "m" : "") << std::abs(dv);
        AllTheHists2[d_it][dv + 3].insert(
            std::make_pair(i, HistReadBlob(mode_dir, ss.str())));
      }
    }
  }

  TCanvas c1("c1", "", 600, 400);

  c1.Print((std::string(argv[3]) + "[").c_str());

  gStyle->SetOptStat(false);

  for (size_t d_it = 0; d_it < NDials; ++d_it) {

    for (std::vector<int>::iterator modeit = RelevantModes[Dials[d_it]].begin();
         modeit != RelevantModes[Dials[d_it]].end(); ++modeit) {
      int mode = *modeit;

      if (!AllTheHists1[d_it][6].count(mode)) {
        continue;
      }

      //*** Q2 ***/
      TPad *p1 = new TPad("p1", "", 0, 0.451, 0.5, 0.85);
      p1->AppendPad();
      p1->cd();

      TH1D *LeftNomCloneQ2 =
          static_cast<TH1D *>(AllTheHists1[d_it][3][mode].Q2->Clone());
      TH1D *RightNomCloneQ2 =
          static_cast<TH1D *>(AllTheHists2[d_it][3][mode].Q2->Clone());

      AllTheHists1[d_it][0][mode].Q2->SetLineColor(kRed);
      AllTheHists1[d_it][1][mode].Q2->SetLineColor(kBlue);
      AllTheHists1[d_it][2][mode].Q2->SetLineColor(kMagenta);
      AllTheHists1[d_it][3][mode].Q2->SetLineColor(kBlack);
      AllTheHists1[d_it][4][mode].Q2->SetLineColor(kGreen - 2);
      AllTheHists1[d_it][5][mode].Q2->SetLineColor(kGray + 2);
      AllTheHists1[d_it][6][mode].Q2->SetLineColor(kAzure);

      AllTheHists1[d_it][6][mode].Q2->GetYaxis()->SetNdivisions(505);
      AllTheHists1[d_it][6][mode].Q2->GetYaxis()->SetTitleSize(0.05);
      AllTheHists1[d_it][6][mode].Q2->GetYaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][6][mode].Q2->GetXaxis()->SetNdivisions(505);
      AllTheHists1[d_it][6][mode].Q2->GetXaxis()->SetTitleSize(0.05);
      AllTheHists1[d_it][6][mode].Q2->GetXaxis()->SetLabelSize(0.05);

      AllTheHists1[d_it][6][mode].Q2->DrawClone("EHIST");
      AllTheHists1[d_it][0][mode].Q2->DrawClone("EHISTSAME");
      AllTheHists1[d_it][1][mode].Q2->DrawClone("EHISTSAME");
      AllTheHists1[d_it][2][mode].Q2->DrawClone("EHISTSAME");
      AllTheHists1[d_it][3][mode].Q2->DrawClone("EHISTSAME");
      AllTheHists1[d_it][4][mode].Q2->DrawClone("EHISTSAME");
      AllTheHists1[d_it][5][mode].Q2->DrawClone("EHISTSAME");

      c1.cd();

      TPad *p1rat = new TPad("p1rat", "", 0, 0, 0.5, 0.449);
      p1rat->AppendPad();
      p1rat->cd();

      AllTheHists1[d_it][6][mode].Q2->Divide(AllTheHists1[d_it][3][mode].Q2);
      AllTheHists1[d_it][0][mode].Q2->Divide(AllTheHists1[d_it][3][mode].Q2);
      AllTheHists1[d_it][1][mode].Q2->Divide(AllTheHists1[d_it][3][mode].Q2);
      AllTheHists1[d_it][2][mode].Q2->Divide(AllTheHists1[d_it][3][mode].Q2);
      AllTheHists1[d_it][4][mode].Q2->Divide(AllTheHists1[d_it][3][mode].Q2);
      AllTheHists1[d_it][5][mode].Q2->Divide(AllTheHists1[d_it][3][mode].Q2);
      AllTheHists1[d_it][3][mode].Q2->Divide(AllTheHists1[d_it][3][mode].Q2);

      double max;

      max = -std::numeric_limits<double>::max();
      for (size_t hit = 0; hit < AllTheHists1[d_it].size(); ++hit) {
        max = std::max(max, AllTheHists1[d_it][hit][mode].Q2->GetMaximum());
      }
      AllTheHists1[d_it][6][mode].Q2->GetYaxis()->SetRangeUser(0, max * 1.1);

      AllTheHists1[d_it][6][mode].Q2->Draw("HIST");
      AllTheHists1[d_it][0][mode].Q2->Draw("HISTSAME");
      AllTheHists1[d_it][1][mode].Q2->Draw("HISTSAME");
      AllTheHists1[d_it][2][mode].Q2->Draw("HISTSAME");
      AllTheHists1[d_it][3][mode].Q2->Draw("HISTSAME");
      AllTheHists1[d_it][4][mode].Q2->Draw("HISTSAME");
      AllTheHists1[d_it][5][mode].Q2->Draw("HISTSAME");

      c1.cd();

      TPad *p2 = new TPad("p2", "", 0.5, 0.451, 1, 0.85);
      p2->AppendPad();
      p2->cd();

      AllTheHists2[d_it][0][mode].Q2->SetLineColor(kRed);
      AllTheHists2[d_it][1][mode].Q2->SetLineColor(kBlue);
      AllTheHists2[d_it][2][mode].Q2->SetLineColor(kMagenta);
      AllTheHists2[d_it][3][mode].Q2->SetLineColor(kBlack);
      AllTheHists2[d_it][4][mode].Q2->SetLineColor(kGreen - 2);
      AllTheHists2[d_it][5][mode].Q2->SetLineColor(kGray + 2);
      AllTheHists2[d_it][6][mode].Q2->SetLineColor(kAzure);

      AllTheHists2[d_it][6][mode].Q2->GetYaxis()->SetNdivisions(505);
      AllTheHists2[d_it][6][mode].Q2->GetYaxis()->SetTitleSize(0.05);
      AllTheHists2[d_it][6][mode].Q2->GetYaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][6][mode].Q2->GetXaxis()->SetNdivisions(505);
      AllTheHists2[d_it][6][mode].Q2->GetXaxis()->SetTitleSize(0.05);
      AllTheHists2[d_it][6][mode].Q2->GetXaxis()->SetLabelSize(0.05);

      AllTheHists2[d_it][6][mode].Q2->DrawClone("EHIST");
      AllTheHists2[d_it][0][mode].Q2->DrawClone("EHISTSAME");
      AllTheHists2[d_it][1][mode].Q2->DrawClone("EHISTSAME");
      AllTheHists2[d_it][2][mode].Q2->DrawClone("EHISTSAME");
      AllTheHists2[d_it][3][mode].Q2->DrawClone("EHISTSAME");
      AllTheHists2[d_it][4][mode].Q2->DrawClone("EHISTSAME");
      AllTheHists2[d_it][5][mode].Q2->DrawClone("EHISTSAME");

      c1.cd();

      TPad *p2rat = new TPad("p2rat", "", 0.5, 0, 1, 0.449);
      p2rat->AppendPad();
      p2rat->cd();

      AllTheHists2[d_it][6][mode].Q2->Divide(AllTheHists2[d_it][3][mode].Q2);
      AllTheHists2[d_it][0][mode].Q2->Divide(AllTheHists2[d_it][3][mode].Q2);
      AllTheHists2[d_it][1][mode].Q2->Divide(AllTheHists2[d_it][3][mode].Q2);
      AllTheHists2[d_it][2][mode].Q2->Divide(AllTheHists2[d_it][3][mode].Q2);
      AllTheHists2[d_it][4][mode].Q2->Divide(AllTheHists2[d_it][3][mode].Q2);
      AllTheHists2[d_it][5][mode].Q2->Divide(AllTheHists2[d_it][3][mode].Q2);
      AllTheHists2[d_it][3][mode].Q2->Divide(AllTheHists2[d_it][3][mode].Q2);

      max = -std::numeric_limits<double>::max();
      for (size_t hit = 0; hit < AllTheHists2[d_it].size(); ++hit) {
        max = std::max(max, AllTheHists2[d_it][hit][mode].Q2->GetMaximum());
      }
      AllTheHists2[d_it][6][mode].Q2->GetYaxis()->SetRangeUser(0, max * 1.1);

      AllTheHists2[d_it][6][mode].Q2->Draw("HIST");
      AllTheHists2[d_it][0][mode].Q2->Draw("HISTSAME");
      AllTheHists2[d_it][1][mode].Q2->Draw("HISTSAME");
      AllTheHists2[d_it][2][mode].Q2->Draw("HISTSAME");
      AllTheHists2[d_it][3][mode].Q2->Draw("HISTSAME");
      AllTheHists2[d_it][4][mode].Q2->Draw("HISTSAME");
      AllTheHists2[d_it][5][mode].Q2->Draw("HISTSAME");

      RightNomCloneQ2->Divide(LeftNomCloneQ2);
      RightNomCloneQ2->SetLineColor(kBlack);
      RightNomCloneQ2->SetLineStyle(2);
      RightNomCloneQ2->Draw("HISTSAME");

      c1.cd();

      TLegend leg(0.1, 0.85, 0.9, 0.95);
      leg.SetBorderSize(0);
      leg.SetNColumns(4);

      leg.AddEntry(AllTheHists2[d_it][0][mode].Q2, "-3 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][1][mode].Q2, "-2 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][2][mode].Q2, "-1 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][3][mode].Q2, "Nominal", "l");
      leg.AddEntry(AllTheHists2[d_it][4][mode].Q2, "+1 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][5][mode].Q2, "+2 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][6][mode].Q2, "+3 #sigma", "l");
      leg.AddEntry(RightNomCloneQ2, (title2 + "/" + title1 + "Nominal").c_str(),
                   "l");
      leg.Draw();

      TLatex ltx;

      std::stringstream ss("");
      ss << "Mode: " << ModeNames[std::abs(mode)]
         << ", Dial: " << DialNames[d_it];

      std::cout << "Drawing " << ss.str() << std::endl;

      ltx.SetTextAlign(22);
      ltx.DrawLatexNDC(0.5, 0.975, ss.str().c_str());
      ltx.DrawLatexNDC(0.25, 0.84, title1.c_str());
      ltx.DrawLatexNDC(0.75, 0.84, title2.c_str());

      c1.Print(argv[3]);
      c1.Clear();

      /** W **/

      TH1D *LeftNomCloneW =
          static_cast<TH1D *>(AllTheHists1[d_it][3][mode].W->Clone());
      TH1D *RightNomCloneW =
          static_cast<TH1D *>(AllTheHists2[d_it][3][mode].W->Clone());

      p1 = new TPad("p1", "", 0, 0.451, 0.5, 0.85);
      p1->AppendPad();
      p1->cd();

      AllTheHists1[d_it][0][mode].W->SetLineColor(kRed);
      AllTheHists1[d_it][1][mode].W->SetLineColor(kBlue);
      AllTheHists1[d_it][2][mode].W->SetLineColor(kMagenta);
      AllTheHists1[d_it][3][mode].W->SetLineColor(kBlack);
      AllTheHists1[d_it][4][mode].W->SetLineColor(kGreen - 2);
      AllTheHists1[d_it][5][mode].W->SetLineColor(kGray + 2);
      AllTheHists1[d_it][6][mode].W->SetLineColor(kAzure);

      AllTheHists1[d_it][6][mode].W->GetYaxis()->SetNdivisions(505);
      AllTheHists1[d_it][6][mode].W->GetYaxis()->SetTitleSize(0.05);
      AllTheHists1[d_it][6][mode].W->GetYaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][6][mode].W->GetXaxis()->SetNdivisions(505);
      AllTheHists1[d_it][6][mode].W->GetXaxis()->SetTitleSize(0.05);
      AllTheHists1[d_it][6][mode].W->GetXaxis()->SetLabelSize(0.05);

      AllTheHists1[d_it][6][mode].W->DrawClone("EHIST");
      AllTheHists1[d_it][0][mode].W->DrawClone("EHISTSAME");
      AllTheHists1[d_it][1][mode].W->DrawClone("EHISTSAME");
      AllTheHists1[d_it][2][mode].W->DrawClone("EHISTSAME");
      AllTheHists1[d_it][3][mode].W->DrawClone("EHISTSAME");
      AllTheHists1[d_it][4][mode].W->DrawClone("EHISTSAME");
      AllTheHists1[d_it][5][mode].W->DrawClone("EHISTSAME");

      c1.cd();

      p1rat = new TPad("p1rat", "", 0, 0, 0.5, 0.449);
      p1rat->AppendPad();
      p1rat->cd();

      AllTheHists1[d_it][6][mode].W->Divide(AllTheHists1[d_it][3][mode].W);
      AllTheHists1[d_it][0][mode].W->Divide(AllTheHists1[d_it][3][mode].W);
      AllTheHists1[d_it][1][mode].W->Divide(AllTheHists1[d_it][3][mode].W);
      AllTheHists1[d_it][2][mode].W->Divide(AllTheHists1[d_it][3][mode].W);
      AllTheHists1[d_it][4][mode].W->Divide(AllTheHists1[d_it][3][mode].W);
      AllTheHists1[d_it][5][mode].W->Divide(AllTheHists1[d_it][3][mode].W);
      AllTheHists1[d_it][3][mode].W->Divide(AllTheHists1[d_it][3][mode].W);

      max = -std::numeric_limits<double>::max();
      for (size_t hit = 0; hit < AllTheHists1[d_it].size(); ++hit) {
        max = std::max(max, AllTheHists1[d_it][hit][mode].W->GetMaximum());
      }
      AllTheHists1[d_it][6][mode].W->GetYaxis()->SetRangeUser(0, max * 1.1);

      AllTheHists1[d_it][6][mode].W->Draw("HIST");
      AllTheHists1[d_it][0][mode].W->Draw("HISTSAME");
      AllTheHists1[d_it][1][mode].W->Draw("HISTSAME");
      AllTheHists1[d_it][2][mode].W->Draw("HISTSAME");
      AllTheHists1[d_it][3][mode].W->Draw("HISTSAME");
      AllTheHists1[d_it][4][mode].W->Draw("HISTSAME");
      AllTheHists1[d_it][5][mode].W->Draw("HISTSAME");

      c1.cd();

      p2 = new TPad("p2", "", 0.5, 0.451, 1, 0.85);
      p2->AppendPad();
      p2->cd();

      AllTheHists2[d_it][0][mode].W->SetLineColor(kRed);
      AllTheHists2[d_it][1][mode].W->SetLineColor(kBlue);
      AllTheHists2[d_it][2][mode].W->SetLineColor(kMagenta);
      AllTheHists2[d_it][3][mode].W->SetLineColor(kBlack);
      AllTheHists2[d_it][4][mode].W->SetLineColor(kGreen - 2);
      AllTheHists2[d_it][5][mode].W->SetLineColor(kGray + 2);
      AllTheHists2[d_it][6][mode].W->SetLineColor(kAzure);

      AllTheHists2[d_it][6][mode].W->GetYaxis()->SetNdivisions(505);
      AllTheHists2[d_it][6][mode].W->GetYaxis()->SetTitleSize(0.05);
      AllTheHists2[d_it][6][mode].W->GetYaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][6][mode].W->GetXaxis()->SetNdivisions(505);
      AllTheHists2[d_it][6][mode].W->GetXaxis()->SetTitleSize(0.05);
      AllTheHists2[d_it][6][mode].W->GetXaxis()->SetLabelSize(0.05);

      AllTheHists2[d_it][6][mode].W->DrawClone("EHIST");
      AllTheHists2[d_it][0][mode].W->DrawClone("EHISTSAME");
      AllTheHists2[d_it][1][mode].W->DrawClone("EHISTSAME");
      AllTheHists2[d_it][2][mode].W->DrawClone("EHISTSAME");
      AllTheHists2[d_it][3][mode].W->DrawClone("EHISTSAME");
      AllTheHists2[d_it][4][mode].W->DrawClone("EHISTSAME");
      AllTheHists2[d_it][5][mode].W->DrawClone("EHISTSAME");

      c1.cd();

      p2rat = new TPad("p2rat", "", 0.5, 0, 1, 0.449);
      p2rat->AppendPad();
      p2rat->cd();

      AllTheHists2[d_it][6][mode].W->Divide(AllTheHists2[d_it][3][mode].W);
      AllTheHists2[d_it][0][mode].W->Divide(AllTheHists2[d_it][3][mode].W);
      AllTheHists2[d_it][1][mode].W->Divide(AllTheHists2[d_it][3][mode].W);
      AllTheHists2[d_it][2][mode].W->Divide(AllTheHists2[d_it][3][mode].W);
      AllTheHists2[d_it][4][mode].W->Divide(AllTheHists2[d_it][3][mode].W);
      AllTheHists2[d_it][5][mode].W->Divide(AllTheHists2[d_it][3][mode].W);
      AllTheHists2[d_it][3][mode].W->Divide(AllTheHists2[d_it][3][mode].W);

      max = -std::numeric_limits<double>::max();
      for (size_t hit = 0; hit < AllTheHists2[d_it].size(); ++hit) {
        max = std::max(max, AllTheHists2[d_it][hit][mode].W->GetMaximum());
      }
      AllTheHists2[d_it][6][mode].W->GetYaxis()->SetRangeUser(0, max * 1.1);

      AllTheHists2[d_it][6][mode].W->Draw("HIST");
      AllTheHists2[d_it][0][mode].W->Draw("HISTSAME");
      AllTheHists2[d_it][1][mode].W->Draw("HISTSAME");
      AllTheHists2[d_it][2][mode].W->Draw("HISTSAME");
      AllTheHists2[d_it][3][mode].W->Draw("HISTSAME");
      AllTheHists2[d_it][4][mode].W->Draw("HISTSAME");
      AllTheHists2[d_it][5][mode].W->Draw("HISTSAME");

      RightNomCloneW->Divide(LeftNomCloneW);
      RightNomCloneW->SetLineColor(kBlack);
      RightNomCloneW->SetLineStyle(2);
      RightNomCloneW->Draw("HISTSAME");

      c1.cd();

      leg.Clear();
      leg.SetBorderSize(0);
      leg.SetNColumns(4);

      leg.AddEntry(AllTheHists2[d_it][0][mode].W, "-3 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][1][mode].W, "-2 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][2][mode].W, "-1 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][3][mode].W, "Nominal", "l");
      leg.AddEntry(AllTheHists2[d_it][4][mode].W, "+1 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][5][mode].W, "+2 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][6][mode].W, "+3 #sigma", "l");
      leg.AddEntry(RightNomCloneW, (title2 + "/" + title1 + "Nominal").c_str(),
                   "l");
      leg.Draw();

      std::cout << "Drawing " << ss.str() << std::endl;

      ltx.SetTextAlign(22);
      ltx.DrawLatexNDC(0.5, 0.975, ss.str().c_str());
      ltx.DrawLatexNDC(0.25, 0.84, title1.c_str());
      ltx.DrawLatexNDC(0.75, 0.84, title2.c_str());

      c1.Print(argv[3]);
      c1.Clear();

      /** TPi **/

      TH1D *LeftNomCloneTpi =
          static_cast<TH1D *>(AllTheHists1[d_it][3][mode].Tpi->Clone());
      TH1D *RightNomCloneTpi =
          static_cast<TH1D *>(AllTheHists2[d_it][3][mode].Tpi->Clone());

      p1 = new TPad("p1", "", 0, 0.451, 0.5, 0.85);
      p1->AppendPad();
      p1->cd();

      AllTheHists1[d_it][0][mode].Tpi->SetLineColor(kRed);
      AllTheHists1[d_it][1][mode].Tpi->SetLineColor(kBlue);
      AllTheHists1[d_it][2][mode].Tpi->SetLineColor(kMagenta);
      AllTheHists1[d_it][3][mode].Tpi->SetLineColor(kBlack);
      AllTheHists1[d_it][4][mode].Tpi->SetLineColor(kGreen - 2);
      AllTheHists1[d_it][5][mode].Tpi->SetLineColor(kGray + 2);
      AllTheHists1[d_it][6][mode].Tpi->SetLineColor(kAzure);

      AllTheHists1[d_it][6][mode].Tpi->GetYaxis()->SetNdivisions(505);
      AllTheHists1[d_it][6][mode].Tpi->GetYaxis()->SetTitleSize(0.05);
      AllTheHists1[d_it][6][mode].Tpi->GetYaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][6][mode].Tpi->GetXaxis()->SetNdivisions(505);
      AllTheHists1[d_it][6][mode].Tpi->GetXaxis()->SetTitleSize(0.05);
      AllTheHists1[d_it][6][mode].Tpi->GetXaxis()->SetLabelSize(0.05);

      AllTheHists1[d_it][6][mode].Tpi->DrawClone("EHIST");
      AllTheHists1[d_it][0][mode].Tpi->DrawClone("EHISTSAME");
      AllTheHists1[d_it][1][mode].Tpi->DrawClone("EHISTSAME");
      AllTheHists1[d_it][2][mode].Tpi->DrawClone("EHISTSAME");
      AllTheHists1[d_it][3][mode].Tpi->DrawClone("EHISTSAME");
      AllTheHists1[d_it][4][mode].Tpi->DrawClone("EHISTSAME");
      AllTheHists1[d_it][5][mode].Tpi->DrawClone("EHISTSAME");

      c1.cd();

      p1rat = new TPad("p1rat", "", 0, 0, 0.5, 0.449);
      p1rat->AppendPad();
      p1rat->cd();

      AllTheHists1[d_it][6][mode].Tpi->Divide(AllTheHists1[d_it][3][mode].Tpi);
      AllTheHists1[d_it][0][mode].Tpi->Divide(AllTheHists1[d_it][3][mode].Tpi);
      AllTheHists1[d_it][1][mode].Tpi->Divide(AllTheHists1[d_it][3][mode].Tpi);
      AllTheHists1[d_it][2][mode].Tpi->Divide(AllTheHists1[d_it][3][mode].Tpi);
      AllTheHists1[d_it][4][mode].Tpi->Divide(AllTheHists1[d_it][3][mode].Tpi);
      AllTheHists1[d_it][5][mode].Tpi->Divide(AllTheHists1[d_it][3][mode].Tpi);
      AllTheHists1[d_it][3][mode].Tpi->Divide(AllTheHists1[d_it][3][mode].Tpi);

      max = -std::numeric_limits<double>::max();
      for (size_t hit = 0; hit < AllTheHists1[d_it].size(); ++hit) {
        max = std::max(max, AllTheHists1[d_it][hit][mode].Tpi->GetMaximum());
      }
      AllTheHists1[d_it][6][mode].Tpi->GetYaxis()->SetRangeUser(0, max * 1.1);

      AllTheHists1[d_it][6][mode].Tpi->Draw("HIST");
      AllTheHists1[d_it][0][mode].Tpi->Draw("HISTSAME");
      AllTheHists1[d_it][1][mode].Tpi->Draw("HISTSAME");
      AllTheHists1[d_it][2][mode].Tpi->Draw("HISTSAME");
      AllTheHists1[d_it][3][mode].Tpi->Draw("HISTSAME");
      AllTheHists1[d_it][4][mode].Tpi->Draw("HISTSAME");
      AllTheHists1[d_it][5][mode].Tpi->Draw("HISTSAME");

      c1.cd();

      p2 = new TPad("p2", "", 0.5, 0.451, 1, 0.85);
      p2->AppendPad();
      p2->cd();

      AllTheHists2[d_it][0][mode].Tpi->SetLineColor(kRed);
      AllTheHists2[d_it][1][mode].Tpi->SetLineColor(kBlue);
      AllTheHists2[d_it][2][mode].Tpi->SetLineColor(kMagenta);
      AllTheHists2[d_it][3][mode].Tpi->SetLineColor(kBlack);
      AllTheHists2[d_it][4][mode].Tpi->SetLineColor(kGreen - 2);
      AllTheHists2[d_it][5][mode].Tpi->SetLineColor(kGray + 2);
      AllTheHists2[d_it][6][mode].Tpi->SetLineColor(kAzure);

      AllTheHists2[d_it][6][mode].Tpi->GetYaxis()->SetNdivisions(505);
      AllTheHists2[d_it][6][mode].Tpi->GetYaxis()->SetTitleSize(0.05);
      AllTheHists2[d_it][6][mode].Tpi->GetYaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][6][mode].Tpi->GetXaxis()->SetNdivisions(505);
      AllTheHists2[d_it][6][mode].Tpi->GetXaxis()->SetTitleSize(0.05);
      AllTheHists2[d_it][6][mode].Tpi->GetXaxis()->SetLabelSize(0.05);

      AllTheHists2[d_it][6][mode].Tpi->DrawClone("EHIST");
      AllTheHists2[d_it][0][mode].Tpi->DrawClone("EHISTSAME");
      AllTheHists2[d_it][1][mode].Tpi->DrawClone("EHISTSAME");
      AllTheHists2[d_it][2][mode].Tpi->DrawClone("EHISTSAME");
      AllTheHists2[d_it][3][mode].Tpi->DrawClone("EHISTSAME");
      AllTheHists2[d_it][4][mode].Tpi->DrawClone("EHISTSAME");
      AllTheHists2[d_it][5][mode].Tpi->DrawClone("EHISTSAME");

      c1.cd();

      p2rat = new TPad("p2rat", "", 0.5, 0, 1, 0.449);
      p2rat->AppendPad();
      p2rat->cd();

      AllTheHists2[d_it][6][mode].Tpi->Divide(AllTheHists2[d_it][3][mode].Tpi);
      AllTheHists2[d_it][0][mode].Tpi->Divide(AllTheHists2[d_it][3][mode].Tpi);
      AllTheHists2[d_it][1][mode].Tpi->Divide(AllTheHists2[d_it][3][mode].Tpi);
      AllTheHists2[d_it][2][mode].Tpi->Divide(AllTheHists2[d_it][3][mode].Tpi);
      AllTheHists2[d_it][4][mode].Tpi->Divide(AllTheHists2[d_it][3][mode].Tpi);
      AllTheHists2[d_it][5][mode].Tpi->Divide(AllTheHists2[d_it][3][mode].Tpi);
      AllTheHists2[d_it][3][mode].Tpi->Divide(AllTheHists2[d_it][3][mode].Tpi);

      max = -std::numeric_limits<double>::max();
      for (size_t hit = 0; hit < AllTheHists2[d_it].size(); ++hit) {
        max = std::max(max, AllTheHists2[d_it][hit][mode].Tpi->GetMaximum());
      }
      AllTheHists2[d_it][6][mode].Tpi->GetYaxis()->SetRangeUser(0, max * 1.1);

      AllTheHists2[d_it][6][mode].Tpi->Draw("HIST");
      AllTheHists2[d_it][0][mode].Tpi->Draw("HISTSAME");
      AllTheHists2[d_it][1][mode].Tpi->Draw("HISTSAME");
      AllTheHists2[d_it][2][mode].Tpi->Draw("HISTSAME");
      AllTheHists2[d_it][3][mode].Tpi->Draw("HISTSAME");
      AllTheHists2[d_it][4][mode].Tpi->Draw("HISTSAME");
      AllTheHists2[d_it][5][mode].Tpi->Draw("HISTSAME");

      RightNomCloneTpi->Divide(LeftNomCloneTpi);
      RightNomCloneTpi->SetLineColor(kBlack);
      RightNomCloneTpi->SetLineStyle(2);
      RightNomCloneTpi->Draw("HISTSAME");

      c1.cd();

      leg.Clear();
      leg.SetBorderSize(0);
      leg.SetNColumns(4);

      leg.AddEntry(AllTheHists2[d_it][0][mode].Tpi, "-3 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][1][mode].Tpi, "-2 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][2][mode].Tpi, "-1 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][3][mode].Tpi, "Nominal", "l");
      leg.AddEntry(AllTheHists2[d_it][4][mode].Tpi, "+1 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][5][mode].Tpi, "+2 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][6][mode].Tpi, "+3 #sigma", "l");
      leg.AddEntry(RightNomCloneTpi,
                   (title2 + "/" + title1 + "Nominal").c_str(), "l");
      leg.Draw();

      std::cout << "Drawing " << ss.str() << std::endl;

      ltx.SetTextAlign(22);
      ltx.DrawLatexNDC(0.5, 0.975, ss.str().c_str());
      ltx.DrawLatexNDC(0.25, 0.84, title1.c_str());
      ltx.DrawLatexNDC(0.75, 0.84, title2.c_str());

      c1.Print(argv[3]);
      c1.Clear();

      /** ENu **/

      TH1D *LeftNomCloneEnu =
          static_cast<TH1D *>(AllTheHists1[d_it][3][mode].Enu->Clone());
      TH1D *RightNomCloneEnu =
          static_cast<TH1D *>(AllTheHists2[d_it][3][mode].Enu->Clone());

      p1 = new TPad("p1", "", 0, 0.451, 0.5, 0.85);
      p1->AppendPad();
      p1->cd();

      AllTheHists1[d_it][0][mode].Enu->SetLineColor(kRed);
      AllTheHists1[d_it][1][mode].Enu->SetLineColor(kBlue);
      AllTheHists1[d_it][2][mode].Enu->SetLineColor(kMagenta);
      AllTheHists1[d_it][3][mode].Enu->SetLineColor(kBlack);
      AllTheHists1[d_it][4][mode].Enu->SetLineColor(kGreen - 2);
      AllTheHists1[d_it][5][mode].Enu->SetLineColor(kGray + 2);
      AllTheHists1[d_it][6][mode].Enu->SetLineColor(kAzure);

      AllTheHists1[d_it][6][mode].Enu->GetYaxis()->SetNdivisions(505);
      AllTheHists1[d_it][6][mode].Enu->GetYaxis()->SetTitleSize(0.05);
      AllTheHists1[d_it][6][mode].Enu->GetYaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][6][mode].Enu->GetXaxis()->SetNdivisions(505);
      AllTheHists1[d_it][6][mode].Enu->GetXaxis()->SetTitleSize(0.05);
      AllTheHists1[d_it][6][mode].Enu->GetXaxis()->SetLabelSize(0.05);

      AllTheHists1[d_it][6][mode].Enu->DrawClone("EHIST");
      AllTheHists1[d_it][0][mode].Enu->DrawClone("EHISTSAME");
      AllTheHists1[d_it][1][mode].Enu->DrawClone("EHISTSAME");
      AllTheHists1[d_it][2][mode].Enu->DrawClone("EHISTSAME");
      AllTheHists1[d_it][3][mode].Enu->DrawClone("EHISTSAME");
      AllTheHists1[d_it][4][mode].Enu->DrawClone("EHISTSAME");
      AllTheHists1[d_it][5][mode].Enu->DrawClone("EHISTSAME");

      c1.cd();

      p1rat = new TPad("p1rat", "", 0, 0, 0.5, 0.449);
      p1rat->AppendPad();
      p1rat->cd();

      AllTheHists1[d_it][6][mode].Enu->Divide(AllTheHists1[d_it][3][mode].Enu);
      AllTheHists1[d_it][0][mode].Enu->Divide(AllTheHists1[d_it][3][mode].Enu);
      AllTheHists1[d_it][1][mode].Enu->Divide(AllTheHists1[d_it][3][mode].Enu);
      AllTheHists1[d_it][2][mode].Enu->Divide(AllTheHists1[d_it][3][mode].Enu);
      AllTheHists1[d_it][4][mode].Enu->Divide(AllTheHists1[d_it][3][mode].Enu);
      AllTheHists1[d_it][5][mode].Enu->Divide(AllTheHists1[d_it][3][mode].Enu);
      AllTheHists1[d_it][3][mode].Enu->Divide(AllTheHists1[d_it][3][mode].Enu);

      max = -std::numeric_limits<double>::max();
      for (size_t hit = 0; hit < AllTheHists1[d_it].size(); ++hit) {
        max = std::max(max, AllTheHists1[d_it][hit][mode].Enu->GetMaximum());
      }
      AllTheHists1[d_it][6][mode].Enu->GetYaxis()->SetRangeUser(0, max * 1.1);

      AllTheHists1[d_it][6][mode].Enu->Draw("HIST");
      AllTheHists1[d_it][0][mode].Enu->Draw("HISTSAME");
      AllTheHists1[d_it][1][mode].Enu->Draw("HISTSAME");
      AllTheHists1[d_it][2][mode].Enu->Draw("HISTSAME");
      AllTheHists1[d_it][3][mode].Enu->Draw("HISTSAME");
      AllTheHists1[d_it][4][mode].Enu->Draw("HISTSAME");
      AllTheHists1[d_it][5][mode].Enu->Draw("HISTSAME");

      c1.cd();

      p2 = new TPad("p2", "", 0.5, 0.451, 1, 0.85);
      p2->AppendPad();
      p2->cd();

      AllTheHists2[d_it][0][mode].Enu->SetLineColor(kRed);
      AllTheHists2[d_it][1][mode].Enu->SetLineColor(kBlue);
      AllTheHists2[d_it][2][mode].Enu->SetLineColor(kMagenta);
      AllTheHists2[d_it][3][mode].Enu->SetLineColor(kBlack);
      AllTheHists2[d_it][4][mode].Enu->SetLineColor(kGreen - 2);
      AllTheHists2[d_it][5][mode].Enu->SetLineColor(kGray + 2);
      AllTheHists2[d_it][6][mode].Enu->SetLineColor(kAzure);

      AllTheHists2[d_it][6][mode].Enu->GetYaxis()->SetNdivisions(505);
      AllTheHists2[d_it][6][mode].Enu->GetYaxis()->SetTitleSize(0.05);
      AllTheHists2[d_it][6][mode].Enu->GetYaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][6][mode].Enu->GetXaxis()->SetNdivisions(505);
      AllTheHists2[d_it][6][mode].Enu->GetXaxis()->SetTitleSize(0.05);
      AllTheHists2[d_it][6][mode].Enu->GetXaxis()->SetLabelSize(0.05);

      AllTheHists2[d_it][6][mode].Enu->DrawClone("EHIST");
      AllTheHists2[d_it][0][mode].Enu->DrawClone("EHISTSAME");
      AllTheHists2[d_it][1][mode].Enu->DrawClone("EHISTSAME");
      AllTheHists2[d_it][2][mode].Enu->DrawClone("EHISTSAME");
      AllTheHists2[d_it][3][mode].Enu->DrawClone("EHISTSAME");
      AllTheHists2[d_it][4][mode].Enu->DrawClone("EHISTSAME");
      AllTheHists2[d_it][5][mode].Enu->DrawClone("EHISTSAME");

      c1.cd();

      p2rat = new TPad("p2rat", "", 0.5, 0, 1, 0.449);
      p2rat->AppendPad();
      p2rat->cd();

      AllTheHists2[d_it][6][mode].Enu->Divide(AllTheHists2[d_it][3][mode].Enu);
      AllTheHists2[d_it][0][mode].Enu->Divide(AllTheHists2[d_it][3][mode].Enu);
      AllTheHists2[d_it][1][mode].Enu->Divide(AllTheHists2[d_it][3][mode].Enu);
      AllTheHists2[d_it][2][mode].Enu->Divide(AllTheHists2[d_it][3][mode].Enu);
      AllTheHists2[d_it][4][mode].Enu->Divide(AllTheHists2[d_it][3][mode].Enu);
      AllTheHists2[d_it][5][mode].Enu->Divide(AllTheHists2[d_it][3][mode].Enu);
      AllTheHists2[d_it][3][mode].Enu->Divide(AllTheHists2[d_it][3][mode].Enu);

      max = -std::numeric_limits<double>::max();
      for (size_t hit = 0; hit < AllTheHists2[d_it].size(); ++hit) {
        max = std::max(max, AllTheHists2[d_it][hit][mode].Enu->GetMaximum());
      }
      AllTheHists2[d_it][6][mode].Enu->GetYaxis()->SetRangeUser(0, max * 1.1);

      AllTheHists2[d_it][6][mode].Enu->Draw("HIST");
      AllTheHists2[d_it][0][mode].Enu->Draw("HISTSAME");
      AllTheHists2[d_it][1][mode].Enu->Draw("HISTSAME");
      AllTheHists2[d_it][2][mode].Enu->Draw("HISTSAME");
      AllTheHists2[d_it][3][mode].Enu->Draw("HISTSAME");
      AllTheHists2[d_it][4][mode].Enu->Draw("HISTSAME");
      AllTheHists2[d_it][5][mode].Enu->Draw("HISTSAME");

      RightNomCloneEnu->Divide(LeftNomCloneEnu);
      RightNomCloneEnu->SetLineColor(kBlack);
      RightNomCloneEnu->SetLineStyle(2);
      RightNomCloneEnu->Draw("HISTSAME");

      c1.cd();

      leg.Clear();
      leg.SetBorderSize(0);
      leg.SetNColumns(4);

      leg.AddEntry(AllTheHists2[d_it][0][mode].Enu, "-3 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][1][mode].Enu, "-2 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][2][mode].Enu, "-1 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][3][mode].Enu, "Nominal", "l");
      leg.AddEntry(AllTheHists2[d_it][4][mode].Enu, "+1 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][5][mode].Enu, "+2 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][6][mode].Enu, "+3 #sigma", "l");
      leg.AddEntry(RightNomCloneEnu,
                   (title2 + "/" + title1 + "Nominal").c_str(), "l");
      leg.Draw();

      std::cout << "Drawing " << ss.str() << std::endl;

      ltx.SetTextAlign(22);
      ltx.DrawLatexNDC(0.5, 0.975, ss.str().c_str());
      ltx.DrawLatexNDC(0.25, 0.84, title1.c_str());
      ltx.DrawLatexNDC(0.75, 0.84, title2.c_str());

      c1.Print(argv[3]);

      /** Topo **/

      TH1D *LeftNomCloneTopo =
          static_cast<TH1D *>(AllTheHists1[d_it][3][mode].Topo->Clone());
      TH1D *RightNomCloneTopo =
          static_cast<TH1D *>(AllTheHists2[d_it][3][mode].Topo->Clone());

      p1 = new TPad("p1", "", 0, 0.451, 0.5, 0.85);
      p1->AppendPad();
      p1->cd();

      AllTheHists1[d_it][0][mode].Topo->SetLineColor(kRed);
      AllTheHists1[d_it][1][mode].Topo->SetLineColor(kBlue);
      AllTheHists1[d_it][2][mode].Topo->SetLineColor(kMagenta);
      AllTheHists1[d_it][3][mode].Topo->SetLineColor(kBlack);
      AllTheHists1[d_it][4][mode].Topo->SetLineColor(kGreen - 2);
      AllTheHists1[d_it][5][mode].Topo->SetLineColor(kGray + 2);
      AllTheHists1[d_it][6][mode].Topo->SetLineColor(kAzure);

      AllTheHists1[d_it][6][mode].Topo->GetYaxis()->SetNdivisions(505);
      AllTheHists1[d_it][6][mode].Topo->GetYaxis()->SetTitleSize(0.05);
      AllTheHists1[d_it][6][mode].Topo->GetYaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][6][mode].Topo->GetXaxis()->SetNdivisions(505);
      AllTheHists1[d_it][6][mode].Topo->GetXaxis()->SetTitleSize(0.05);
      AllTheHists1[d_it][6][mode].Topo->GetXaxis()->SetLabelSize(0.05);

      AllTheHists1[d_it][6][mode].Topo->DrawClone("EHIST");
      AllTheHists1[d_it][0][mode].Topo->DrawClone("EHISTSAME");
      AllTheHists1[d_it][1][mode].Topo->DrawClone("EHISTSAME");
      AllTheHists1[d_it][2][mode].Topo->DrawClone("EHISTSAME");
      AllTheHists1[d_it][3][mode].Topo->DrawClone("EHISTSAME");
      AllTheHists1[d_it][4][mode].Topo->DrawClone("EHISTSAME");
      AllTheHists1[d_it][5][mode].Topo->DrawClone("EHISTSAME");

      c1.cd();

      p1rat = new TPad("p1rat", "", 0, 0, 0.5, 0.449);
      p1rat->AppendPad();
      p1rat->cd();

      AllTheHists1[d_it][6][mode].Topo->Divide(
          AllTheHists1[d_it][3][mode].Topo);
      AllTheHists1[d_it][0][mode].Topo->Divide(
          AllTheHists1[d_it][3][mode].Topo);
      AllTheHists1[d_it][1][mode].Topo->Divide(
          AllTheHists1[d_it][3][mode].Topo);
      AllTheHists1[d_it][2][mode].Topo->Divide(
          AllTheHists1[d_it][3][mode].Topo);
      AllTheHists1[d_it][4][mode].Topo->Divide(
          AllTheHists1[d_it][3][mode].Topo);
      AllTheHists1[d_it][5][mode].Topo->Divide(
          AllTheHists1[d_it][3][mode].Topo);
      AllTheHists1[d_it][3][mode].Topo->Divide(
          AllTheHists1[d_it][3][mode].Topo);

      max = -std::numeric_limits<double>::max();
      for (size_t hit = 0; hit < AllTheHists1[d_it].size(); ++hit) {
        max = std::max(max, AllTheHists1[d_it][hit][mode].Topo->GetMaximum());
      }
      AllTheHists1[d_it][6][mode].Topo->GetYaxis()->SetRangeUser(0, max * 1.1);

      AllTheHists1[d_it][6][mode].Topo->Draw("HIST");
      AllTheHists1[d_it][0][mode].Topo->Draw("HISTSAME");
      AllTheHists1[d_it][1][mode].Topo->Draw("HISTSAME");
      AllTheHists1[d_it][2][mode].Topo->Draw("HISTSAME");
      AllTheHists1[d_it][3][mode].Topo->Draw("HISTSAME");
      AllTheHists1[d_it][4][mode].Topo->Draw("HISTSAME");
      AllTheHists1[d_it][5][mode].Topo->Draw("HISTSAME");

      c1.cd();

      p2 = new TPad("p2", "", 0.5, 0.451, 1, 0.85);
      p2->AppendPad();
      p2->cd();

      AllTheHists2[d_it][0][mode].Topo->SetLineColor(kRed);
      AllTheHists2[d_it][1][mode].Topo->SetLineColor(kBlue);
      AllTheHists2[d_it][2][mode].Topo->SetLineColor(kMagenta);
      AllTheHists2[d_it][3][mode].Topo->SetLineColor(kBlack);
      AllTheHists2[d_it][4][mode].Topo->SetLineColor(kGreen - 2);
      AllTheHists2[d_it][5][mode].Topo->SetLineColor(kGray + 2);
      AllTheHists2[d_it][6][mode].Topo->SetLineColor(kAzure);

      AllTheHists2[d_it][6][mode].Topo->GetYaxis()->SetNdivisions(505);
      AllTheHists2[d_it][6][mode].Topo->GetYaxis()->SetTitleSize(0.05);
      AllTheHists2[d_it][6][mode].Topo->GetYaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][6][mode].Topo->GetXaxis()->SetNdivisions(505);
      AllTheHists2[d_it][6][mode].Topo->GetXaxis()->SetTitleSize(0.05);
      AllTheHists2[d_it][6][mode].Topo->GetXaxis()->SetLabelSize(0.05);

      AllTheHists2[d_it][6][mode].Topo->DrawClone("EHIST");
      AllTheHists2[d_it][0][mode].Topo->DrawClone("EHISTSAME");
      AllTheHists2[d_it][1][mode].Topo->DrawClone("EHISTSAME");
      AllTheHists2[d_it][2][mode].Topo->DrawClone("EHISTSAME");
      AllTheHists2[d_it][3][mode].Topo->DrawClone("EHISTSAME");
      AllTheHists2[d_it][4][mode].Topo->DrawClone("EHISTSAME");
      AllTheHists2[d_it][5][mode].Topo->DrawClone("EHISTSAME");

      c1.cd();

      p2rat = new TPad("p2rat", "", 0.5, 0, 1, 0.449);
      p2rat->AppendPad();
      p2rat->cd();

      AllTheHists2[d_it][6][mode].Topo->Divide(
          AllTheHists2[d_it][3][mode].Topo);
      AllTheHists2[d_it][0][mode].Topo->Divide(
          AllTheHists2[d_it][3][mode].Topo);
      AllTheHists2[d_it][1][mode].Topo->Divide(
          AllTheHists2[d_it][3][mode].Topo);
      AllTheHists2[d_it][2][mode].Topo->Divide(
          AllTheHists2[d_it][3][mode].Topo);
      AllTheHists2[d_it][4][mode].Topo->Divide(
          AllTheHists2[d_it][3][mode].Topo);
      AllTheHists2[d_it][5][mode].Topo->Divide(
          AllTheHists2[d_it][3][mode].Topo);
      AllTheHists2[d_it][3][mode].Topo->Divide(
          AllTheHists2[d_it][3][mode].Topo);

      max = -std::numeric_limits<double>::max();
      for (size_t hit = 0; hit < AllTheHists2[d_it].size(); ++hit) {
        max = std::max(max, AllTheHists2[d_it][hit][mode].Topo->GetMaximum());
      }
      AllTheHists2[d_it][6][mode].Topo->GetYaxis()->SetRangeUser(0, max * 1.1);

      AllTheHists2[d_it][6][mode].Topo->Draw("HIST");
      AllTheHists2[d_it][0][mode].Topo->Draw("HISTSAME");
      AllTheHists2[d_it][1][mode].Topo->Draw("HISTSAME");
      AllTheHists2[d_it][2][mode].Topo->Draw("HISTSAME");
      AllTheHists2[d_it][3][mode].Topo->Draw("HISTSAME");
      AllTheHists2[d_it][4][mode].Topo->Draw("HISTSAME");
      AllTheHists2[d_it][5][mode].Topo->Draw("HISTSAME");

      RightNomCloneTopo->Divide(LeftNomCloneTopo);
      RightNomCloneTopo->SetLineColor(kBlack);
      RightNomCloneTopo->SetLineStyle(2);
      RightNomCloneTopo->Draw("HISTSAME");

      c1.cd();

      leg.Clear();
      leg.SetBorderSize(0);
      leg.SetNColumns(4);

      leg.AddEntry(AllTheHists2[d_it][0][mode].Topo, "-3 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][1][mode].Topo, "-2 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][2][mode].Topo, "-1 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][3][mode].Topo, "Nominal", "l");
      leg.AddEntry(AllTheHists2[d_it][4][mode].Topo, "+1 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][5][mode].Topo, "+2 #sigma", "l");
      leg.AddEntry(AllTheHists2[d_it][6][mode].Topo, "+3 #sigma", "l");
      leg.AddEntry(RightNomCloneTopo,
                   (title2 + "/" + title1 + "Nominal").c_str(), "l");
      leg.Draw();

      std::cout << "Drawing " << ss.str() << std::endl;

      ltx.SetTextAlign(22);
      ltx.DrawLatexNDC(0.5, 0.975, ss.str().c_str());
      ltx.DrawLatexNDC(0.25, 0.84, title1.c_str());
      ltx.DrawLatexNDC(0.75, 0.84, title2.c_str());

      c1.Print(argv[3]);

      /*** pmuthetamu ***/
      c1.cd();
      c1.Clear();

      AllTheHists1[d_it][4][mode].pthetalep->GetYaxis()->SetNdivisions(505);
      AllTheHists1[d_it][4][mode].pthetalep->GetYaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][4][mode].pthetalep->GetYaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][4][mode].pthetalep->GetXaxis()->SetNdivisions(505);
      AllTheHists1[d_it][4][mode].pthetalep->GetXaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][4][mode].pthetalep->GetXaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][4][mode].pthetalep->GetZaxis()->SetNdivisions(505);
      AllTheHists1[d_it][4][mode].pthetalep->GetZaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][4][mode].pthetalep->GetZaxis()->SetLabelSize(0.05);

      AllTheHists1[d_it][2][mode].pthetalep->GetYaxis()->SetNdivisions(505);
      AllTheHists1[d_it][2][mode].pthetalep->GetYaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][2][mode].pthetalep->GetYaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][2][mode].pthetalep->GetXaxis()->SetNdivisions(505);
      AllTheHists1[d_it][2][mode].pthetalep->GetXaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][2][mode].pthetalep->GetXaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][2][mode].pthetalep->GetZaxis()->SetNdivisions(505);
      AllTheHists1[d_it][2][mode].pthetalep->GetZaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][2][mode].pthetalep->GetZaxis()->SetLabelSize(0.05);

      AllTheHists2[d_it][4][mode].pthetalep->GetYaxis()->SetNdivisions(505);
      AllTheHists2[d_it][4][mode].pthetalep->GetYaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][4][mode].pthetalep->GetYaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][4][mode].pthetalep->GetXaxis()->SetNdivisions(505);
      AllTheHists2[d_it][4][mode].pthetalep->GetXaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][4][mode].pthetalep->GetXaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][4][mode].pthetalep->GetZaxis()->SetNdivisions(505);
      AllTheHists2[d_it][4][mode].pthetalep->GetZaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][4][mode].pthetalep->GetZaxis()->SetLabelSize(0.05);

      AllTheHists2[d_it][2][mode].pthetalep->GetYaxis()->SetNdivisions(505);
      AllTheHists2[d_it][2][mode].pthetalep->GetYaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][2][mode].pthetalep->GetYaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][2][mode].pthetalep->GetXaxis()->SetNdivisions(505);
      AllTheHists2[d_it][2][mode].pthetalep->GetXaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][2][mode].pthetalep->GetXaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][2][mode].pthetalep->GetZaxis()->SetNdivisions(505);
      AllTheHists2[d_it][2][mode].pthetalep->GetZaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][2][mode].pthetalep->GetZaxis()->SetLabelSize(0.05);

      AllTheHists1[d_it][4][mode].pthetalep->Divide(
          AllTheHists1[d_it][3][mode].pthetalep);
      AllTheHists2[d_it][4][mode].pthetalep->Divide(
          AllTheHists2[d_it][3][mode].pthetalep);

      double maxp =
          std::max(AllTheHists1[d_it][4][mode].pthetalep->GetMaximum(),
                   AllTheHists2[d_it][4][mode].pthetalep->GetMaximum());
      AllTheHists1[d_it][4][mode].pthetalep->GetZaxis()->SetRangeUser(
          0, maxp * 1.1);
      AllTheHists2[d_it][4][mode].pthetalep->GetZaxis()->SetRangeUser(
          0, maxp * 1.1);

      AllTheHists1[d_it][2][mode].pthetalep->Divide(
          AllTheHists1[d_it][3][mode].pthetalep);
      AllTheHists2[d_it][2][mode].pthetalep->Divide(
          AllTheHists2[d_it][3][mode].pthetalep);

      double maxl =
          std::max(AllTheHists1[d_it][2][mode].pthetalep->GetMaximum(),
                   AllTheHists2[d_it][2][mode].pthetalep->GetMaximum());
      AllTheHists1[d_it][2][mode].pthetalep->GetZaxis()->SetRangeUser(
          0, maxl * 1.1);
      AllTheHists2[d_it][2][mode].pthetalep->GetZaxis()->SetRangeUser(
          0, maxl * 1.1);

      p1 = new TPad("p1", "", 0, 0.451, 0.5, 0.85);
      p1->AppendPad();
      p1->cd();
      p1->SetLeftMargin(0.15);
      p1->SetRightMargin(0.2);
      p1->SetBottomMargin(0.2);

      AllTheHists1[d_it][4][mode].pthetalep->Draw("COLZ");
      c1.cd();

      p1rat = new TPad("p1rat", "", 0, 0, 0.5, 0.449);
      p1rat->AppendPad();
      p1rat->cd();
      p1rat->SetLeftMargin(0.15);
      p1rat->SetRightMargin(0.2);
      p1rat->SetBottomMargin(0.2);
      AllTheHists1[d_it][2][mode].pthetalep->Draw("COLZ");
      c1.cd();

      p2 = new TPad("p2", "", 0.5, 0.451, 1, 0.85);
      p2->AppendPad();
      p2->cd();
      p2->SetLeftMargin(0.15);
      p2->SetRightMargin(0.2);
      p2->SetBottomMargin(0.2);
      AllTheHists2[d_it][4][mode].pthetalep->Draw("COLZ");
      c1.cd();

      p2rat = new TPad("p2rat", "", 0.5, 0, 1, 0.449);
      p2rat->AppendPad();
      p2rat->cd();
      p2rat->SetLeftMargin(0.15);
      p2rat->SetRightMargin(0.2);
      p2rat->SetBottomMargin(0.2);
      AllTheHists2[d_it][2][mode].pthetalep->Draw("COLZ");
      c1.cd();

      ltx.DrawLatexNDC(0.5, 0.975, ss.str().c_str());

      ltx.DrawLatexNDC(0.25, 0.84, (title1 + " +1 #sigma").c_str());
      ltx.DrawLatexNDC(0.75, 0.84, (title2 + " +1 #sigma").c_str());

      ltx.DrawLatexNDC(0.25, 0.44, (title1 + " -1 #sigma").c_str());
      ltx.DrawLatexNDC(0.75, 0.44, (title2 + " -1 #sigma").c_str());

      c1.Print(argv[3]);

      /*** pmuthetamu ***/
      c1.cd();
      c1.Clear();

      AllTheHists1[d_it][4][mode].pthetahmfspi->GetYaxis()->SetNdivisions(505);
      AllTheHists1[d_it][4][mode].pthetahmfspi->GetYaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][4][mode].pthetahmfspi->GetYaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][4][mode].pthetahmfspi->GetXaxis()->SetNdivisions(505);
      AllTheHists1[d_it][4][mode].pthetahmfspi->GetXaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][4][mode].pthetahmfspi->GetXaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][4][mode].pthetahmfspi->GetZaxis()->SetNdivisions(505);
      AllTheHists1[d_it][4][mode].pthetahmfspi->GetZaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][4][mode].pthetahmfspi->GetZaxis()->SetLabelSize(0.05);

      AllTheHists1[d_it][2][mode].pthetahmfspi->GetYaxis()->SetNdivisions(505);
      AllTheHists1[d_it][2][mode].pthetahmfspi->GetYaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][2][mode].pthetahmfspi->GetYaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][2][mode].pthetahmfspi->GetXaxis()->SetNdivisions(505);
      AllTheHists1[d_it][2][mode].pthetahmfspi->GetXaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][2][mode].pthetahmfspi->GetXaxis()->SetLabelSize(0.05);
      AllTheHists1[d_it][2][mode].pthetahmfspi->GetZaxis()->SetNdivisions(505);
      AllTheHists1[d_it][2][mode].pthetahmfspi->GetZaxis()->SetTitleSize(0.075);
      AllTheHists1[d_it][2][mode].pthetahmfspi->GetZaxis()->SetLabelSize(0.05);

      AllTheHists2[d_it][4][mode].pthetahmfspi->GetYaxis()->SetNdivisions(505);
      AllTheHists2[d_it][4][mode].pthetahmfspi->GetYaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][4][mode].pthetahmfspi->GetYaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][4][mode].pthetahmfspi->GetXaxis()->SetNdivisions(505);
      AllTheHists2[d_it][4][mode].pthetahmfspi->GetXaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][4][mode].pthetahmfspi->GetXaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][4][mode].pthetahmfspi->GetZaxis()->SetNdivisions(505);
      AllTheHists2[d_it][4][mode].pthetahmfspi->GetZaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][4][mode].pthetahmfspi->GetZaxis()->SetLabelSize(0.05);

      AllTheHists2[d_it][2][mode].pthetahmfspi->GetYaxis()->SetNdivisions(505);
      AllTheHists2[d_it][2][mode].pthetahmfspi->GetYaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][2][mode].pthetahmfspi->GetYaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][2][mode].pthetahmfspi->GetXaxis()->SetNdivisions(505);
      AllTheHists2[d_it][2][mode].pthetahmfspi->GetXaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][2][mode].pthetahmfspi->GetXaxis()->SetLabelSize(0.05);
      AllTheHists2[d_it][2][mode].pthetahmfspi->GetZaxis()->SetNdivisions(505);
      AllTheHists2[d_it][2][mode].pthetahmfspi->GetZaxis()->SetTitleSize(0.075);
      AllTheHists2[d_it][2][mode].pthetahmfspi->GetZaxis()->SetLabelSize(0.05);

      AllTheHists1[d_it][4][mode].pthetahmfspi->Divide(
          AllTheHists1[d_it][3][mode].pthetahmfspi);
      AllTheHists2[d_it][4][mode].pthetahmfspi->Divide(
          AllTheHists2[d_it][3][mode].pthetahmfspi);

      maxp = std::max(AllTheHists1[d_it][4][mode].pthetahmfspi->GetMaximum(),
                      AllTheHists2[d_it][4][mode].pthetahmfspi->GetMaximum());
      AllTheHists1[d_it][4][mode].pthetahmfspi->GetZaxis()->SetRangeUser(
          0, maxp * 1.1);
      AllTheHists2[d_it][4][mode].pthetahmfspi->GetZaxis()->SetRangeUser(
          0, maxp * 1.1);

      AllTheHists1[d_it][2][mode].pthetahmfspi->Divide(
          AllTheHists1[d_it][3][mode].pthetahmfspi);
      AllTheHists2[d_it][2][mode].pthetahmfspi->Divide(
          AllTheHists2[d_it][3][mode].pthetahmfspi);

      maxl = std::max(AllTheHists1[d_it][2][mode].pthetahmfspi->GetMaximum(),
                      AllTheHists2[d_it][2][mode].pthetahmfspi->GetMaximum());
      AllTheHists1[d_it][2][mode].pthetahmfspi->GetZaxis()->SetRangeUser(
          0, maxl * 1.1);
      AllTheHists2[d_it][2][mode].pthetahmfspi->GetZaxis()->SetRangeUser(
          0, maxl * 1.1);

      p1 = new TPad("p1", "", 0, 0.451, 0.5, 0.85);
      p1->AppendPad();
      p1->cd();
      p1->SetLeftMargin(0.15);
      p1->SetRightMargin(0.2);
      p1->SetBottomMargin(0.2);

      AllTheHists1[d_it][4][mode].pthetahmfspi->Draw("COLZ");
      c1.cd();

      p1rat = new TPad("p1rat", "", 0, 0, 0.5, 0.449);
      p1rat->AppendPad();
      p1rat->cd();
      p1rat->SetLeftMargin(0.15);
      p1rat->SetRightMargin(0.2);
      p1rat->SetBottomMargin(0.2);
      AllTheHists1[d_it][2][mode].pthetahmfspi->Draw("COLZ");
      c1.cd();

      p2 = new TPad("p2", "", 0.5, 0.451, 1, 0.85);
      p2->AppendPad();
      p2->cd();
      p2->SetLeftMargin(0.15);
      p2->SetRightMargin(0.2);
      p2->SetBottomMargin(0.2);
      AllTheHists2[d_it][4][mode].pthetahmfspi->Draw("COLZ");
      c1.cd();

      p2rat = new TPad("p2rat", "", 0.5, 0, 1, 0.449);
      p2rat->AppendPad();
      p2rat->cd();
      p2rat->SetLeftMargin(0.15);
      p2rat->SetRightMargin(0.2);
      p2rat->SetBottomMargin(0.2);
      AllTheHists2[d_it][2][mode].pthetahmfspi->Draw("COLZ");
      c1.cd();

      ltx.DrawLatexNDC(0.5, 0.975, ss.str().c_str());

      ltx.DrawLatexNDC(0.25, 0.84, (title1 + " +1 #sigma").c_str());
      ltx.DrawLatexNDC(0.75, 0.84, (title2 + " +1 #sigma").c_str());

      ltx.DrawLatexNDC(0.25, 0.44, (title1 + " -1 #sigma").c_str());
      ltx.DrawLatexNDC(0.75, 0.44, (title2 + " -1 #sigma").c_str());

      c1.Print(argv[3]);

      c1.Clear();
    }
  }
  c1.Print((std::string(argv[3]) + "]").c_str());
  // Can't be bothered letting root shut the files nicely
  std::cout << "[INFO]: Aborting otherwise ROOT carefully deletes thousands of "
               "histograms rather than just exiting and giving the VM back."
            << std::endl;
  abort();
}
