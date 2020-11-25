#include "spline_definitions.h"

#include "TFile.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TH3D.h"

struct EMuHistPair {
  TH2D *mu;
  TH2D *e;

  TH2D *&get(bool is_ele) { return is_ele ? e : mu; }
};

EMuHistPair Prob_HABkwd[7];
EMuHistPair Prob_Fwd[7];
EMuHistPair ProbInt_HABkwd[6];
EMuHistPair ProbInt_Fwd[6];

EMuHistPair ProbEP_HABkwd[7][24];
EMuHistPair ProbEP_Fwd[7][24];

EMuHistPair ProbEPInt_HABkwd[6][24];
EMuHistPair ProbEPInt_Fwd[6][24];

// Script to convert Ko's original histogram inputs to a more wieldy format
void convertinputs() {
  TFile *MCfile = new TFile("RadnumuCCQE.root");
  TFile *MCfileE = new TFile("RadnueCCQE.root");

  std::stringstream ss("");

  for (int is_ele = 0; is_ele < 2; ++is_ele) {
    TFile *MCF = (is_ele ? MCfileE : MCfile);
    for (int i = 0; i < 7; i++) {
      ss.str("");
      ss << "Prob" << (is_ele ? "E" : "M") << i;
      MCF->GetObject(ss.str().c_str(), Prob_HABkwd[i].get(is_ele));

      ss.str("");
      ss << "ProbA" << (is_ele ? "E" : "M") << i;
      MCF->GetObject(ss.str().c_str(), Prob_Fwd[i].get(is_ele));

      if (i < 6) {
        ss.str("");
        ss << "ProbInt" << (is_ele ? "E" : "M") << i;
        MCF->GetObject(ss.str().c_str(), ProbInt_HABkwd[i].get(is_ele));

        ss.str("");
        ss << "ProbIntA" << (is_ele ? "E" : "M") << i;
        MCF->GetObject(ss.str().c_str(), ProbInt_Fwd[i].get(is_ele));
      }

      for (int j = 0; j < 24; j++) {
        ss.str("");
        ss << "EP" << (is_ele ? "E" : "M") << i << "_" << (j + 1);
        MCF->GetObject(ss.str().c_str(), ProbEP_HABkwd[i][j].get(is_ele));

        ss.str("");
        ss << "EPA" << (is_ele ? "E" : "M") << i << "_" << (j + 1);
        MCF->GetObject(ss.str().c_str(), ProbEP_Fwd[i][j].get(is_ele));

        if (i < 6) {
          ss.str("");
          ss << "EPInt" << (is_ele ? "E" : "M") << i << "_" << (j + 1);
          MCF->GetObject(ss.str().c_str(), ProbEPInt_HABkwd[i][j].get(is_ele));

          ss.str("");
          ss << "EPIntA" << (is_ele ? "E" : "M") << i << "_" << (j + 1);
          MCF->GetObject(ss.str().c_str(), ProbEPInt_Fwd[i][j].get(is_ele));
        }
      }
    }
  }

  // The inputs that have been loaded above are a series of TH2Ds covering
  // costheta_lep and E_nu space. The following code stitches their content
  // together into a set of TGraph2Ds

  // The first Enu spline bin for the first set of Enu histograms
  Int_t INTEnuA[5] = {0, 13, 18, 24, 30};
  // The last Enu spline bin for the first set of Enu histograms
  Int_t INTEnuB[5] = {12, 17, 22, 28, 36};
  // The index of the relevant histogram in the first histogram list
  Int_t INTHist[5] = {0, 1, 2, 4, 6};
  // The first Enu spline bin for the second set of Enu histograms
  Int_t INTEnuIntA[6] = {12, 17, 22, 23, 28, 29};
  // The last Enu spline bin for the second set of Enu histograms
  Int_t INTEnuIntB[6] = {13, 18, 23, 24, 29, 30};
  // The index of the relevant histogram in the second histogram list
  Int_t INTHistInt[6] = {0, 1, 2, 3, 4, 5};

  // This list specifies which histograms to read in order to get a single
  // contiguous enu axis. The first index says whether to read from the first
  // list (INTHist) or the second list (INTHistInt) The second index is the
  std::vector<std::pair<int, int> > EnuHists = {{0, 0}, {1, 0}, {0, 1}, {1, 1},
                                                {0, 2}, {1, 2}, {1, 3}, {0, 4},
                                                {1, 4}, {1, 5}, {0, 6}};

  TFile *fout = TFile::Open("radcorr_inputs.root", "RECREATE");

  for (int is_ele = 0; is_ele < 2; ++is_ele) {

    TGraph2D colinear_gamma_total_prob;
    std::array<TGraph2D, 24> colinear_Egamma_prob;
    double lastbc;
    for (auto &h : EnuHists) {
      TH2 *hist = (h.first ? ProbInt_HABkwd[h.second].get(is_ele)
                           : Prob_HABkwd[h.second].get(is_ele));
      TH2 *hist_fwd = (h.first ? ProbInt_Fwd[h.second].get(is_ele)
                               : Prob_Fwd[h.second].get(is_ele));
      for (int i = 0; i < hist->GetXaxis()->GetNbins(); ++i) {

        // Skip if this bin has content that has already been put into a spline
        // point
        if (colinear_gamma_total_prob.GetN() &&
            (std::fabs(hist->GetXaxis()->GetBinCenter(i + 1) - lastbc) <
             1E-4)) {
          continue;
        }

        lastbc = hist->GetXaxis()->GetBinCenter(i + 1);

        bool known = false;
        for (auto bc : enu_knots) {
          if (std::fabs(bc - lastbc) < 1E-4) {
            known = true;
          }
        }
        if (!known) {
          std::cout
              << "[WARN]: When translating to spline points found unexpected "
                 "bin center: "
              << lastbc << std::endl;
        }

        for (int j = 0; j < hist->GetYaxis()->GetNbins(); ++j) {
          if (!j) {
            // std::cout << hist->GetXaxis()->GetBinCenter(i + 1) << " "
            //           << hist->GetYaxis()->GetBinCenter(j + 1) << " "
            //           << hist->GetBinContent(i + 1, j + 1) << std::endl;
          }
          colinear_gamma_total_prob.SetPoint(
              colinear_gamma_total_prob.GetN(),
              hist->GetXaxis()->GetBinCenter(i + 1),
              hist->GetYaxis()->GetBinCenter(j + 1),
              hist->GetBinContent(i + 1, j + 1));

          for (int e = 0; e < 24; ++e) {
            TH2 *hist_e = (h.first ? ProbEPInt_HABkwd[h.second][e].get(is_ele)
                                   : ProbEP_HABkwd[h.second][e].get(is_ele));

            colinear_Egamma_prob[e].SetPoint(
                colinear_Egamma_prob[e].GetN(),
                hist_e->GetXaxis()->GetBinCenter(i + 1),
                hist_e->GetYaxis()->GetBinCenter(j + 1),
                hist_e->GetBinContent(i + 1, j + 1));
          }

          // std::cout << hist->GetXaxis()->GetBinCenter(i + 1) << ", "
          //           << hist->GetYaxis()->GetBinCenter(j + 1) << ": "
          //           << hist->GetBinContent(i + 1, j + 1) << std::endl;
        }
        for (int j = 0; j < hist_fwd->GetYaxis()->GetNbins(); ++j) {
          colinear_gamma_total_prob.SetPoint(
              colinear_gamma_total_prob.GetN(),
              hist_fwd->GetXaxis()->GetBinCenter(i + 1),
              hist_fwd->GetYaxis()->GetBinCenter(j + 1),
              hist_fwd->GetBinContent(i + 1, j + 1));
          // std::cout << hist->GetXaxis()->GetBinCenter(i + 1) << ", "
          //           << hist_fwd->GetYaxis()->GetBinCenter(j + 1) << ": "
          //           << hist_fwd->GetBinContent(i + 1, j + 1) << std::endl;

          for (int e = 0; e < 24; ++e) {
            TH2 *hist_fwd_e = (h.first ? ProbEPInt_Fwd[h.second][e].get(is_ele)
                                       : ProbEP_Fwd[h.second][e].get(is_ele));

            colinear_Egamma_prob[e].SetPoint(
                colinear_Egamma_prob[e].GetN(),
                hist_fwd_e->GetXaxis()->GetBinCenter(i + 1),
                hist_fwd_e->GetYaxis()->GetBinCenter(j + 1),
                hist_fwd_e->GetBinContent(i + 1, j + 1));
          }
        }
      }
    }
    
    std::stringstream ss("");
    ss << (is_ele ? "electron" : "muon") << "_colinear_gamma_total_prob";
    fout->WriteObject(&colinear_gamma_total_prob, ss.str().c_str());

    for (int e = 0; e < 24; ++e) {
      ss.str("");
      ss <<  (is_ele ? "electron" : "muon") << "_colinear_gamma_prob_" << egamma_knots[e] << "_GeV";
      fout->WriteObject(&colinear_Egamma_prob[e], ss.str().c_str());
    }

  }
  fout->Write();
  fout->Close();
}