#include "CommonBlockIFace.h"

#include "NReWeightFactory.h"
#include "NSyst.h"

std::string dial_name = "";
std::string neut_card = "";
double tweak_value = 0;
bool is_tweak = true;

void Usage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << " --neut-card <neut.card> [--dial <name>] "
               "[--one-up|--one-down|--get-value-for-tweak "
               "<value>|--get-tweak-for-value <value>]"
            << std::endl;
}

void ParseOpts(int argc, char const *argv[]) {
  for (int opt_it = 1; opt_it < argc; opt_it++) {
    std::string arg = argv[opt_it];
    if ((arg == "-h") || (arg == "-?") || (arg == "--help")) {
      Usage(argv);
      exit(0);
    } else if ((opt_it + 1) < argc) {
      if (arg == "--neut-card") {
        neut_card = argv[++opt_it];
      } else if (arg == "--dial") {
        dial_name = argv[++opt_it];
      } else if (arg == "--get-value-for-tweak") {
        tweak_value = std::stod(argv[++opt_it]);
      } else if (arg == "--get-tweak-for-value") {
        tweak_value = std::stod(argv[++opt_it]);
        is_tweak = false;
      }
    } else if (arg == "--one-up") {
      tweak_value = 1;
    } else if (arg == "--one-down") {
      tweak_value = -1;
    }
  }
  if (!neut_card.size()) {
    std::cout << "[ERROR]: Need to be passed a NEUT card via --neut-card"
              << std::endl;

    Usage(argv);
    exit(1);
  }
}

int main(int argc, char const *argv[]) {

  ParseOpts(argc, argv);
  if (dial_name.size()) {
    silencelibgfortran_();
  }
  neut::CommonBlockIFace::Initialize(neut_card, bool(dial_name.size()));
  if (dial_name.size()) {
    unsilencelibgfortran_();
  }

  auto NRW = neut::rew::MakeNReWeightInstance();

  if (dial_name.size()) {
    auto dial = NRW->DialFromString(dial_name);

    if (is_tweak) {
      NRW->SetDial_NumberOfSigmas(dial, tweak_value);
      NRW->Reconfigure();
      std::cout << NRW->GetDial_To_Value(dial) << std::endl;
      return 0;
    } else {
      NRW->SetDial_To_Value(dial, tweak_value);
      std::cout << NRW->GetDial_Tweak(dial) << std::endl;
      return 0;
    }
  }

  for (auto const &d : neut::NSyst::Dials) {
    if (NRW->DialIsHandled(d.second)) {
      double FromValue = NRW->GetDial_From_Value(d.second);
      double OneSigmaDown = NRW->GetDial_OneSigma(d.second, -1);
      double OneSigmaUp = NRW->GetDial_OneSigma(d.second, 1);
      if (!OneSigmaUp && !OneSigmaDown) {
        std::cout << "\t" << d.first << " (id = " << d.second << "): "
                  << " = { From: " << FromValue << " }" << std::endl;
      } else {
        std::cout << "\t" << d.first << " (id = " << d.second << "): "
                  << " = { From: " << FromValue
                  << ", SigmaDown: " << OneSigmaDown
                  << ", SigmaUp: " << OneSigmaUp << " }" << std::endl;
      }
      std::string const &doc = neut::NSyst::GetDialDocumentation(d.second);
      if (doc.size()) {
        std::cout << "\t\t" << doc << std::endl;
      }
    }
  }
}