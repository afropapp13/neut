#include "CommonBlockIFace.h"

#include "NReWeightFactory.h"
#include "NSyst.h"

int main(int argc, char const *argv[]) {

  if (argc < 2) {
    std::cout << "[ERROR]: Need to be passed a NEUT card." << std::endl;
    return 1;
  }

  neut::CommonBlockIFace::Initialize(argv[1]);

  auto NRW = neut::rew::MakeNReWeightInstance();

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
      if(doc.size()){
        std::cout << "\t\t" << doc << std::endl;
      }
    }
  }
}