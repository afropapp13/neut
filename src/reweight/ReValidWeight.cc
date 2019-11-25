#include "CommonBlockIFace.h"
#include "NReWeightFactory.h"
#include "neutrootTreeSingleton.h"

#include <iostream>

using namespace neut;
using namespace neut::rew;

int main(int argc, char const *argv[]) {

  if (argc != 3) {
    std::cout
        << "[ERROR]: Expected to be passed: <NEUTCard.card> <NEUTVector.root>"
        << std::endl;
    return 1;
  }
  std::cout << CommonBlockIFace::ParamsToString();

  CommonBlockIFace::Initialize(argv[1]);
  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();

  std::cout << CommonBlockIFace::ParamsToString(1);

  NReWeight *NRW = DefaultConfiguredNReWeightFactory();

  NeutrootTreeSingleton &neutRdr = NeutrootTreeSingleton::Initialize(argv[2]);

  Long64_t Nents = neutRdr.GetEntries();

  NRW->Systematics().Init(kXSecTwkDial_MaCCQE, 1);
  NRW->Reconfigure();

  NeutVect *nvct = neutRdr.GetNeutVectAddress();
  for (Long64_t ent = 0; ent < Nents; ++ent) {
    std::cout << "Getting entry: " << ent << std::endl;
    neutRdr.GetEntry(ent);

    CommonBlockIFace::ReadVect(nvct);

    std::cout << NRW->CalcWeight() << std::endl;
  }
}
