#include "neutvtx.h"

#include <iostream>


NeutVtx::NeutVtx(Int_t nvtx)
{
  fPos = NULL;
  SetNvtx(nvtx);

}

void
NeutVtx::SetNvtx(Int_t nvtx)
{
  fNvtx = nvtx;

  if (nvtx > 0){
	if (fPos == NULL){
	  fPos = new TObjArray(fNvtx);
	  fPos->SetOwner(kTRUE);
	}else{
	  fPos->Expand(fNvtx);
	}
  }else{
	delete fPos;
	fPos = NULL;
	fNvtx = 0;
  }

}

NeutVtx::~NeutVtx()
{
  SetNvtx(0);
}

void
NeutVtx::SetPos(int vtxid, TLorentzVector pos)
{

  TLorentzVector *pos_p;

  if (fNvtx<=vtxid){
  	SetNvtx(vtxid);
  }

  pos_p = new TLorentzVector();

  *pos_p = pos;

  fPos->AddAt(pos_p,vtxid);

}

void
NeutVtx::SetPos(Int_t nvtx, TLorentzVector *pos_array)
{

  Int_t i;

  for (i = 0 ; i < nvtx ; i++){
	SetPos(i,pos_array[i]);
  }

};

void
NeutVtx::Dump()
{
  //  if (neutcard_.quiet) return;

  int i;
  std::cout << "Event #" << EventNo << "\n";
  for (i = 0 ; i < Nvtx() ; i++){
    std::cout << "i=" << i << "\n";

    std::cout << "Vertex Pos(1)=(" << (Pos(i))->X() << ","
	      << (Pos(i))->Y() << ","
	      << (Pos(i))->Z() << ","
	      << (Pos(i))->T()  << ")"
	      << "\n";
  }
}


ClassImp(NeutVtx)
