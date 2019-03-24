#include <iostream>
#include <strings.h>

#include "NeutRootHandlers.h"

void
NeutRootHandlers::nulltermstr(char *str, int len)
{
  int i;
  for (i = 0 ; i < len ; i++){
	if (str[i] == ' '){
	  str[i]='\0';
	  break;
	}
  }
}

int
NeutRootHandlers::open(char *filename, char *opt)
{
  TObject *o;
  TFile   *f;

  if (filename == NULL) return 0;

  std::cout << "Opening file " << filename << "\n";
  std::cout << "Checking whether the file is already opened or not." 
			<< "\n";

  if (!gDirectory->cd("root:/")){
	gDirectory->cd("Rint:/");
  }
  gDirectory->pwd();

  //  o = gDirectory->FindObject(filename);
  o = gROOT->GetListOfFiles()->FindObject(filename);
  if (o != NULL){
	o->Print();
  }

  if (o != NULL){
	if (o->InheritsFrom("TFile")){
	  std::cout << "File " << filename << " is already opened with option = "
				<< o->GetOption() << "\n";
	  return 0;
	}else{
	  std::cout << "Keyword " << filename 
				<< " has been used but it is not TFile.\n";
	  std::cout << "Current contents of " << filename << " are " << "\n";
	  o->Print();
	  return -1;
	}
  }

  if (opt != NULL){
	std::cout << "Opening file " << filename << "\n";
	f = new TFile(filename, opt);
  }else{
	f = new TFile(filename);
  }
  if (!(f->IsOpen())){
	std::cout << "Failed to open file " << filename << "\n";
	return -1;
  }
  return 0;
}

  
/***********************************************************************/

int
NeutRootHandlers::write(char *filename)
{
  TObject *o;

  if (filename == NULL) return 0;

  std::cout << "Writing file " << filename << "\n";
  std::cout << "Checking whether the file is already opened or not." 
			<< "\n";

  if (!gDirectory->cd("root:/")){
	gDirectory->cd("Rint:/");
  }
  gDirectory->pwd();

  //  o = gDirectory->FindObject(filename);
  o = gROOT->GetListOfFiles()->FindObject(filename);
  if (o != NULL){
	o->Print();
  }

  if (o != NULL){
	if (o->InheritsFrom("TFile")){
	  ((TFile *)o)->Write();
	  return 0;
	}else{
	  std::cout << "Keyword " << filename
				<< " has been used but it is not TFile.\n";
	  std::cout << "Current contents of " << filename << " are " << "\n";
	  o->Print();
	  return -1;
	}
  }else{ /* o == NULL */
	std::cout << "File " << filename << " has not been opened." << "\n";
	return -1;
  }

}


/***********************************************************************/

int
NeutRootHandlers::close(char *filename)
{
  TObject *o;

  if (filename == NULL) return 0;

  std::cout << "Closing file " << filename << "\n";
  std::cout << "Checking whether the file is already opened or not." 
			<< "\n";

  if (!gDirectory->cd("root:/")){
	gDirectory->cd("Rint:/");
  }
  gDirectory->pwd();

  //  o = gDirectory->FindObject(filename);
  o = gROOT->GetListOfFiles()->FindObject(filename);
  if (o != NULL){
	o->Print();
  }

  if (o != NULL){
	if (o->InheritsFrom("TFile")){
	  ((TFile *)o)->Close();
	  return 0;
	}else{
	  std::cout << "Keyword " << filename
				<< " has been used but it is not TFile.\n";
	  std::cout << "Current contents of " << filename << " are " << "\n";
	  o->Print();
	  return -1;
	}
  }else{ /* o == NULL */
	std::cout << "File " << filename << " has not been opened." << "\n";
	return -1;
  }

}

/***********************************************************************/

TTree *
NeutRootHandlers::maketree(char *filename, char *treename, char *title)
{
  TObject *o;
  TTree   *t;

  char dirname[1024];

  if (filename == NULL) return 0;

  std::cout << "Opening TTree " << treename 
			<< "in the file " << filename
			<< "." << "\n";
  std::cout << "Checking whether the file is already opened or not."
			<< "\n";

  strncpy(dirname, filename, 1024);
  strncat(dirname, ":/",  1024);
  
  if (!gDirectory->cd(dirname)){
	std::cout << "The specified file " << filename << "has not been opened."
			  << "\n";
	return NULL;
  }
  gDirectory->pwd();

  o = gDirectory->FindObject(treename);
  if (o != NULL){
	o->Print();
  }

  if (o != NULL){
	if (o->InheritsFrom("TTree")){
	  std::cout << "TTree " << treename << " is already created." << "\n";
	  return (TTree *)o;
	}else{
	  std::cout << "Keyword " << treename
				<< " has been used but it is not TTree.\n";
	  std::cout << "Current contents of " << treename << " are " << "\n";
	  o->Print();
	  return NULL;
	}
  }

  std::cout << "Creating Tree " << treename << "\n";
  t = new TTree(treename, title);
  
  if (t == NULL){
	std::cout << "Failed to create tree " << treename << "\n";
  }
  return t;
}



/***********************************************************************/

TTree *
NeutRootHandlers::attachtree(char *filename, char *treename)
{
  TObject *o;

  char dirname[1024];

  if (filename == NULL) return 0;

  std::cout << "Opening TTree " << treename 
			<< "in the file " << filename
			<< "." << "\n";
  std::cout << "Checking whether the file is already opened or not."
			<< "\n";

  strncpy(dirname, filename, 1024);
  strncat(dirname, ":/",  1024);
  
  if (!gDirectory->cd(dirname)){
	std::cout << "The specified file " << filename << "has not been opened."
			  << "\n";
	return NULL;
  }
  gDirectory->pwd();

  o = gDirectory->FindObject(treename);
  if (o != NULL){
	o->Print();
  }

  if (o != NULL){
	if (o->InheritsFrom("TTree")){
	  std::cout << "TTree " << treename << " is already created." << "\n";
	  return (TTree *)o;
	}else{
	  std::cout << "Keyword " << treename
				<< " has been used but it is not TTree.\n";
	  std::cout << "Current contents of " << treename << " are " << "\n";
	  o->Print();
	  return NULL;
	}
  }else{
	std::cout << "There is no tree named " << treename << "." << "\n";
  }
  return 0;
}
