#include "cardnameC.h"

#include <string>
#include <cstring>

extern "C" {

void resolvecardname_();

int cardistoml_(){
  if(!necardname_.isset){
    resolvecardname_();
  }
  char cstr[1025];
  cstr[1024] = '\0';

  std::strncpy(cstr, necardname_.fcard, 1024);

  return(std::string(cstr).find(".toml") != std::string::npos);
}

}