#include "toml/toml_helper.h"

#include "cardnameC.h"
#include "nrcardC.h"

#include <cstring>
#include <iostream>
#include <map>
#include <string>

extern "C" {
void nrcardtoml_();
void resolvecardname_();
}

inline std::string get_card_name() {
  if (!necardname_.isset) {
    resolvecardname_();
  }
  char card_name_cstr[1025];

  int card_name_len;
  for (int i = 1023; i >= 0; --i) {
    if (necardname_.fcard[i] != ' ') {
      card_name_len = i;
      break;
    }
  }

  std::strncpy(card_name_cstr, necardname_.fcard, card_name_len + 1);
  card_name_cstr[card_name_len + 1] = '\0';
  return card_name_cstr;
}

void nrcardtoml_() {

  std::string card_name = get_card_name();

  auto const &card_toml = toml_h::parse_card(card_name);

  auto const &neut_table = toml_h::find(card_toml, "neut");
  if (toml_h::contains_rec(neut_table, {"target", "hadron_rescattering"})) {

    auto const &hadron_rescattering_table =
        toml_h::find_rec(neut_table, {"target", "hadron_rescattering"});

    toml_h::set_option_if_present(nucres_.nucrescat, hadron_rescattering_table,
                                  "nucleon_rescattering",
                                  std::map<std::string, int>{
                                      {"on", 1},
                                      {"off", 0},
                                  });
    toml_h::set_if_present_rec(nucres_.xnucfact, hadron_rescattering_table,
                               {"nucleon_cross_sections", "mean_free_path"});

    toml_h::set_option_if_present(nucres_.nucresflg, hadron_rescattering_table,
                                  "nucleon_rescattering_nuclear_model",
                                  std::map<std::string, int>{
                                      {"lfg", 2},
                                      {"gfg", 1},
                                  });
  }
}