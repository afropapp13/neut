#include "toml/toml_helper.h"

#include "cardnameC.h"
#include "necardevC.h"

#include <cstring>
#include <iostream>
#include <map>
#include <string>

extern "C" {
void necardevtoml_();
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

void necardevtoml_() {

  std::string card_name = get_card_name();

  auto const &card_toml = toml_h::parse_card(card_name);

  auto const neut_table = toml_h::find(card_toml, "neut");

  nevccard_.nevtevct = toml_h::find<int>(neut_table, "number_of_events");

  auto const neut_vertex_table = toml_h::find(neut_table, "vertex");

  std::string position_type =
      toml_h::find<std::string>(neut_vertex_table, "position_type");
  nevccard_.mposevct = toml_h::ensure_valid_string_option(
      position_type,
      std::map<std::string, int>{{"origin", 1}, {"fixed", 1}, {"random", 2}});

  if (position_type == "origin") {
    std::copy_n(std::vector<float>{0.0, 0.0, 0.0}.begin(), 3,
                nevccard_.posevct);
  } else if (position_type == "fixed") {
    auto const position =
        toml_h::find<std::array<float, 3> >(neut_vertex_table, "position");
    std::copy_n(position.begin(), 3, nevccard_.posevct);
  } else if (position_type == "random") {
    nevccard_.radevct = toml_h::find<float>(neut_vertex_table, "radius");
  }

  auto const neut_probe_table = toml_h::find(neut_table, "probe");

  auto species = toml_h::find<std::string>(neut_probe_table, "species");

  nevccard_.idptevct =
      toml_h::ensure_valid_string_option(species, std::map<std::string, int>{
                                                      {"numu", 14},
                                                      {"numubar", -14},
                                                      {"nue", 12},
                                                      {"nuebar", -12},
                                                      {"nutau", 16},
                                                      {"nutaubar", -16},
                                                  });

  auto direction_type =
      toml_h::find<std::string>(neut_probe_table, "direction_type");
  nevccard_.mdirevct = toml_h::ensure_valid_string_option(
      direction_type,
      std::map<std::string, int>{{"z", 1}, {"fixed", 1}, {"random", 2}});

  if (direction_type == "z") {
    std::copy_n(std::vector<float>{0.0, 0.0, 1.0}.begin(), 3,
                nevccard_.direvct);
  } else if (direction_type == "fixed") {
    auto const direction =
        toml_h::find<std::array<float, 3> >(neut_probe_table, "direction");
    std::copy_n(direction.begin(), 3, nevccard_.direvct);
  }

  auto energy_spectrum_type =
      toml_h::find<std::string>(neut_probe_table, "energy_spectrum_type");
  nevccard_.mpvevct = toml_h::ensure_valid_string_option(
      energy_spectrum_type,
      std::map<std::string, int>{{"fixed", 1}, {"uniform", 2}, {"root", 3}});

  if (energy_spectrum_type == "fixed") {
    nevccard_.pvevct[0] = toml_h::find<float>(neut_probe_table, "energy");
  } else if (energy_spectrum_type == "uniform") {
    auto const energy_range =
        toml_h::find<std::array<float, 2> >(neut_probe_table, "energy_range");
    std::copy_n(energy_range.begin(), 2, nevccard_.pvevct);
  } else if (energy_spectrum_type == "root") {

    auto const energy_histogram = toml_h::find<std::array<std::string, 3> >(
        neut_probe_table, "energy_histogram");

    std::string filenmevct = energy_histogram[0];
    std::string histnmevct = energy_histogram[1];
    std::memset(nevccard_.filenmevct, ' ', 1000);
    std::strncpy(nevccard_.filenmevct, filenmevct.c_str(),
                 std::min(1000, int(filenmevct.size())));
    std::memset(nevccard_.histnmevct, ' ', 1000);
    std::strncpy(nevccard_.histnmevct, histnmevct.c_str(),
                 std::min(1000, int(histnmevct.size())));

    nevccard_.inmevevct = toml_h::ensure_valid_string_option(
        energy_histogram[2],
        std::map<std::string, int>{{"mev", 1}, {"gev", 0}});
  }

  /// Done reading now just writing out for the user/logs to see

  char filenmevct_cstr[1001];
  for (int i = 999; i >= 0; --i) {
    if (nevccard_.filenmevct[i] != ' ') {
      std::strncpy(filenmevct_cstr, nevccard_.filenmevct, i + 1);
      filenmevct_cstr[i + 1] = '\0';
      break;
    }
  }

  char histnmevct_cstr[1001];
  for (int i = 999; i >= 0; --i) {
    if (nevccard_.histnmevct[i] != ' ') {
      std::strncpy(histnmevct_cstr, nevccard_.histnmevct, i + 1);
      histnmevct_cstr[i + 1] = '\0';
      break;
    }
  }

  toml_h::set_option_if_present(nevccard_.outputformat, neut_table,
                                "output_format",
                                std::map<std::string, int>{
                                    {"neutroot", 1},
                                    {"NuHepMC", 2},
                                });
}