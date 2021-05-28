#include "toml/toml_helper.h"

#include "cardnameC.h"
#include "necardC.h"
#include "neutfilepathC.h"
#include "neutmodelC.h"
#include "neutparamsC.h"
#include "nieves1p1hC.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <string>

extern "C" {
void necardtoml_();
void resolvecardname_();
void nesetfgparams_();
void mcmassgv_(int *, float *);
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

void read_target_params(toml::value const &target_table) {

  // these options are required
  neuttarget_.numatom = toml_h::find<int>(target_table, "num_nucleons");
  neuttarget_.numbndp = toml_h::find<int>(target_table, "num_bound_protons");
  neuttarget_.numfrep = toml_h::find<int>(target_table, "num_free_protons");
  neuttarget_.numbndn = neuttarget_.numatom - neuttarget_.numbndp;
  // set up nuclear model defaults

  nesetfgparams_();

  // calculate default vnuini
  float neutron_mass_gev;
  int neutron_pdg = 2112;
  mcmassgv_(&neutron_pdg, &neutron_mass_gev);
  nenupr_.vnuini = -1.0 * (std::sqrt(neutron_mass_gev * neutron_mass_gev +
                                     nenupr_.pfsurf * nenupr_.pfsurf) -
                           neutron_mass_gev);

  if (toml_h::contains(target_table, "hadron_rescattering")) {
    auto const &hadron_rescattering_table =
        toml_h::find(target_table, "hadron_rescattering");

    toml_h::set_option_if_present(
        neutcard_.nenefmodl, hadron_rescattering_table, "pion_model_low_energy",
        std::map<std::string, int>{
            {"salcedo", 0},
            {"tuned", 1},
        });

    toml_h::set_option_if_present(neutcard_.nenefmodh,
                                  hadron_rescattering_table,
                                  "pion_model_high_energy",
                                  std::map<std::string, int>{
                                      {"salcedo", 0},
                                      {"proton_neutron_separated", 1},
                                  });

    toml_h::set_option_if_present(
        neutcard_.nenefkinh, hadron_rescattering_table, "pion_kinematic_model",
        std::map<std::string, int>{
            {"isotropic_decay", 0},
            {"SAID", 1},
        });

    toml_h::set_if_present(neffpr_.fefcoul, hadron_rescattering_table,
                           "pion_coulomb_correction");

    if (toml_h::contains(hadron_rescattering_table, "pion_cross_sections")) {

      auto const &pion_cross_sections_table =
          toml_h::find(hadron_rescattering_table, "pion_cross_sections");

      toml_h::set_if_present(neffpr_.fefqe, pion_cross_sections_table, "qe");
      toml_h::set_if_present(neffpr_.fefqeh, pion_cross_sections_table,
                             "qe_high");

      toml_h::set_if_present(neffpr_.fefinel, pion_cross_sections_table,
                             "inelastic");

      toml_h::set_if_present(neffpr_.fefabs, pion_cross_sections_table,
                             "absorption");
      toml_h::set_if_present(neffpr_.fefcoh, pion_cross_sections_table,
                             "coherent");

      toml_h::set_if_present(neffpr_.fefcx, pion_cross_sections_table,
                             "charge_exchange");

      toml_h::set_if_present(neffpr_.fefcxh, pion_cross_sections_table,
                             "charge_exchange_high");

      toml_h::set_if_present(neffpr_.fefqehf, pion_cross_sections_table,
                             "qe_fraction");

      toml_h::set_if_present(neffpr_.fefcohf, pion_cross_sections_table,
                             "coh_fraction");
      toml_h::set_if_present(neffpr_.fefcxhf, pion_cross_sections_table,
                             "charge_exchange_fraction");

      toml_h::set_if_present(neffpr_.fefall, pion_cross_sections_table,
                             "mean_free_path");
    }

    toml_h::set_option_if_present(neutcard_.nenefo16, hadron_rescattering_table,
                                  "pion_rescattering",
                                  std::map<std::string, int>{
                                      {"on", 0},
                                      {"off", 1},
                                  });

    toml_h::set_option_if_present(nuceffver_.nefkinver,
                                  hadron_rescattering_table,
                                  "pion_rescattering_nuclear_model",
                                  std::map<std::string, int>{
                                      {"data_driven", 0},
                                      {"lfg", 1},
                                  });
  }

  if (toml_h::contains(target_table, ("nuclear_effects"))) {
    auto const &nuclear_effects_table =
        toml_h::find(target_table, "nuclear_effects");

    toml_h::set_option_if_present(neutcard_.nefrmflg, nuclear_effects_table,
                                  "fermi_motion",
                                  std::map<std::string, int>{
                                      {"on", 0},
                                      {"off", 1},
                                  });

    toml_h::set_option_if_present(neutcard_.nepauflg, nuclear_effects_table,
                                  "pauli_blocking",
                                  std::map<std::string, int>{
                                      {"on", 0},
                                      {"off", 1},
                                  });

    toml_h::set_option_if_present(neutpiabs_.neabspiemit, nuclear_effects_table,
                                  "nucleon_ejection",
                                  std::map<std::string, int>{
                                      {"off", 0},
                                      {"on", 1},
                                  });

    toml_h::set_option_if_present(neutcard_.nucdexite, nuclear_effects_table,
                                  "gamma_from_O16_deexcitation",
                                  std::map<std::string, int>{
                                      {"off", 0},
                                      {"on", 1},
                                  });

    toml_h::set_if_present(nenupr_.pfsurf, nuclear_effects_table,
                           "fermi_surface_momentum_gev");
    toml_h::set_if_present(nenupr_.pfmax, nuclear_effects_table,
                           "fermi_surface_energy_gev");
    toml_h::set_if_present(nenupr_.vnuini, nuclear_effects_table,
                           "initial_potential_energy_gev");
    toml_h::set_if_present(nenupr_.vnufin, nuclear_effects_table,
                           "final_potential_energy_gev");

    toml_h::set_option_if_present(nenupr_.iformlen, nuclear_effects_table,
                                  "formation_length",
                                  std::map<std::string, int>{
                                      {"on", 1},
                                      {"off", 0},
                                      {"on_except_qe", 110},
                                      {"on_multpi_DIS", 100},
                                  });
    toml_h::set_if_present(nenupr_.fzmu2, nuclear_effects_table,
                           "formation_zone_scat_parameter");
  }
}

void read_QE_params(toml::value const &QE_table) {

  int MDLQE = 0;

  std::string model_choice = toml_h::find<std::string>(QE_table, "model");
  (void)toml_h::ensure_valid_string_option(model_choice,
                                           std::vector<std::string>{
                                               {"nieves"},
                                               {"benhar"},
                                               {"effective_sf"},
                                               {"TEM"},
                                               {"llewelyn-smith"},
                                           });
  if (model_choice == "nieves") {
    MDLQE = 2000;
  } else if (model_choice == "benhar") {
    MDLQE = 400;
  } else if (model_choice == "effective_sf") {
    MDLQE = 600;
  } else if (model_choice == "TEM") {
    MDLQE = 700;
  }

  toml_h::set_if_present_rec(nemdls_.pfsf, QE_table,
                             {"behnar", "fermi_surface_for_pauli_mev"});
  toml_h::set_if_present_rec(nenupr_.sfebshift, QE_table,
                             {"behnar", "separation_energy_shift_gev"});
  toml_h::set_option_if_present_rec(
      nenupr_.sfebnegbeh, QE_table,
      {"behnar", "separation_energy_negative_behavior"},
      std::map<std::string, int>{
          {"pin_to_0", 0},
          {"redraw_from_sf", 1},
      });

  toml_h::set_if_present(nemdls_.kapp, QE_table, "kappa");

  int erpa = 0;
  toml_h::set_option_if_present(erpa, QE_table, "erpa",
                                std::map<std::string, int>{
                                    {"off", 0},
                                    {"on", 1},
                                });
  if (erpa) {
    if (model_choice == "nieves") {
      std::cout
          << "[ERROR]: neut.interaction_model.QE.erpa = \"on\" is incompatible "
             "with neut.interaction_model.QE.model = \"nieves\"."
          << std::endl;
      abort();
    }
    MDLQE += 1000;
  }

  toml_h::set_if_present(nemdls_.sccfv, QE_table,
                         "second_class_currents_vector_scale");
  toml_h::set_if_present(nemdls_.sccfa, QE_table,
                         "second_class_currents_axial_scale");
  toml_h::set_if_present(nemdls_.fpqe, QE_table, "pseudo_scalar_scale");

  toml_h::set_if_present(neutpiless_.ipilessdcy, QE_table,
                         "pionless_delta_decay");
  toml_h::set_if_present(neutpiless_.rpilessdcy, QE_table,
                         "pionless_delta_decay_fraction");

  toml_h::set_option_if_present_rec(
      nievesqepar_.nvqerfg, QE_table, {"nieves", "nuclear_model"},
      std::map<std::string, int>{{"rfg", 1}, {"lfg", 0}});

  toml_h::set_option_if_present_rec(nievesqepar_.nvqerpa, QE_table,
                                    {"nieves", "rpa"},
                                    std::map<std::string, int>{
                                        {"off", 0},
                                        {"on", 1},
                                    });

  toml_h::set_if_present_rec(nievesqepar_.xnvrpafp0in, QE_table,
                             {"nieves", "rpa_parameters", "fp0in"});
  toml_h::set_if_present_rec(nievesqepar_.xnvrpapf0ex, QE_table,
                             {"nieves", "rpa_parameters", "fp0ex"});
  toml_h::set_if_present_rec(nievesqepar_.xnvrpafstar, QE_table,
                             {"nieves", "rpa_parameters", "fstar"});
  toml_h::set_if_present_rec(nievesqepar_.xnvrpaf, QE_table,
                             {"nieves", "rpa_parameters", "af"});
  toml_h::set_if_present_rec(nievesqepar_.xnvrpapilambda, QE_table,
                             {"nieves", "rpa_parameters", "pilambda"});
  toml_h::set_if_present_rec(nievesqepar_.xnvrpacr0, QE_table,
                             {"nieves", "rpa_parameters", "cr0"});
  toml_h::set_if_present_rec(nievesqepar_.xnvrparholambda, QE_table,
                             {"nieves", "rpa_parameters", "rholambda"});
  toml_h::set_if_present_rec(nievesqepar_.xnvrpagp, QE_table,
                             {"nieves", "rpa_parameters", "gp"});
  toml_h::set_if_present_rec(nievesqepar_.xnvrpaxmpi, QE_table,
                             {"nieves", "rpa_parameters", "xmpi"});
  toml_h::set_if_present_rec(nievesqepar_.xnvrpaxmrho, QE_table,
                             {"nieves", "rpa_parameters", "xmrho"});
  toml_h::set_if_present_rec(nievesqepar_.xnvrpairel, QE_table,
                             {"nieves", "rpa_parameters", "irel"});

  toml_h::set_option_if_present(nievesqepar_.nvqebind, QE_table,
                                {"nieves", "binding_energy"},
                                std::map<std::string, int>{
                                    {"off", 0},
                                    {"on", 1},
                                });
  toml_h::set_if_present_rec(nievesqepar_.nvbindfermicor, QE_table,
                             {"nieves", "binding_energy_shift"});

  if (toml_h::contains(QE_table, "vector_form_factors")) {
    auto const &vector_form_factors_table =
        toml_h::find(QE_table, "vector_form_factors");

    std::string cc_form =
        toml_h::find<std::string>(vector_form_factors_table, "cc_form");
    (void)toml_h::ensure_valid_string_option(cc_form, std::vector<std::string>{
                                                          {"dipole"},
                                                          {"BBBA05"},
                                                          {"BBBA07"},
                                                      });

    if (model_choice == "dipole") {
      MDLQE += 1;
    } else if (model_choice == "BBBA05") {
      MDLQE += 2;
    } else if (model_choice == "BBBA07") {
      MDLQE += 3;
    }

    std::string nc_form =
        toml_h::find<std::string>(vector_form_factors_table, "nc_form");
    (void)toml_h::ensure_valid_string_option(nc_form, std::vector<std::string>{
                                                          {"default"},
                                                          {"dipole"},
                                                          {"BBBA05"},
                                                      });

    if (model_choice == "dipole") {
      MDLQE += 10;
    } else if (model_choice == "BBBA05") {
      MDLQE += 20;
    }

    toml_h::set_if_present(nemdls_.xmvqe, vector_form_factors_table,
                           "vector_mass");
  }

  nemdls_.mdlqe = MDLQE;

  if (toml_h::contains(QE_table, "axial_form_factors")) {
    auto const &axial_form_factors_table =
        toml_h::find(QE_table, "axial_form_factors");

    std::string form =
        toml_h::find<std::string>(axial_form_factors_table, "form");
    nemdls_.mdlqeaf =
        toml_h::ensure_valid_string_option(form, std::map<std::string, int>{
                                                     {"dipole", 1},
                                                     {"BBBA07", 2},
                                                     {"two_component", 3},
                                                     {"three_component", 4},
                                                     {"z_expansions", 5},
                                                 });

    toml_h::set_if_present_rec(nemdls_.xmaqe, axial_form_factors_table,
                               {"dipole", "axial_mass"});
    toml_h::set_if_present_rec(nemdls_.xmancel, axial_form_factors_table,
                               {"dipole", "nc_axial_mass"});

    toml_h::set_if_present_rec(nemdls_.axffalpha, axial_form_factors_table,
                               {"n_component", "alpha"});
    toml_h::set_if_present_rec(nemdls_.axffbeta, axial_form_factors_table,
                               {"n_component", "beta"});
    toml_h::set_if_present_rec(nemdls_.axffgamma, axial_form_factors_table,
                               {"n_component", "gamma"});
    toml_h::set_if_present_rec(nemdls_.axfftheta, axial_form_factors_table,
                               {"n_component", "theta"});

    toml_h::set_if_present_rec(nemdls_.axzexpq4, axial_form_factors_table,
                               {"z_expansion", "q4"});
    toml_h::set_if_present_rec(nemdls_.axzexpnt, axial_form_factors_table,
                               {"z_expansion", "num_terms"});
    toml_h::set_if_present_rec(nemdls_.axzexpt0, axial_form_factors_table,
                               {"z_expansion", "t0"});
    toml_h::set_if_present_rec(nemdls_.axzexptc, axial_form_factors_table,
                               {"z_expansion", "t_cut"});
    toml_h::set_if_present_rec(nemdls_.axzexpa0, axial_form_factors_table,
                               {"z_expansion", "a0"});
    toml_h::set_if_present_rec(nemdls_.axzexpa1, axial_form_factors_table,
                               {"z_expansion", "a1"});
    toml_h::set_if_present_rec(nemdls_.axzexpa2, axial_form_factors_table,
                               {"z_expansion", "a2"});
    toml_h::set_if_present_rec(nemdls_.axzexpa3, axial_form_factors_table,
                               {"z_expansion", "a3"});
    toml_h::set_if_present_rec(nemdls_.axzexpa4, axial_form_factors_table,
                               {"z_expansion", "a4"});
    toml_h::set_if_present_rec(nemdls_.axzexpa5, axial_form_factors_table,
                               {"z_expansion", "a5"});
    toml_h::set_if_present_rec(nemdls_.axzexpa6, axial_form_factors_table,
                               {"z_expansion", "a6"});
    toml_h::set_if_present_rec(nemdls_.axzexpa7, axial_form_factors_table,
                               {"z_expansion", "a7"});
    toml_h::set_if_present_rec(nemdls_.axzexpa8, axial_form_factors_table,
                               {"z_expansion", "a8"});
    toml_h::set_if_present_rec(nemdls_.axzexpa9, axial_form_factors_table,
                               {"z_expansion", "a9"});
  }
}

void read_2p2h_params(toml::value const &table_2p2h) {
  toml_h::set_option_if_present(nemdls_.mdl2p2h, table_2p2h, "model",
                                std::map<std::string, int>{
                                    {"nieves_table", 1},
                                    {"nieves_hadron_tensor", 2},
                                });
}

void read_RES_params(toml::value const &RES_table) {
  toml_h::set_option_if_present(
      nemdls_.mdlspi, RES_table, "model",
      std::map<std::string, int>{{"rein_sehgal", 1}, {"kabirnezhad", 3}});

  if (nemdls_.mdlspi == 1) { // RS
    // GetRS Table
    toml::value RES_RS_table; // if its empty, you're going to get defaults.
    if (toml_h::contains(RES_table, "rein_sehgal")) {
      RES_RS_table = toml_h::find(RES_table, "rein_sehgal");
    }

    toml_h::set_if_present(neut1pi_.rnebgscl, RES_RS_table,
                           "nonresonant_background_scale");

    toml_h::set_option_if_present(
        nemdls_.mdlspiej, RES_RS_table, "ejection_model",
        std::map<std::string, int>{{"isotropic", 0},
                                   {"delta_only", 1},
                                   {"all_resonances", 2},
                                   {"delta_plus_isotropic", 3}});

    toml_h::set_option_if_present(
        neut1pi_.neiff, RES_RS_table, "form_factors",
        std::map<std::string, int>{{"originalFF", 0}, {"graczyk_sobczyk", 1}});

    toml_h::set_option_if_present_rec(neut1pi_.nenrtype, RES_RS_table,
                                      {"originalFF", "nr_type"},
                                      std::map<std::string, int>{
                                          {"original", 0},
                                          {"halved", 1},
                                      });

    toml_h::set_if_present_rec(neut1pi_.xmarsres, RES_RS_table,
                               {"originalFF", "axial_mass"});
    toml_h::set_if_present_rec(neut1pi_.xmvrsres, RES_RS_table,
                               {"originalFF", "vector_mass"});

    toml_h::set_if_present_rec(neut1pi_.xmanffres, RES_RS_table,
                               {"graczyk_sobczyk", "axial_mass"});
    toml_h::set_if_present_rec(neut1pi_.xmvnffres, RES_RS_table,
                               {"graczyk_sobczyk", "vector_mass"});
    toml_h::set_if_present_rec(neut1pi_.rneca5i, RES_RS_table,
                               {"graczyk_sobczyk", "ca5"});

  } else if (nemdls_.mdlspi == 3) {
    // GetMK Table
    toml::value RES_MK_table; // if its empty, you're going to get defaults.
    if (toml_h::contains(RES_table, "kabirnezhad")) {
      RES_MK_table = toml_h::find(RES_table, "kabirnezhad");
    }
    toml_h::set_if_present(nemdls_.xmabkgm, RES_MK_table, "bkg_axial_mass");
    toml_h::set_if_present(neut1pi_.xmanffres, RES_MK_table, "axial_mass");
    toml_h::set_if_present(neut1pi_.xmvnffres, RES_MK_table, "vector_mass");
    toml_h::set_if_present(neut1pi_.rneca5i, RES_MK_table, "ca5");
  }
}

void read_COH_params(toml::value const &COH_table) {
  toml_h::set_option_if_present(neutcoh_.necohepi, COH_table, "model",
                                std::map<std::string, int>{
                                    {"rein_sehgal", 0},
                                    {"kartavtsev", 1},
                                    {"berger_sehgal", 2},
                                });
  // reweight only
  // toml_h::set_if_present(nemdls_.xmacoh, COH_table, "axial_mass");
  // toml_h::set_if_present(nemdls_.rad0nu, COH_table, "radius");
  // toml_h::set_if_present(nemdls_.fa1coh, COH_table, "fa1");
  // toml_h::set_if_present(nemdls_.fb1coh, COH_table, "fa2");
}

void read_DIS_params(toml::value const &DIS_table) {
  toml_h::set_option_if_present(neutdis_.nepdf, DIS_table, "pdf_set",
                                std::map<std::string, int>{
                                    {"grv94_di", 7},
                                    {"grv98_lo", 12},
                                });
  toml_h::set_option_if_present(neutdis_.nebodek, DIS_table,
                                "bodek_yang_correction",
                                std::map<std::string, int>{
                                    {"BY19", 2},
                                    {"BY05", 1},
                                    {"off", 0},
                                });

  toml_h::set_option_if_present(neutdis_.nebylw, DIS_table,
                                "bodek_yang_low_W_K_factor",
                                std::map<std::string, int>{
                                    {"on", 1},
                                    {"off", 0},
                                });

  toml_h::set_option_if_present(neutdis_.nebyforcet1, DIS_table,
                                "bodek_yang_force_ax_vec_equal",
                                std::map<std::string, int>{
                                    {"on", 1},
                                    {"off", 0},
                                });

  toml_h::set_option_if_present(neutdis_.nemult, DIS_table,
                                "multiplicity_model",
                                std::map<std::string, int>{
                                    {"neut", 0},
                                    {"fit_to_bc", 1},
                                    {"AGKY", 2},
                                });
}

void read_DIF_params(toml::value const &DIF_table) {
  toml_h::set_option_if_present(neutdif_.nedifpi, DIF_table, "model",
                                std::map<std::string, int>{
                                    {"rein_w_cut", 0},
                                    {"rein", 1},
                                });
  // reweight only
  // toml_h::set_if_present(nemdls_.xmadif, DIF_table, "axial_mass");
  // toml_h::set_if_present(nemdls_.nucvoldif, DIF_table, "nuclear_volume");
}

std::map<std::string, int> static const ModeNameMapping = {
    // from nemodsel.F on 2020/10/21
    // DATA MODNEU /  1, 11, 12, 13, 21, 31, 32, 33, 34, 41,
    // &              51, 51, 52, 16, 36, 22, 42, 43, 23, 44,
    // &              45,  0, 26, 46, 17, 38, 39, 2,  15, 35/

    {"CC_QE_nu", 1}, // mode = 1

    {"CC_RES_ppi+_nu", 2}, // mode = 11
    {"CC_RES_ppi0_nu", 3}, // mode = 12
    {"CC_RES_npi+_nu", 4}, // mode = 13

    {"CC_multi_pi_nu", 5}, // mode = 21

    {"NC_RES_npi0_nu", 6}, // mode = 31
    {"NC_RES_ppi0_nu", 7}, // mode = 32
    {"NC_RES_ppi-_nu", 8}, // mode = 33
    {"NC_RES_npi+_nu", 9}, // mode = 34

    {"NC_multi_pi_nu", 10}, // mode = 41

    {"NC_elastic_free_p_nu", 11},  // mode = 51
    {"NC_elastic_bound_p_nu", 12}, // mode = 51
    {"NC_elastic_n_nu", 13},       // mode = 52

    {"CC_COH_nu", 14}, // mode = 16

    {"NC_COH_nu", 15}, // mode = 36

    {"CC_eta_nu", 16}, // mode = 22

    {"NC_eta_n_nu", 17}, // mode = 42
    {"NC_eta_p_nu", 18}, // mode = 43

    {"CC_kaon_nu", 19}, // mode = 23

    {"NC_kaon_n_nu", 20}, // mode = 44
    {"NC_kaon_p_nu", 21}, // mode = 45

    {"CC_DIS_nu", 23}, // mode = 26

    {"NC_DIS_nu", 24}, // mode = 46

    {"CC_1gamma_nu", 25}, // mode = 17

    {"NC_1gamma_n_nu", 26}, // mode = 38
    {"NC_1gamma_p_nu", 27}, // mode = 39

    {"CC_2p2h_nu", 28}, // mode = 2

    {"CC_DIF_nu", 29}, // mode = 15

    {"NC_DIF_nu", 30}, // mode = 35

    // from nemodsel.F on 2020/10/21
    //  DATA MODNEUB/ -1,-11,-12,-13,-21,-31,-32,-33,-34,-41,
    // &              -1,-51,-51,-52,-16,-36,-22,-42,-43,-23,
    // &              -44,-45,-26,-46,-17,-38,-39,-2,-15,-35/

    {"CC_QE_free_proton_nubar", 1}, // mode = -1

    {"CC_RES_npi-_nubar", 2}, // mode = -11
    {"CC_RES_ppi0_nubar", 3}, // mode = -12
    {"CC_RES_ppi-_nubar", 4}, // mode = -13

    {"CC_multi_pi_nubar", 5}, // mode = -21

    {"NC_RES_npi0_nubar", 6}, // mode = -31
    {"NC_RES_ppi0_nubar", 7}, // mode = -32
    {"NC_RES_ppi-_nubar", 8}, // mode = -33
    {"NC_RES_npi+_nubar", 9}, // mode = -34

    {"NC_multi_pi_nubar", 10}, // mode = -41

    {"CC_QE_bound_proton_nubar", 11}, // mode = -1

    {"NC_elastic_free_p_nubar", 12},  // mode = -51
    {"NC_elastic_bound_p_nubar", 13}, // mode = -51
    {"NC_elastic_n_nubar", 14},       // mode = -52

    {"CC_COH_nubar", 15}, // mode = -16

    {"NC_COH_nubar", 16}, // mode = -36

    {"CC_eta_nubar", 17}, // mode = -22

    {"NC_eta_n_nubar", 18}, // mode = -42
    {"NC_eta_p_nubar", 19}, // mode = -43

    {"CC_kaon_nubar", 20}, // mode = -23

    {"NC_kaon_n_nubar", 21}, // mode = -44
    {"NC_kaon_p_nubar", 22}, // mode = -45

    {"CC_DIS_nubar", 23}, // mode = -26

    {"NC_DIS_nubar", 24}, // mode = -46

    {"CC_1gamma_nubar", 25}, // mode = -17

    {"NC_1gamma_n_nubar", 26}, // mode = -38
    {"NC_1gamma_p_nubar", 27}, // mode = -39

    {"CC_2p2h_nubar", 28}, // mode = -2

    {"CC_DIF_nubar", 29}, // mode = -15

    {"NC_DIF_nubar", 30}, // mode = -35

};

void necardtoml_() {

  std::string card_name = get_card_name();

  auto const &card_toml = toml_h::parse_card(card_name);

  auto const &neut_table = toml_h::find(card_toml, "neut");

  if (toml_h::contains(neut_table, "cross_section_table_path")) {
    std::string cross_section_table_path =
        toml_h::find<std::string>(neut_table, "cross_section_table_path");
    std::memset(neutfilepath_.crstblpath, ' ', 1024);
    std::strncpy(neutfilepath_.crstblpath, cross_section_table_path.c_str(),
                 std::min(1024, int(cross_section_table_path.size())));
  }

  auto const &target_table = toml_h::find(neut_table, "target");
  read_target_params(target_table);

  if (toml_h::contains(neut_table, "interaction_model")) {

    auto const &interaction_model_table =
        toml_h::find(neut_table, "interaction_model");

    std::string mode_choice_type =
        toml_h::find<std::string>(interaction_model_table, "mode_choice_type");
    (void)toml_h::ensure_valid_string_option(
        mode_choice_type,
        std::vector<std::string>{{"normal"}, {"rescaled"}, {"single"}});

    if (mode_choice_type == "normal") {
      neutcard_.nemodflg = 0;
    } else if (mode_choice_type == "rescaled") {
      neutcard_.nemodflg = -1;

      if (toml_h::contains(interaction_model_table, "mode_rescaling")) {
        auto const &mode_rescaling_table =
            toml_h::find(interaction_model_table, "mode_rescaling");

        for (auto const &mode : ModeNameMapping) {

          if (!toml_h::contains(mode_rescaling_table, mode.first)) {
            continue;
          }

          if (mode.first.find("nubar") != std::string::npos) {
            neutcard_.crsneutb[mode.second - 1] =
                toml_h::find<float>(mode_rescaling_table, mode.first);
          } else {
            neutcard_.crsneut[mode.second - 1] =
                toml_h::find<float>(mode_rescaling_table, mode.first);
          }
        }
      }

    } else {
      neutcard_.nemodflg = toml_h::ensure_valid_string_option(
          interaction_model_table, "mode", ModeNameMapping);
    }

    toml_h::set_option_if_present(
        neutradcorr_.iradcorr, interaction_model_table, "radiative_corrections",
        std::map<std::string, int>{
            {"on", 1},
            {"off", 0},
        });

    if (toml_h::contains(interaction_model_table, "QE")) {
      read_QE_params(toml_h::find(interaction_model_table, "QE"));
    }
    if (toml_h::contains(interaction_model_table, "2p2h")) {
      read_2p2h_params(toml_h::find(interaction_model_table, "2p2h"));
    }
    if (toml_h::contains(interaction_model_table, "RES")) {
      read_RES_params(toml_h::find(interaction_model_table, "RES"));
    }
    if (toml_h::contains(interaction_model_table, "COH")) {
      read_COH_params(toml_h::find(interaction_model_table, "COH"));
    }
    if (toml_h::contains(interaction_model_table, "DIS")) {
      read_DIS_params(toml_h::find(interaction_model_table, "DIS"));
    }
    if (toml_h::contains(interaction_model_table, "DIF")) {
      read_DIF_params(toml_h::find(interaction_model_table, "DIF"));
    }
  }
}