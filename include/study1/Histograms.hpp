#pragma once

#include <TH1.h>
#include <TH2.h>
#include <ROOT/TThreadedObject.hxx>

using TTO_TH1D = ROOT::TThreadedObject<TH1D>;
using TTO_TH2D = ROOT::TThreadedObject<TH2D>;
using up_TTO_TH1D = const std::unique_ptr<TTO_TH1D>;
using up_TTO_TH2D = const std::unique_ptr<TTO_TH2D>;

namespace study1 {

struct Electron {
    up_TTO_TH1D hist1D_p = std::make_unique<TTO_TH1D>("p_e", "", 200, 0, 11);
    up_TTO_TH1D hist1D_phi = std::make_unique<TTO_TH1D>("phi_e", "", 200, -190, 190);
    up_TTO_TH1D hist1D_theta = std::make_unique<TTO_TH1D>("theta_e", "", 200, 0, 35);
    up_TTO_TH1D hist1D_chi2 = std::make_unique<TTO_TH1D>("chi2_e", "", 200, -6, 6);
    up_TTO_TH1D hist1D_vz = std::make_unique<TTO_TH1D>("vz_e", "", 200, -50, 40);

    up_TTO_TH1D hist1D_p_cut = std::make_unique<TTO_TH1D>("p_e_cut", "", 200, 0, 11);
    up_TTO_TH1D hist1D_phi_cut = std::make_unique<TTO_TH1D>("phi_e_cut", "", 200, -190, 190);
    up_TTO_TH1D hist1D_theta_cut = std::make_unique<TTO_TH1D>("theta_e_cut", "", 200, 0, 35);
    up_TTO_TH1D hist1D_chi2_cut = std::make_unique<TTO_TH1D>("chi2_e_cut", "", 200, -6, 6);
    up_TTO_TH1D hist1D_vz_cut = std::make_unique<TTO_TH1D>("vz_e_cut", "", 200, -50, 40);

    up_TTO_TH1D hist1D_Q2 = std::make_unique<TTO_TH1D>("Q2", "", 200, 0, 1);
    up_TTO_TH1D hist1D_W = std::make_unique<TTO_TH1D>("W", "", 200, 0, 5.5);
    up_TTO_TH1D hist1D_xB = std::make_unique<TTO_TH1D>("xB", "", 200, 0, 1.2);
    up_TTO_TH1D hist1D_MM_eP = std::make_unique<TTO_TH1D>("MM_eP", "", 200, -5, 50); 

    up_TTO_TH1D hist1D_pp_calc = std::make_unique<TTO_TH1D>("pp_calc", "", 200, 0, 0.4);

    up_TTO_TH1D hist1D_delta_phi = std::make_unique<TTO_TH1D>("delta_phi", "", 200, -180, 180);
    up_TTO_TH1D hist1D_delta_phi_ATOF = std::make_unique<TTO_TH1D>("delta_phi_ATOF", "", 200, -360, 360);
    up_TTO_TH1D hist1D_delta_phi_cut = std::make_unique<TTO_TH1D>("delta_phi_cut", "", 200, -360, 360);

    up_TTO_TH2D hist2D_delta_phi_vs_phiElec= std::make_unique<TTO_TH2D>("delta_phi_vs_phiElec", "", 200, -360, 360, 200, -360, 360);

    up_TTO_TH2D hist2D_dEdx_vs_p = std::make_unique<TTO_TH2D>("dEdx_vs_p", "", 200, 0, 11, 200, 0, 1000000);

    up_TTO_TH2D hist2D_ATOF_energy_vs_time = std::make_unique<TTO_TH2D>("ATOF_energy_vs_time", "", 200, 170, 220, 200, 0, 6);
    up_TTO_TH2D hist2D_pe_calc_vs_ATOF_time = std::make_unique<TTO_TH2D>("pe_calc_vs_ATOF_time", "", 200, 0, 1, 200, 150, 200);
    up_TTO_TH2D hist2D_pe_calc_vs_ATOF_energy = std::make_unique<TTO_TH2D>("pe_calc_vs_ATOF_energy", "", 200, 0, 1, 200, 0, 10);
    up_TTO_TH2D hist2D_tof_vs_ATOF_energy = std::make_unique<TTO_TH2D>("tof_vs_ATOF_energy", "", 200, 100, 300, 200, 0, 6);
    up_TTO_TH2D hist2D_tof_vs_ATOF_time = std::make_unique<TTO_TH2D>("tof_vs_ATOF_time", "", 200, 100, 300, 200, 170, 220);


    up_TTO_TH1D hist1D_adc = std::make_unique<TTO_TH1D>("adc", "", 200, 0, 10000);
    up_TTO_TH1D hist1D_t = std::make_unique<TTO_TH1D>("t", "", 200, 0, 1000);
    up_TTO_TH1D hist1D_tot = std::make_unique<TTO_TH1D>("tot", "", 200, 0, 1000);
    up_TTO_TH1D hist1D_ped = std::make_unique<TTO_TH1D>("ped", "", 200, 0, 1000);

    up_TTO_TH1D hist1D_ped_cut = std::make_unique<TTO_TH1D>("ped_cut", "", 200, 0, 1000);
    up_TTO_TH1D hist1D_t_cut = std::make_unique<TTO_TH1D>("t_cut", "", 200, 0, 1000);
    up_TTO_TH1D hist1D_tot_cut = std::make_unique<TTO_TH1D>("tot_cut", "", 200, 0, 1000);
    up_TTO_TH1D hist1D_adc_cut = std::make_unique<TTO_TH1D>("adc_cut", "", 200, 0, 10000);

    up_TTO_TH1D hist1D_ATOF_tdc = std::make_unique<TTO_TH1D>("ATOF_tdc", "", 200, 11000, 14000);
    up_TTO_TH1D hist1D_ATOF_tot = std::make_unique<TTO_TH1D>("ATOF_tot", "", 200, 0, 10000);
    up_TTO_TH2D hist2D_ATOF_tdc_vs_tot = std::make_unique<TTO_TH2D>("ATOF_tdc_vs_tot", "", 200, 11000, 14000, 200, 0, 10000);

};

struct Histograms {
    Electron electron;
};

}  // namespace study1