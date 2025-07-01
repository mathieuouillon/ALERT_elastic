#include <study1/Drawing.hpp>

namespace study1 {

Drawing::Drawing(Histograms& histograms, const toml::parse_result& config)
    : m_histograms(histograms), m_config(config) {
}

auto Drawing::draw_electron_kinematics() -> void {

    // Cuts ---------------------------------------------------------------------------------------
    const double vz_min_e = m_config["electron"]["vz_min"].value_or(NaN);
    const double vz_max_e = m_config["electron"]["vz_max"].value_or(NaN);
    const double chi2_min_e = m_config["electron"]["chi2_min"].value_or(NaN);
    const double chi2_max_e = m_config["electron"]["chi2_max"].value_or(NaN);
    //---------------------------------------------------------------------------------------------

    // Merge 1D histograms ------------------------------------------------------------------------
    const std::shared_ptr<TH1D> hist1D_p = m_histograms.electron.hist1D_p->Merge();
    const std::shared_ptr<TH1D> hist1D_chi2 = m_histograms.electron.hist1D_chi2->Merge();
    const std::shared_ptr<TH1D> hist1D_phi = m_histograms.electron.hist1D_phi->Merge();
    const std::shared_ptr<TH1D> hist1D_theta = m_histograms.electron.hist1D_theta->Merge();
    const std::shared_ptr<TH1D> hist1D_vz = m_histograms.electron.hist1D_vz->Merge();

    const std::shared_ptr<TH1D> hist1D_p_cut = m_histograms.electron.hist1D_p_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_chi2_cut = m_histograms.electron.hist1D_chi2_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_phi_cut = m_histograms.electron.hist1D_phi_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_theta_cut = m_histograms.electron.hist1D_theta_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_vz_cut = m_histograms.electron.hist1D_vz_cut->Merge();
    // --------------------------------------------------------------------------------------------

    Plotting::draw_hist1D(hist1D_p, hist1D_p_cut, m_path_electron, {.label = "p_{e} [GeV/c]"});
    Plotting::draw_hist1D(hist1D_chi2, hist1D_chi2_cut, m_path_electron, {.cuts = {chi2_min_e, chi2_max_e}, .label = "chi2pid_{e}"});
    Plotting::draw_hist1D(hist1D_phi, hist1D_phi_cut, m_path_electron, {.label = "#phi_{e} [deg.]"});
    Plotting::draw_hist1D(hist1D_theta, hist1D_theta_cut, m_path_electron, {.label = "#theta_{e} [deg.]"});
    Plotting::draw_hist1D(hist1D_vz, hist1D_vz_cut, m_path_electron, {.cuts = {vz_min_e, vz_max_e}, .label = "Vz_{e} [cm]"});


    // Electron based kinematics ------------------------------------------------------------------
    const std::shared_ptr<TH1D> hist1D_Q2 = m_histograms.electron.hist1D_Q2->Merge();
    const std::shared_ptr<TH1D> hist1D_W = m_histograms.electron.hist1D_W->Merge();
    const std::shared_ptr<TH1D> hist1D_xB = m_histograms.electron.hist1D_xB->Merge();
    const std::shared_ptr<TH1D> hist1D_MM_eP = m_histograms.electron.hist1D_MM_eP->Merge();

    Plotting::draw_hist1D(hist1D_Q2, m_path_electron, {.label = "Q2 [GeV^2]"});
    Plotting::draw_hist1D(hist1D_W, m_path_electron, {.label = "W [GeV]"});
    Plotting::draw_hist1D(hist1D_xB, m_path_electron, {.label = "xB"});
    Plotting::draw_hist1D(hist1D_MM_eP, m_path_electron, {.label = "MM_{eP} [GeV^2]"});
    // --------------------------------------------------------------------------------------------


    // Delta phi ----------------------------------------------------------------------------------
    const std::shared_ptr<TH1D> hist1D_delta_phi = m_histograms.electron.hist1D_delta_phi->Merge();
    const std::shared_ptr<TH2D> hist2D_delta_phi_vs_phiElec = m_histograms.electron.hist2D_delta_phi_vs_phiElec->Merge();
    
    Plotting::draw_hist1D(hist1D_delta_phi, m_path_electron, {.label = "#Delta #phi [deg.]"});
    Plotting::draw_hist2D(hist2D_delta_phi_vs_phiElec, m_path_electron, {.label_x = "#Delta #phi [deg.]", .label_y = "#phi_{e} [deg.]"});
    // --------------------------------------------------------------------------------------------


    // Proton variables ------------------------------------------------------------------
    const std::shared_ptr<TH1D> hist1D_pp_calc = m_histograms.electron.hist1D_pp_calc->Merge();
    Plotting::draw_hist1D(hist1D_pp_calc, m_path_electron, {.label = "p_{p} [GeV/c]"});
    // --------------------------------------------------------------------------------------------


    // AHDC variables ------------------------------------------------------------------
    const std::shared_ptr<TH1D> hist1D_adc = m_histograms.electron.hist1D_adc->Merge();
    const std::shared_ptr<TH1D> hist1D_t = m_histograms.electron.hist1D_t->Merge();
    const std::shared_ptr<TH1D> hist1D_tot = m_histograms.electron.hist1D_tot->Merge();
    const std::shared_ptr<TH1D> hist1D_ped = m_histograms.electron.hist1D_ped->Merge();

    const std::shared_ptr<TH1D> hist1D_adc_cut = m_histograms.electron.hist1D_adc_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_t_cut = m_histograms.electron.hist1D_t_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_tot_cut = m_histograms.electron.hist1D_tot_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_ped_cut = m_histograms.electron.hist1D_ped_cut->Merge();

    Plotting::draw_hist1D(hist1D_adc, m_path_electron, {.label = "ADC"});
    Plotting::draw_hist1D(hist1D_t, m_path_electron, {.label = "t [ns]"});
    Plotting::draw_hist1D(hist1D_tot, m_path_electron, {.label = "TOT [ns]"});
    Plotting::draw_hist1D(hist1D_ped, m_path_electron, {.label = "ped"});

    Plotting::draw_hist1D(hist1D_adc_cut, m_path_electron, {.label = "ADC"});
    Plotting::draw_hist1D(hist1D_t_cut, m_path_electron, {.label = "t [ns]"});
    Plotting::draw_hist1D(hist1D_tot_cut, m_path_electron, {.label = "TOT [ns]"});
    Plotting::draw_hist1D(hist1D_ped_cut, m_path_electron, {.label = "ped"});


    // ATOF variables ------------------------------------------------------------------

    const std::shared_ptr<TH2D> hist2D_ATOF_energy_vs_time = m_histograms.electron.hist2D_ATOF_energy_vs_time->Merge();
    Plotting::draw_hist2D(hist2D_ATOF_energy_vs_time, m_path_electron, {.label_x = "time [ns]", .label_y = "energy [MeV]"});


    const std::shared_ptr<TH1D> hist1D_ATOF_tdc = m_histograms.electron.hist1D_ATOF_tdc->Merge();
    const std::shared_ptr<TH1D> hist1D_ATOF_tot = m_histograms.electron.hist1D_ATOF_tot->Merge();
    const std::shared_ptr<TH2D> hist2D_ATOF_tdc_vs_tot = m_histograms.electron.hist2D_ATOF_tdc_vs_tot->Merge();
    Plotting::draw_hist1D(hist1D_ATOF_tdc, m_path_electron, {.label = "ATOF TDC [ns]"});
    Plotting::draw_hist1D(hist1D_ATOF_tot, m_path_electron, {.label = "ATOF ToT [ns]"});
    Plotting::draw_hist2D(hist2D_ATOF_tdc_vs_tot, m_path_electron, {.label_x = "ATOF TDC [ns]", .label_y = "ATOF ToT [ns]"});


    const std::shared_ptr<TH2D> hist2D_pe_calc_vs_ATOF_time = m_histograms.electron.hist2D_pe_calc_vs_ATOF_time->Merge();
    const std::shared_ptr<TH2D> hist2D_pe_calc_vs_ATOF_energy = m_histograms.electron.hist2D_pe_calc_vs_ATOF_energy->Merge();
    const std::shared_ptr<TH2D> hist2D_tof_vs_ATOF_energy = m_histograms.electron.hist2D_tof_vs_ATOF_energy->Merge();
    const std::shared_ptr<TH2D> hist2D_tof_vs_ATOF_time = m_histograms.electron.hist2D_tof_vs_ATOF_time->Merge();
    Plotting::draw_hist2D(hist2D_pe_calc_vs_ATOF_time, m_path_electron, {.label_x = "p_{4He}^{calc} [GeV/c]", .label_y = "ATOF time [ns]"});
    Plotting::draw_hist2D(hist2D_pe_calc_vs_ATOF_energy, m_path_electron, {.label_x = "p_{4He}^{calc} [GeV/c]", .label_y = "ATOF energy [MeV]"});
    Plotting::draw_hist2D(hist2D_tof_vs_ATOF_energy, m_path_electron, {.label_x = "expected time [ns]", .label_y = "ATOF energy [MeV]"});
    Plotting::draw_hist2D(hist2D_tof_vs_ATOF_time, m_path_electron, {.label_x = "expected time [ns]", .label_y = "ATOF time [ns]"});


}


}  // namespace study1