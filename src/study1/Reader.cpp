#include <set>
#include <study1/Reader.hpp>

namespace study1 {

Reader::Reader(Histograms& histograms, const toml::parse_result& config, Counter& counter, const std::vector<int>& pids)
    : m_histograms(histograms), m_config(config), m_counter(counter), m_pids(pids) {
}

Reader::~Reader() = default;

auto Reader::operator()(const std::string& file) -> void {
    bool write = true;
    // fmt::print("Processing file: {}\n", file);
    hipo::hipoeventfile events(file);

    auto writer = hipo::writer();
    auto& dictionary = events.get_dictionary();

    for (const std::string& s : dictionary.getSchemaList())
        writer.getDictionary().addSchema(dictionary.getSchema(s.c_str()));

    std::vector<std::string> dst_schema = {"AHDC::adc", "RUN::config", "RUN::scaler", "REC::Event", "REC::Particle", "REC::Calorimeter",
                                           "REC::ForwardTagger", "REC::Scintillator", "REC::Track", "REC::CovMat", "REC::Traj", "REC::Cherenkov", "RUN::config", "RUN::rf"};

    const std::filesystem::path path = file;
    const std::string name = path.filename().string();
    const std::string base_name = name.substr(0, name.find_last_of('.'));
    if (write) writer.open(std::string("../hipo_out/" + base_name + ".hipo").c_str());

    std::unordered_map<int, std::vector<Core::Particle>> particle_collections;
    for (const auto& pid : m_pids) {
        particle_collections[pid] = std::vector<Core::Particle>();
    }

    // Open a txt file to write the results
    std::string file_name = file.substr(file.find_last_of("/") + 1);
    std::ofstream output_file("../output/" + file_name + ".txt");
    if (!output_file.is_open()) {
        std::cerr << "Error opening output file: " << file << ".txt" << std::endl;
        return;
    }
    int nb_of_good_events = 0;
    for (auto event : events) {

        for (auto& [key, vec] : particle_collections) {
            vec.clear();
        }

        hipo::bank& REC_Particle = event.get_bank("REC::Particle");
        hipo::bank& REC_Calorimeter = event.get_bank("REC::Calorimeter");
        hipo::bank& REC_Cherenkov = event.get_bank("REC::Cherenkov");

        if (REC_Particle.getRows() == 0) continue;
        m_counter.nb_of_event++;

        bool is_a_good_event = false;

        get_topology(REC_Particle, particle_collections);
        std::vector<Core::Particle>& electrons = particle_collections[11];
        std::vector<Core::Particle>& photons = particle_collections[22];

        // Check if photons are in the Forward Tagger
        std::vector<Core::Particle> photons_in_ft;
        for (auto& photon : photons) {
            if (1000 <= photon.status() && photon.status() < 2000) photons_in_ft.push_back(photon);
        }

        // Check if electrons are in the Forward Detector
        std::vector<Core::Particle> electrons_in_fd;
        for (auto& electron : electrons) {
            if (2000 <= std::abs(electron.status()) && std::abs(electron.status()) < 4000) electrons_in_fd.push_back(electron);
        }

        // electrons.insert(electrons.end(), photons_in_ft.begin(), photons_in_ft.end());

        int nb_elastic = 0;
        std::vector<Core::Particle> elastic_electrons;

        for (auto& electron : electrons) {
            if (select_electron(electron)) {
                double E = electron.E();
                double Ebeam = 2.23951;
                double Q2 = 2.0 * E * Ebeam * (1.0 - std::cos(electron.theta() * M_PI / 180.0));
                double W = std::sqrt(0.938272 * 0.938272 + 2.0 * 0.938272 * (Ebeam - E) - Q2);
                double xB = Q2 / (2.0 * 0.938272 * (Ebeam - E));

                ROOT::Math::PxPyPzEVector beam(0.0, 0.0, Ebeam, Ebeam);
                ROOT::Math::PxPyPzEVector k2 = electron.PxPyPzEVector();
                ROOT::Math::PxPyPzEVector p1 = {0, 0, 0, 0.938272};
                double MM_eP = (beam + p1 - k2).M2();

                m_histograms.electron.hist1D_W->Get()->Fill(W);

                if (0.9 < W && W < 1.12) {
                    nb_elastic++;
                    elastic_electrons.push_back(electron);

                    m_histograms.electron.hist1D_Q2->Get()->Fill(Q2);
                    m_histograms.electron.hist1D_xB->Get()->Fill(xB);
                    m_histograms.electron.hist1D_MM_eP->Get()->Fill(MM_eP);

                    // Expected proton variables
                    double mass = 0.938272;
                    double El_Theta = electron.theta();
                    double pe_calc = Ebeam / (1 + (Ebeam / mass) * (1 - cos(El_Theta * TMath::Pi() / 180.)));
                    double pp_calc = pe_calc * TMath::Sqrt(pow((1 + Ebeam / mass), 2) * pow((1 - cos(El_Theta * TMath::Pi() / 180.)), 2) + pow(sin(El_Theta * TMath::Pi() / 180.), 2));

                    m_histograms.electron.hist1D_pp_calc->Get()->Fill(pp_calc);
                }
            }
        }

        if (nb_elastic != 1) continue;
        Core::Particle& elastic_electron = elastic_electrons[0];
        m_counter.nb_of_event_with_elastic_electron++;

        // todo try with adc
        hipo::bank& AHDC_adc = event.get_bank("AHDC::adc");
        hipo::bank& AHDC_hits = event.get_bank("AHDC::hits");
        hipo::bank& AHDC_aiprediction = event.get_bank("AHDC::ai:prediction");
        hipo::bank& AHDC_Track = event.get_bank("AHDC::track");
        hipo::bank& ATOF_hits = event.get_bank("ATOF::hits");

        if (AHDC_adc.getRows() == 0) continue;
        for (int i = 0; i < AHDC_adc.getRows(); i++) {
            int layer = AHDC_adc.get<int>("layer", i);
            int component = AHDC_adc.get<int>("component", i);
            int ADC = AHDC_adc.get<int>("ADC", i);
            int integral = AHDC_adc.get<int>("integral", i);
            double time = AHDC_adc.get<float>("time", i);
            double leadingEdgeTime = AHDC_adc.get<float>("leadingEdgeTime", i);
            double timeOverThreshold = AHDC_adc.get<float>("timeOverThreshold", i);
            double constantFractionTime = AHDC_adc.get<float>("constantFractionTime", i);

            int ahdc_superlayer = layer / 10 - 1;
            int ahdc_layer = layer % 10 - 1;
            int ahdc_wire = component - 1;

            double adcOffset = AHDC_adc.get<float>("ped", i);

            double adc = ADC;
            double adc_min = 0;
            double adc_max = 4095;
            double t_min = 200;
            double t_max = 500;
            double tot_min = 350;
            double tot_max = 600;
            double ped_min = 180;
            double ped_max = 360;
            if (!((adc >= adc_min) && (adc <= adc_max) && (leadingEdgeTime >= t_min) && (leadingEdgeTime <= t_max) && (timeOverThreshold >= tot_min) && (timeOverThreshold <= tot_max) && (adcOffset >= ped_min) && (adcOffset <= ped_max))) continue;
            //if (ahdc_superlayer != 0) continue;

            auto [wire_nb, layer_radius] = Core::get_info_superlayer(ahdc_superlayer);

            double dR = 4.0;
            double R_layer = layer_radius + dR * ahdc_layer;
            double wirePhiIndex = ahdc_wire + 0.5 * (wire_nb % 2) + 0.5 * ahdc_layer * (1 - 2 * (wire_nb % 2));
            double thster = Core::to_radians(-20.0);
            double alphaW_layer = Core::to_radians(360.0 / wire_nb);

            double x_start = R_layer * std::cos(alphaW_layer * wirePhiIndex);
            double y_start = R_layer * std::sin(alphaW_layer * wirePhiIndex);
            double x_end = R_layer * std::cos(alphaW_layer * wirePhiIndex + thster * (std::pow(-1, ahdc_superlayer)));
            double y_end = R_layer * std::sin(alphaW_layer * wirePhiIndex + thster * (std::pow(-1, ahdc_superlayer)));

            Line3D line({x_start, y_start, -300.0 / 2}, {x_end, y_end, 300.0 / 2});
            auto [x, y, z] = line.lerp((elastic_electron.vz() * 10 + 230) / 300);
            double phi = x == 0.0 && y == 0.0 ? 0.0 : TMath::ATan2(y, x);
            double phiX = Core::phi_to_range(Core::to_radians(elastic_electron.phi()) + TMath::Pi());
            double delta_phi = Core::phi_to_range(phi - phiX);

            m_histograms.electron.hist1D_delta_phi->Get()->Fill(Core::to_degrees(delta_phi));
            m_histograms.electron.hist2D_delta_phi_vs_phiElec->Get()->Fill(Core::to_degrees(delta_phi), elastic_electron.phi());

            if (std::abs(delta_phi) < 220.0) {
                m_histograms.electron.hist1D_delta_phi_cut->Get()->Fill(Core::to_degrees(delta_phi));

                is_a_good_event = true;
                nb_of_good_events++;
            }
        }

        if (is_a_good_event) {

            std::vector<Hit> hits;
            for (int i = 0; i < AHDC_adc.getRows(); i++) {
                int layer = AHDC_adc.get<int>("layer", i);
                int component = AHDC_adc.get<int>("component", i);
                int ADC = AHDC_adc.get<int>("ADC", i);
                int integral = AHDC_adc.get<int>("integral", i);
                double time = AHDC_adc.get<float>("time", i);
                double leadingEdgeTime = AHDC_adc.get<float>("leadingEdgeTime", i);
                double timeOverThreshold = AHDC_adc.get<float>("timeOverThreshold", i);
                double constantFractionTime = AHDC_adc.get<float>("constantFractionTime", i);

                int ahdc_superlayer = layer / 10 - 1;
                int ahdc_layer = layer % 10 - 1;
                int ahdc_wire = component - 1;
                hits.emplace_back(i, ahdc_superlayer, ahdc_layer, ahdc_wire);
            }

            /*
            std::set<int> ahdc_superlayers;
            for (int i = 0; i < AHDC_adc.getRows(); i++) {
                int layer = AHDC_adc.get<int>("layer", i);
                int ahdc_superlayer = layer / 10 - 1;
                ahdc_superlayers.insert(ahdc_superlayer);
            }

            // fmt::println("AHDCHits: {}", ahdc_superlayers);
            if (ahdc_superlayers.size() == 5) {
                m_counter.nb_of_event_with_ahdc_hits++;
            }




            if (AHDC_aiprediction.getRows() != 0) {
                m_counter.nb_of_event_with_ahdc_aiprediction++;

                for (int i = 0; i < AHDC_aiprediction.getRows(); i++) {
                    double pred = AHDC_aiprediction.get<double>("pred", i);
                    if (pred > 0.2) {
                        m_counter.nb_of_event_with_ahdc_track_found++;
                    }
                }
            } 
            */
        }

        if (is_a_good_event && write) {
            hipo::event outEvent = *event.get_event_ptr();
            outEvent.reset();

            for (const std::string& s : dst_schema) {
                auto bank = hipo::bank(dictionary.getSchema(s.c_str()));
                event.get_event_ptr()->getStructure(bank);
                if (bank.getRows() > 0) {
                    outEvent.addStructure(bank);
                }
            }

            writer.addEvent(outEvent);
        }
    }
    if (write) writer.close();

    if (nb_of_good_events == 0) std::filesystem::remove("../hipo_out/" + base_name + ".hipo");

    output_file.close();
}

auto Reader::get_topology(const hipo::bank& REC_Particle, std::unordered_map<int, std::vector<Core::Particle>>& particle_collections) -> void {
    const int rows = REC_Particle.getRows();

    // Pre-allocate for each particle types
    for (auto& [key, vec] : particle_collections) {
        vec.reserve(rows);
    }

    for (int i = 0; i < REC_Particle.getRows(); i++) {
        const int pid = REC_Particle.get<int>("pid", i);

        if (!particle_collections.contains(pid)) continue;

        // Get particle properties only after we know we want this particle
        const int status = REC_Particle.get<int>("status", i);
        const int charge = REC_Particle.get<int>("charge", i);
        const double px = REC_Particle.get<double>("px", i);
        const double py = REC_Particle.get<double>("py", i);
        const double pz = REC_Particle.get<double>("pz", i);
        const double vx = REC_Particle.get<double>("vx", i);
        const double vy = REC_Particle.get<double>("vy", i);
        const double vz = REC_Particle.get<double>("vz", i);
        const double vt = REC_Particle.get<double>("vt", i);
        const double beta = REC_Particle.get<double>("beta", i);
        const double chi2pid = REC_Particle.get<double>("chi2pid", i);

        const double mass = Core::get_mass(pid);
        const double E = Core::compute_energy(px, py, pz, pid);

        particle_collections[pid].emplace_back(pid, status, i, charge, mass, px, py, pz, E, vx, vy, vz, vt, beta, chi2pid);
    }
}

auto Reader::select_electron(const Core::Particle& electron) const -> bool {
    bool pass = false;

    // Cuts ---------------------------------------------------------------------------------------
    const double vz_min_e = m_config["electron"]["vz_min"].value_or(NaN);
    const double vz_max_e = m_config["electron"]["vz_max"].value_or(NaN);
    const double chi2_min_e = m_config["electron"]["chi2_min"].value_or(NaN);
    const double chi2_max_e = m_config["electron"]["chi2_max"].value_or(NaN);
    //---------------------------------------------------------------------------------------------

    // Get variables for electrons kinematics to apply the cuts -----------------------------------
    const double p_e = electron.p();
    const double chi2_e = electron.chi2pid();
    const double vz_e = electron.vz();
    const double phi_e = electron.phi();
    const double theta_e = electron.theta();
    //---------------------------------------------------------------------------------------------

    // Cuts ---------------------------------------------------------------------------------------
    const bool cut_chi2 = chi2_min_e < chi2_e && chi2_e < chi2_max_e;
    const bool cut_vz = vz_min_e <= vz_e && vz_e <= vz_max_e;
    const bool cuts = cut_chi2 && cut_vz;
    // --------------------------------------------------------------------------------------------

    // Fill the histograms before any cuts are applied --------------------------------------------
    m_histograms.electron.hist1D_p->Get()->Fill(p_e);
    m_histograms.electron.hist1D_phi->Get()->Fill(phi_e);
    m_histograms.electron.hist1D_theta->Get()->Fill(theta_e);
    m_histograms.electron.hist1D_chi2->Get()->Fill(chi2_e);
    m_histograms.electron.hist1D_vz->Get()->Fill(vz_e);

    //---------------------------------------------------------------------------------------------

    // Fill histograms after cuts: ----------------------------------------------------------------
    if (cut_chi2) m_histograms.electron.hist1D_vz_cut->Get()->Fill(vz_e);
    if (cut_vz) m_histograms.electron.hist1D_chi2_cut->Get()->Fill(chi2_e);

    if (cuts) {
        m_histograms.electron.hist1D_p_cut->Get()->Fill(p_e);
        m_histograms.electron.hist1D_phi_cut->Get()->Fill(phi_e);
        m_histograms.electron.hist1D_theta_cut->Get()->Fill(theta_e);
        pass = true;
    }
    //---------------------------------------------------------------------------------------------

    return pass;
}

}  // namespace study1