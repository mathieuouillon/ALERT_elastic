#pragma once

#include <study2/Reader.hpp>

namespace study2 {

Reader::Reader(Histograms& histograms)
    : m_histograms(histograms) {
}

Reader::~Reader() = default;

auto Reader::operator()(const std::string& file) -> void {

    for (hipo::hipoeventfile events(file); auto event : events) {

        hipo::bank& REC_Particle = event.get_bank("REC::Particle");
        hipo::bank& REC_Event = event.get_bank("REC::Event");

        double RF_time = REC_Event.get<double>("RFTime", 0); // ns
        
        for (int i = 0; i < REC_Particle.getRows(); i++) {
            const int pid = REC_Particle.get<int>("pid", i);
            const double px = REC_Particle.get<double>("px", i);
            const double py = REC_Particle.get<double>("py", i);
            const double pz = REC_Particle.get<double>("pz", i);
            const double vx = REC_Particle.get<double>("vx", i);
            const double vy = REC_Particle.get<double>("vy", i);
            const double vz = REC_Particle.get<double>("vz", i);
            const double vt = REC_Particle.get<double>("vt", i);

            if (pid == 11) {
                m_histograms.electron.hist1D_vt->Get()->Fill(vt);

                Core::ScintillatorBank scintillator_bank = Core::read_Scintillator_bank(event.get_bank("REC::Scintillator"), i);
                
                double time = scintillator_bank.time;
                double path = scintillator_bank.path;
                double energy = scintillator_bank.energy;
                double speedOfLight = 29.9792458;
                double p = std::sqrt(px*px + py*py + pz*pz);
                double m_e = 0.000511; // GeV

                double beta = p / std::sqrt(p*p + std::pow(m_e,2));
                double vtx = time - path/speedOfLight/beta;
                
                m_histograms.electron.hist1D_vtx->Get()->Fill(vtx);

                m_histograms.electron.hist1D_corrected_vt->Get()->Fill(vt - RF_time);
                m_histograms.electron.hist1D_corrected_vtx->Get()->Fill(vtx - RF_time);


            }
        }



    }
}
}