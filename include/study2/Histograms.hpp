#pragma once

#include <TH1.h>
#include <TH2.h>
#include <ROOT/TThreadedObject.hxx>

using TTO_TH1D = ROOT::TThreadedObject<TH1D>;
using TTO_TH2D = ROOT::TThreadedObject<TH2D>;
using up_TTO_TH1D = const std::unique_ptr<TTO_TH1D>;
using up_TTO_TH2D = const std::unique_ptr<TTO_TH2D>;

namespace study2 {

struct Electron {
    up_TTO_TH1D hist1D_vt = std::make_unique<TTO_TH1D>("vt", "", 200, 0, 200);
    up_TTO_TH1D hist1D_vtx = std::make_unique<TTO_TH1D>("vtx", "", 200, 0, 200);

    up_TTO_TH1D hist1D_corrected_vt = std::make_unique<TTO_TH1D>("corrected_vt", "", 400, -200, 200);
    up_TTO_TH1D hist1D_corrected_vtx = std::make_unique<TTO_TH1D>("corrected_vtx", "", 400, 0, 50);

};

struct Histograms {
    Electron electron;
};

}  // namespace study2