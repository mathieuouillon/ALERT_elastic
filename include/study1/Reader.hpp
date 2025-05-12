#pragma once

// C++ headers
#include <string>

// <fmt> headers
#include <fmt/core.h>

// toml++ headers
#include <toml++/toml.h>

// Project headers
#include <hipo4/hipoeventiterator.h>
#include <hipo4/reader.h>
#include <hipo4/writer.h>
#include <Core/Constantes.hpp>
#include <Core/Helpers.hpp>
#include <Core/Particle.hpp>
#include <Core/ReadBank.hpp>
#include <study1/Counter.hpp>
#include <study1/Histograms.hpp>
#include <study1/Hit.hpp>
#include <study1/Line3D.hpp>
#include <vector>

namespace study1 {

class Reader {
   private:
    // ****** private variables

    // Define Constants here
    static constexpr bool DEBUG = false;
    static constexpr bool PRINT = true;

    Histograms& m_histograms;
    const toml::parse_result& m_config;
    Counter& m_counter;
    std::vector<int> m_pids;

    // ****** private constants
    static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

    // ****** private methods
    auto select_electron(const Core::Particle& electron) const -> bool;
    auto get_topology(const hipo::bank& REC_Particle, std::unordered_map<int, std::vector<Core::Particle>>& particle_collections) -> void;

   public:
    // ****** constructors and destructor
    explicit Reader(Histograms& histograms, const toml::parse_result& config, Counter& counter, const std::vector<int>& pids = {11});
    ~Reader();

    // ****** public methods
    auto operator()(const std::string& file) -> void;
};

}  // namespace study1