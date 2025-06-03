#pragma once

// C++ headers
#include <string>
#include <vector>

// <fmt> headers
#include <fmt/core.h>

// toml++ headers
#include <toml++/toml.h>

// Project headers
#include <hipo4/hipoeventiterator.h>
#include <Core/Constantes.hpp>
#include <Core/Helpers.hpp>
#include <Core/Particle.hpp>
#include <Core/ReadBank.hpp>
#include <study2/Histograms.hpp>


namespace study2 {

class Reader {
   private:
    // ****** private constants
    static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

    // ****** private variables
    Histograms& m_histograms;

   public:
    // ****** constructors and destructor
    explicit Reader(Histograms& histograms);
    ~Reader();

    // ****** public methods
    auto operator()(const std::string& file) -> void;
};

}  // namespace study1