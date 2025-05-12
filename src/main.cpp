// C++ headers
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <ranges>
#include <vector>

// fmt headers
#include <fmt/core.h>
#include <fmt/ranges.h>

#include <toml++/toml.hpp>

// ROOT headers
#include <TCanvas.h>
#include <TH1D.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>

// Project headers
#include <Core/Constantes.hpp>
#include <Core/Helpers.hpp>
#include <study1/Drawing.hpp>
#include <study1/Histograms.hpp>
#include <study1/Reader.hpp>
#include <study1/Counter.hpp>
#include <thread_pool/multi_thread.hpp>

int main() {
    ROOT::EnableThreadSafety();  // To stop random errors in multithread mode

    std::vector<double> times;

    auto start_time = std::chrono::high_resolution_clock::now();

    // Read the configuration file
    const toml::table config = toml::parse_file("../config/study1.toml");

    // Read the files
    std::vector<std::string> files = Core::read_recursive_file_in_directory("/volatile/clas12/ouillon/AI_ALERT/test_workflow-ALERT-data-2/recon/");
    // std::vector<std::string> files = Core::read_recursive_file_in_directory("/volatile/clas12/rg-l/production/auto-pass0/v1/mon/recon/021400/");
    // files.resize(static_cast<int>(20.f / 100 * files.size()));

    fmt::print("Number of files: {}\n", files.size());

    // Process the data
    study1::Counter counter;
    study1::Histograms histograms;
    study1::Reader reader(histograms, config, counter, {11, 22});
    multithread_reader(reader, files, 20);

    study1::Drawing drawing(histograms, config);
    drawing.draw_electron_kinematics();

    fmt::print("Number of events: {}\n", counter.nb_of_event.load());
    fmt::print("Number of events with elastic electron: {}, percentage; {}\n", counter.nb_of_event_with_elastic_electron.load(),  counter.nb_of_event_with_elastic_electron.load() / static_cast<double>(counter.nb_of_event.load()) * 100);
    fmt::print("Number of events with AHDCHits: {}, percentage; {}\n", counter.nb_of_event_with_ahdc_hits.load(), counter.nb_of_event_with_ahdc_hits.load() / static_cast<double>(counter.nb_of_event.load()) * 100);
    fmt::print("Number of events with AHDCAiprediction: {}, percentage; {}\n", counter.nb_of_event_with_ahdc_aiprediction.load(), counter.nb_of_event_with_ahdc_aiprediction.load() / static_cast<double>(counter.nb_of_event.load()) * 100);
    fmt::print("Number of events with AHDCAiprediction and track found: {}, percentage; {}\n", counter.nb_of_event_with_ahdc_track_found.load(), counter.nb_of_event_with_ahdc_track_found.load() / static_cast<double>(counter.nb_of_event.load()) * 100);



    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    fmt::println("Time take: {} milliseconds", duration.count());

    return EXIT_SUCCESS;
}
