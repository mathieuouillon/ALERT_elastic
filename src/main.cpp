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
#include <study1/Counter.hpp>
#include <study1/Drawing.hpp>
#include <study1/Histograms.hpp>
#include <study1/Reader.hpp>
#include <thread_pool/multi_thread.hpp>

#include <study2/Histograms.hpp>
#include <study2/Reader.hpp>
#include "fmt/base.h"

void study2_main() {
    ROOT::EnableThreadSafety();  // To stop random errors in multithread mode2

    std::vector<std::string> files = Core::read_recursive_file_in_directory("/volatile/clas12/rg-l/production/auto-pass0/v1/mon/recon/021787/");

    study2::Histograms histograms;
    study2::Reader reader(histograms);
    multithread_reader(reader, files, 1);

    const std::shared_ptr<TH1D> hist1D_vt = histograms.electron.hist1D_vt->Merge();
    const std::shared_ptr<TH1D> hist1D_vtx = histograms.electron.hist1D_vtx->Merge();
    Plotting::draw_hist1D(hist1D_vt, hist1D_vtx, "../plots/study2/", {.legend1 = "vt", .legend2 = "vtx", .log_y = true, .label = "vt [ns]"});

    const std::shared_ptr<TH1D> hist1D_corrected_vt = histograms.electron.hist1D_corrected_vt->Merge();
    const std::shared_ptr<TH1D> hist1D_corrected_vtx = histograms.electron.hist1D_corrected_vtx->Merge();
    Plotting::draw_hist1D(hist1D_corrected_vtx, "../plots/study2/", {.legend1 = "corrected_vt", .log_y = true, .label = "corrected_vt [ns]"});



}

int main() {
    ROOT::EnableThreadSafety();  // To stop random errors in multithread mode

    std::vector<double> times;

    auto start_time = std::chrono::high_resolution_clock::now();

    // Read the configuration file
    const toml::table config = toml::parse_file("../config/study1.toml");

    // Read the files
    std::vector<std::string> files = Core::read_recursive_file_in_directory("/volatile/clas12/ouillon/AI_ALERT/test_workflow-ALERT-data-3/recon/");
    // std::vector<std::string> files = Core::read_recursive_file_in_directory("/volatile/clas12/rg-l/production/auto-pass0/v1/mon/recon/021400/");
    // files.resize(static_cast<int>(20.f / 100 * files.size()));

    fmt::println("First ten files:");
    for (const auto& file : files | std::views::take(10)) {
        fmt::print("{}\n", file);
    }

    fmt::print("Number of files: {}\n", files.size());

    // Process the data
    study1::Counter counter;
    study1::Histograms histograms;
    study1::Reader reader(histograms, config, counter, {11, 22});
    multithread_reader(reader, files, 20);

    study1::Drawing drawing(histograms, config);
    drawing.draw_electron_kinematics();

    fmt::print("Number of events: {}\n", counter.nb_of_event.load());
    fmt::print("Number of events with elastic electron: {}, percentage; {}\n", counter.nb_of_event_with_elastic_electron.load(), counter.nb_of_event_with_elastic_electron.load() / static_cast<double>(counter.nb_of_event.load()) * 100);
    fmt::print("Number of events with AHDCHits: {}, percentage; {}\n", counter.nb_of_event_with_ahdc_hits.load(), counter.nb_of_event_with_ahdc_hits.load() / static_cast<double>(counter.nb_of_event.load()) * 100);
    fmt::print("Number of events with AHDCAiprediction: {}, percentage; {}\n", counter.nb_of_event_with_ahdc_aiprediction.load(), counter.nb_of_event_with_ahdc_aiprediction.load() / static_cast<double>(counter.nb_of_event.load()) * 100);
    fmt::print("Number of events with AHDCAiprediction and track found: {}, percentage; {}\n", counter.nb_of_event_with_ahdc_track_found.load(), counter.nb_of_event_with_ahdc_track_found.load() / static_cast<double>(counter.nb_of_event.load()) * 100);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    fmt::println("Time take: {} milliseconds", duration.count());

    return EXIT_SUCCESS;
}
