#pragma once

#include <atomic>

namespace study1 {

struct Counter {
    std::atomic_int nb_of_event{};
    std::atomic_int nb_of_event_with_elastic_electron{};
    std::atomic_int nb_of_event_with_ahdc_hits{};
    std::atomic_int nb_of_event_with_ahdc_aiprediction{};
    std::atomic_int nb_of_event_with_ahdc_track_found{};
};

}  // namespace study1