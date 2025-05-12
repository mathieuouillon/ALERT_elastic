#pragma once
// C++ headers
#include <cmath>
#include <utility>
// Project headers
#include <Core/Constantes.hpp>
#include <Core/Helpers.hpp>
#include <study1/Line3D.hpp>

namespace study1 {

class Hit {
   public:
    bool m_used = false;
    bool m_visited = false;

    int m_id;
    int m_super_layer;
    int m_layer;
    int m_wire;
    int m_wire_nb;
    int m_clusterId = -1;  // -1 means unclassified

    double m_x;
    double m_y;
    double m_r;
    Line3D m_line;

    Hit(int id, int superlayer, int layer, int wire);

    Hit() = default;                                 // Default constructor
    Hit(const Hit& other) = default;                 // Copy constructor
    Hit(Hit&& other) noexcept = default;             // Move constructor
    Hit& operator=(const Hit& other) = default;      // Copy assignment operator
    Hit& operator=(Hit&& other) noexcept = default;  // Move assignment operator
    ~Hit() = default;
};

inline Hit::Hit(int id, int superlayer, int layer, int wire)
    : m_id(id), m_super_layer(superlayer), m_layer(layer), m_wire(wire) {
    auto [wire_nb, layer_radius] = Core::get_info_superlayer(superlayer);
    double dR = 4.0;
    double R_layer = layer_radius + dR * layer;
    double wirePhiIndex = wire + 0.5 * (wire % 2) + 0.5 * layer * (1 - 2 * (wire_nb % 2));
    double thster = Core::to_radians(-20.0);
    double alphaW_layer = Core::to_radians(360.0 / wire_nb);

    double x_start = R_layer * std::cos(alphaW_layer * wirePhiIndex);
    double y_start = R_layer * std::sin(alphaW_layer * wirePhiIndex);
    double x_end = R_layer * std::cos(alphaW_layer * wirePhiIndex + thster * (std::pow(-1, superlayer)));
    double y_end = R_layer * std::sin(alphaW_layer * wirePhiIndex + thster * (std::pow(-1, superlayer)));

    m_line = Line3D({x_start, y_start, -300.0 / 2}, {x_end, y_end, 300.0 / 2});
    m_wire_nb = wire_nb;
    m_x = x_start;
    m_y = y_start;
    m_r = std::hypot(m_x, m_y);
}

}  // namespace study1