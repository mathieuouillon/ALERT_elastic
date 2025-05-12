#pragma once

// C++ headers
#include <array>

namespace study1 {

class Line3D {
   public:
    using Point3D = std::array<double, 3>;

    Line3D() = default;

    Line3D(const Point3D& start, const Point3D& end)
        : m_start(start), m_end(end) {}

    Point3D get_point_at_z(double z) const {
        double t = (z - m_start[2]) / (m_end[2] - m_start[2]);
        return Point3D{
            m_start[0] + t * (m_end[0] - m_start[0]),
            m_start[1] + t * (m_end[1] - m_start[1]),
            z};
    }

    Point3D lerp(double t) const {
        return Point3D{
            m_start[0] + t * (m_end[0] - m_start[0]),
            m_start[1] + t * (m_end[1] - m_start[1]),
            m_start[2] + t * (m_end[2] - m_start[2])};
    }

   private:
    Point3D m_start;
    Point3D m_end;
};
}  // namespace study1