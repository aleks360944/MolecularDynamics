#pragma once

#include <cmath>
#include <stdexcept>
#include "Vec2.hpp"

inline void validateLennardJonesParameters(double epsilon, double sigma) {
    if (!std::isfinite(epsilon) || epsilon < 0.0 ||
        !std::isfinite(sigma) || sigma <= 0.0) {
        throw std::invalid_argument("Lennard-Jones epsilon must be non-negative and sigma must be positive");
    }
}

inline Vec2 calculateLennardJonesForceFromDelta(
    double dx,
    double dy,
    double distanceSquared,
    double epsilon,
    double sigma
) {
    validateLennardJonesParameters(epsilon, sigma);
    if (!std::isfinite(distanceSquared) || distanceSquared <= 0.0) {
        throw std::domain_error("Lennard-Jones force requires distinct finite positions");
    }

    const double ratio2 = (sigma * sigma) / distanceSquared;
    const double ratio6 = ratio2 * ratio2 * ratio2;
    const double ratio12 = ratio6 * ratio6;
    const double factor = (24.0 * epsilon / distanceSquared)
                        * (ratio6 - 2.0 * ratio12);

    return {factor * dx, factor * dy};
}

inline Vec2 calculateLennardJonesForce(
    const Vec2& positionA,
    const Vec2& positionB,
    double epsilon,
    double sigma
) {
    const double dx = positionB.x - positionA.x;
    const double dy = positionB.y - positionA.y;

    const double distanceSquared = dx * dx + dy * dy;
    return calculateLennardJonesForceFromDelta(
        dx, dy, distanceSquared, epsilon, sigma
    );
}
