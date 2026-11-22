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

inline double calculateForceShiftedLennardJonesPotential(
    double distanceSquared,
    double epsilon,
    double sigma,
    double cutoffFactor = 2.5
) {
    validateLennardJonesParameters(epsilon, sigma);
    if (!std::isfinite(distanceSquared) || distanceSquared <= 0.0 ||
        !std::isfinite(cutoffFactor) || cutoffFactor <= 0.0) {
        throw std::domain_error("Force-shifted Lennard-Jones potential requires valid distance and cutoff");
    }

    const double distance = std::sqrt(distanceSquared);
    const double cutoff = cutoffFactor * sigma;
    if (!std::isfinite(cutoff)) {
        throw std::overflow_error("Lennard-Jones cutoff is not finite");
    }
    if (distance >= cutoff) {
        return 0.0;
    }

    const double ratio = sigma / distance;
    const double ratio2 = ratio * ratio;
    const double ratio6 = ratio2 * ratio2 * ratio2;
    const double ratio12 = ratio6 * ratio6;
    const double potential = 4.0 * epsilon * (ratio12 - ratio6);

    const double cutoffRatio = sigma / cutoff;
    const double cutoffRatio2 = cutoffRatio * cutoffRatio;
    const double cutoffRatio6 = cutoffRatio2 * cutoffRatio2 * cutoffRatio2;
    const double cutoffRatio12 = cutoffRatio6 * cutoffRatio6;
    const double cutoffPotential =
        4.0 * epsilon * (cutoffRatio12 - cutoffRatio6);
    const double cutoffDerivative =
        (24.0 * epsilon / cutoff) * (cutoffRatio6 - 2.0 * cutoffRatio12);

    return potential - cutoffPotential - (distance - cutoff) * cutoffDerivative;
}

inline Vec2 calculateForceShiftedLennardJonesForceFromDelta(
    double dx,
    double dy,
    double distanceSquared,
    double epsilon,
    double sigma,
    double cutoffFactor = 2.5
) {
    validateLennardJonesParameters(epsilon, sigma);
    if (!std::isfinite(distanceSquared) || distanceSquared <= 0.0 ||
        !std::isfinite(cutoffFactor) || cutoffFactor <= 0.0) {
        throw std::domain_error("Force-shifted Lennard-Jones force requires valid distance and cutoff");
    }

    const double distance = std::sqrt(distanceSquared);
    const double cutoff = cutoffFactor * sigma;
    if (!std::isfinite(cutoff)) {
        throw std::overflow_error("Lennard-Jones cutoff is not finite");
    }
    if (distance >= cutoff) {
        return {0.0, 0.0};
    }

    const double ratio2 = (sigma * sigma) / distanceSquared;
    const double ratio6 = ratio2 * ratio2 * ratio2;
    const double ratio12 = ratio6 * ratio6;
    const double radialDerivative =
        (24.0 * epsilon / distance) * (ratio6 - 2.0 * ratio12);

    const double cutoffRatio = sigma / cutoff;
    const double cutoffRatio2 = cutoffRatio * cutoffRatio;
    const double cutoffRatio6 = cutoffRatio2 * cutoffRatio2 * cutoffRatio2;
    const double cutoffRatio12 = cutoffRatio6 * cutoffRatio6;
    const double cutoffDerivative =
        (24.0 * epsilon / cutoff) * (cutoffRatio6 - 2.0 * cutoffRatio12);

    const double factor = (radialDerivative - cutoffDerivative) / distance;
    return {factor * dx, factor * dy};
}
