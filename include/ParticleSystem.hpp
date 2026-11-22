#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Particle.hpp"
#include "LennardJones.hpp"

namespace SimulationWorld {
inline constexpr double Width = 19.0;
inline constexpr double Height = 10.0;
inline constexpr double PixelsPerUnit = 100.0;
}

struct ParticleSystem {
    struct Pair {
        std::size_t first;
        std::size_t second;
    };

    std::vector<Particle> particles;
    std::vector<Vec2> previousAccelerations;
    std::vector<Vec2> velocitySigns;
    bool forcesInitialized = false;

    void addParticle(const Particle& particle) {
        if (!particle.isValid() || particle.position.x < 0.0 ||
            particle.position.x > SimulationWorld::Width ||
            particle.position.y < 0.0 ||
            particle.position.y > SimulationWorld::Height) {
            throw std::invalid_argument(
                "Particle must be finite, have positive mass, and start inside the simulation world"
            );
        }

        particles.push_back(particle);
        invalidateNeighborList();
    }

    // Добавление частицы с проверкой минимального расстояния.
    bool tryAddParticle(const Particle& particle, double minDistance) {
        if (!particle.isValid() || !std::isfinite(minDistance) || minDistance < 0.0 ||
            particle.position.x < 0.0 || particle.position.x > SimulationWorld::Width ||
            particle.position.y < 0.0 || particle.position.y > SimulationWorld::Height) {
            return false;
        }

        const double minDistanceSquared = minDistance * minDistance;
        for (const Particle& existing : particles) {
            const double dx = existing.position.x - particle.position.x;
            const double dy = existing.position.y - particle.position.y;
            if (dx * dx + dy * dy < minDistanceSquared) {
                return false;
            }
        }

        particles.push_back(particle);
        invalidateNeighborList();
        return true;
    }

    std::size_t size() const {
        return particles.size();
    }

    double distanceBetween(const Particle& a, const Particle& b) const {
        return (b.position - a.position).magnitude();
    }

    // Кинетическая энергия.
    double calculateKineticEnergy() const {
        double energy = 0.0;
        for (const Particle& particle : particles) {
            const double speedSquared =
                particle.velocity.x * particle.velocity.x +
                particle.velocity.y * particle.velocity.y;
            energy += 0.5 * particle.mass * speedSquared;
        }
        return energy;
    }

    // Потенциальная энергия для force-shifted LJ с тем же cutoff, что и у силы.
    double calculatePotentialEnergy(double epsilon, double sigma) {
        validateLennardJonesParameters(epsilon, sigma);
        ensureNeighborList(sigma);

        if (forcesInitialized &&
            (cachedForceEpsilon != epsilon || cachedForceSigma != sigma)) {
            forcesInitialized = false;
        }

        double energy = 0.0;
        for (const Pair& pair : neighborPairs) {
            const Particle& a = particles[pair.first];
            const Particle& b = particles[pair.second];
            const double dx = b.position.x - a.position.x;
            const double dy = b.position.y - a.position.y;
            const double distanceSquared = dx * dx + dy * dy;

            if (distanceSquared < MinimumDistanceSquared) {
                continue;
            }

            energy += calculateForceShiftedLennardJonesPotential(
                distanceSquared, epsilon, sigma, CutoffFactor
            );
        }
        return energy;
    }

    // Средняя полная энергия на одну частицу.
    double calculateAverageEnergy(double epsilon, double sigma) {
        if (particles.empty()) {
            return 0.0;
        }

        const double totalEnergy =
            calculateKineticEnergy() + calculatePotentialEnergy(epsilon, sigma);
        return totalEnergy / static_cast<double>(particles.size());
    }

    // Расчёт ускорений по списку пар, построенному клеточной сеткой.
    void calculateForces(double epsilon, double sigma) {
        validateLennardJonesParameters(epsilon, sigma);
        ensureNeighborList(sigma);

        for (Particle& particle : particles) {
            particle.acceleration = {0.0, 0.0};
        }

        for (const Pair& pair : neighborPairs) {
            Particle& a = particles[pair.first];
            Particle& b = particles[pair.second];
            const double dx = b.position.x - a.position.x;
            const double dy = b.position.y - a.position.y;
            const double distanceSquared = dx * dx + dy * dy;

            if (distanceSquared < MinimumDistanceSquared) {
                continue;
            }

            const Vec2 force = calculateForceShiftedLennardJonesForceFromDelta(
                dx, dy, distanceSquared, epsilon, sigma, CutoffFactor
            );
            a.applyForce(force);
            b.applyForce({-force.x, -force.y});
        }

        cachedForceEpsilon = epsilon;
        cachedForceSigma = sigma;
        forcesInitialized = true;
    }

    // Интегратор Velocity Verlet.
    void step(double dt, double epsilon, double sigma) {
        if (!std::isfinite(dt) || dt <= 0.0) {
            throw std::invalid_argument("Simulation time step must be finite and positive");
        }
        validateLennardJonesParameters(epsilon, sigma);

        const bool neighborListRebuilt = ensureNeighborList(sigma);
        if (!forcesInitialized || neighborListRebuilt ||
            cachedForceEpsilon != epsilon || cachedForceSigma != sigma) {
            calculateForces(epsilon, sigma);
        }

        previousAccelerations.resize(particles.size());
        for (std::size_t i = 0; i < particles.size(); ++i) {
            previousAccelerations[i] = particles[i].acceleration;
        }

        // Предсказываем позиции и отражаем их от границ до расчёта новых сил.
        velocitySigns.resize(particles.size(), Vec2{1.0, 1.0});
        for (std::size_t i = 0; i < particles.size(); ++i) {
            velocitySigns[i] = {1.0, 1.0};
            Particle& particle = particles[i];
            particle.predictPosition(dt);
            if (reflectCoordinate(particle.position.x, SimulationWorld::Width)) {
                velocitySigns[i].x = -1.0;
            }
            if (reflectCoordinate(particle.position.y, SimulationWorld::Height)) {
                velocitySigns[i].y = -1.0;
            }
        }

        calculateForces(epsilon, sigma);

        for (std::size_t i = 0; i < particles.size(); ++i) {
            particles[i].correctVelocity(previousAccelerations[i], dt);
            particles[i].velocity.x *= velocitySigns[i].x;
            particles[i].velocity.y *= velocitySigns[i].y;
        }
    }

    // Обработка границ симуляции для внешних вызовов.
    bool handleWallCollisions(double width, double height) {
        if (!std::isfinite(width) || width <= 0.0 ||
            !std::isfinite(height) || height <= 0.0) {
            throw std::invalid_argument("Wall dimensions must be finite and positive");
        }

        bool collided = false;
        for (Particle& particle : particles) {
            if (!particle.isValid()) {
                throw std::invalid_argument("Cannot reflect an invalid particle state");
            }
            const bool flipX = reflectCoordinate(particle.position.x, width);
            const bool flipY = reflectCoordinate(particle.position.y, height);
            particle.velocity.x *= flipX ? -1.0 : 1.0;
            particle.velocity.y *= flipY ? -1.0 : 1.0;
            collided = collided || flipX || flipY;
        }

        if (collided) {
            forcesInitialized = false;
        }
        return collided;
    }

private:
    static constexpr double CutoffFactor = 2.5;
    static constexpr double SkinFactor = 0.3;
    static constexpr double MinimumDistanceSquared = 1e-18;

    std::vector<Pair> neighborPairs;
    std::vector<Vec2> neighborReferencePositions;
    std::vector<std::vector<std::size_t>> cells;
    std::size_t cellCountX = 0;
    std::size_t cellCountY = 0;
    double cellWidth = 0.0;
    double cellHeight = 0.0;
    double neighborListSigma = 0.0;
    double cachedForceEpsilon = 0.0;
    double cachedForceSigma = 0.0;
    bool neighborListValid = false;

    void invalidateNeighborList() {
        neighborListValid = false;
        forcesInitialized = false;
    }

    static std::size_t chooseCellCount(
        double extent,
        double listRadius,
        std::size_t upperBound
    ) {
        const double requested = std::floor(extent / listRadius);
        if (!std::isfinite(requested) || requested >= static_cast<double>(upperBound)) {
            return upperBound;
        }
        return std::max<std::size_t>(1, static_cast<std::size_t>(requested));
    }

    std::size_t cellIndex(double coordinate, double cellSize, std::size_t count) const {
        const double rawIndex = std::floor(coordinate / cellSize);
        if (rawIndex <= 0.0) {
            return 0;
        }
        if (rawIndex >= static_cast<double>(count)) {
            return count - 1;
        }
        return static_cast<std::size_t>(rawIndex);
    }

    bool ensureNeighborList(double sigma) {
        const double cutoff = CutoffFactor * sigma;
        const double skin = SkinFactor * sigma;
        const double listRadius = cutoff + skin;
        if (!std::isfinite(listRadius) || listRadius <= 0.0) {
            throw std::invalid_argument("Lennard-Jones neighbor radius must be finite and positive");
        }

        bool rebuild = !neighborListValid || neighborListSigma != sigma ||
                       neighborReferencePositions.size() != particles.size();
        if (!rebuild) {
            const double rebuildDistanceSquared = 0.25 * skin * skin;
            for (std::size_t i = 0; i < particles.size(); ++i) {
                if (!particles[i].isValid()) {
                    throw std::invalid_argument("Particle state became invalid during simulation");
                }
                const double dx = particles[i].position.x - neighborReferencePositions[i].x;
                const double dy = particles[i].position.y - neighborReferencePositions[i].y;
                if (dx * dx + dy * dy >= rebuildDistanceSquared) {
                    rebuild = true;
                    break;
                }
            }
        }
        if (!rebuild) {
            return false;
        }

        for (const Particle& particle : particles) {
            if (!particle.isValid() || particle.position.x < 0.0 ||
                particle.position.x > SimulationWorld::Width ||
                particle.position.y < 0.0 || particle.position.y > SimulationWorld::Height) {
                throw std::invalid_argument("Neighbor grid requires valid particles inside the simulation world");
            }
        }

        const std::size_t maxSize = std::numeric_limits<std::size_t>::max();
        const std::size_t maxCells = particles.size() > maxSize / 4
            ? maxSize
            : std::max<std::size_t>(1, particles.size() * 4);
        cellCountX = chooseCellCount(SimulationWorld::Width, listRadius, maxCells);
        const std::size_t maxRows = std::max<std::size_t>(1, maxCells / cellCountX);
        cellCountY = chooseCellCount(SimulationWorld::Height, listRadius, maxRows);
        cellWidth = SimulationWorld::Width / static_cast<double>(cellCountX);
        cellHeight = SimulationWorld::Height / static_cast<double>(cellCountY);

        const std::size_t requiredCellCount = cellCountX * cellCountY;
        if (cells.size() != requiredCellCount) {
            cells.resize(requiredCellCount);
        }
        for (auto& cell : cells) {
            cell.clear();
        }

        for (std::size_t i = 0; i < particles.size(); ++i) {
            const std::size_t x = cellIndex(particles[i].position.x, cellWidth, cellCountX);
            const std::size_t y = cellIndex(particles[i].position.y, cellHeight, cellCountY);
            cells[y * cellCountX + x].push_back(i);
        }

        neighborPairs.clear();
        const double listRadiusSquared = listRadius * listRadius;
        for (std::size_t i = 0; i < particles.size(); ++i) {
            const std::size_t cellX = cellIndex(particles[i].position.x, cellWidth, cellCountX);
            const std::size_t cellY = cellIndex(particles[i].position.y, cellHeight, cellCountY);
            const std::size_t minX = cellX == 0 ? 0 : cellX - 1;
            const std::size_t maxX = std::min(cellCountX - 1, cellX + 1);
            const std::size_t minY = cellY == 0 ? 0 : cellY - 1;
            const std::size_t maxY = std::min(cellCountY - 1, cellY + 1);

            for (std::size_t y = minY; y <= maxY; ++y) {
                for (std::size_t x = minX; x <= maxX; ++x) {
                    const auto& cell = cells[y * cellCountX + x];
                    for (const std::size_t j : cell) {
                        if (j <= i) {
                            continue;
                        }
                        const double dx = particles[j].position.x - particles[i].position.x;
                        const double dy = particles[j].position.y - particles[i].position.y;
                        if (dx * dx + dy * dy < listRadiusSquared) {
                            neighborPairs.push_back({i, j});
                        }
                    }
                }
            }
        }

        neighborReferencePositions.resize(particles.size());
        for (std::size_t i = 0; i < particles.size(); ++i) {
            neighborReferencePositions[i] = particles[i].position;
        }
        neighborListSigma = sigma;
        neighborListValid = true;
        forcesInitialized = false;
        return true;
    }

    static bool reflectCoordinate(double& position, double length) {
        if (!std::isfinite(position)) {
            throw std::overflow_error("Particle position became non-finite during wall reflection");
        }
        if (position >= 0.0 && position <= length) {
            return false;
        }

        double wrapped = std::fmod(position, 2.0 * length);
        if (wrapped < 0.0) {
            wrapped += 2.0 * length;
        }

        const bool flipVelocity = wrapped >= length;
        position = flipVelocity ? 2.0 * length - wrapped : wrapped;
        return flipVelocity;
    }
};
