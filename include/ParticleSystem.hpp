#pragma once

#include <cmath>
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
    std::vector<Particle> particles;
    std::vector<Vec2> previousAccelerations;
    std::vector<Vec2> velocitySigns;
    bool forcesInitialized = false;

    void addParticle(const Particle& particle) {
        if (!particle.isValid() || particle.position.x < 0.0 ||
            particle.position.x > SimulationWorld::Width ||
            particle.position.y < 0.0 ||
            particle.position.y > SimulationWorld::Height) {
            throw std::invalid_argument("Particle must be finite, have positive mass, and start inside the simulation world");
        }
        particles.push_back(particle);
        forcesInitialized = false;
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

            const double distanceSquared = dx * dx + dy * dy;

            if (distanceSquared < minDistanceSquared) {
                return false;
            }
        }

        particles.push_back(particle);
        forcesInitialized = false;
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

    // Потенциальная энергия.
    double calculatePotentialEnergy(double epsilon, double sigma) const {
        validateLennardJonesParameters(epsilon, sigma);
        double energy = 0.0;

        for (std::size_t i = 0; i < particles.size(); ++i) {
            for (std::size_t j = i + 1; j < particles.size(); ++j) {
                const double dx = particles[j].position.x - particles[i].position.x;
                const double dy = particles[j].position.y - particles[i].position.y;
                const double distanceSquared = dx * dx + dy * dy;

                if (distanceSquared < 1e-18) {
                    continue;
                }

                const double ratio2 = (sigma * sigma) / distanceSquared;
                const double ratio6 = ratio2 * ratio2 * ratio2;
                const double ratio12 = ratio6 * ratio6;

                energy += 4.0 * epsilon * (ratio12 - ratio6);
            }
        }

        return energy;
    }

    // Средняя полная энергия на одну частицу.
    double calculateAverageEnergy(double epsilon, double sigma) const {
        if (particles.empty()) {
            return 0.0;
        }

        const double totalEnergy =
            calculateKineticEnergy() +
            calculatePotentialEnergy(epsilon, sigma);

        return totalEnergy / static_cast<double>(particles.size());
    }

    // Расчёт ускорений всех частиц.
    void calculateForces(double epsilon, double sigma) {
        validateLennardJonesParameters(epsilon, sigma);

        // Сбрасываем старые ускорения.
        for (Particle& particle : particles) {
            particle.acceleration = {0.0, 0.0};
        }

        // Рассчитываем взаимодействия между парами.
        for (std::size_t i = 0; i < particles.size(); ++i) {
            for (std::size_t j = i + 1; j < particles.size(); ++j) {

                Particle& a = particles[i];
                Particle& b = particles[j];

                const double dx = b.position.x - a.position.x;
                const double dy = b.position.y - a.position.y;
                const double distanceSquared = dx * dx + dy * dy;

                if (distanceSquared < 1e-18) {
                    continue;
                }

                const Vec2 force = calculateLennardJonesForceFromDelta(
                    dx,
                    dy,
                    distanceSquared,
                    epsilon,
                    sigma
                );

                a.applyForce(force);

                const Vec2 oppositeForce{
                    -force.x,
                    -force.y
                };

                b.applyForce(oppositeForce);
            }
        }
    }

    // Интегратор Velocity Verlet.
    void step(double dt, double epsilon, double sigma) {
        if (!std::isfinite(dt) || dt <= 0.0) {
            throw std::invalid_argument("Simulation time step must be finite and positive");
        }
        validateLennardJonesParameters(epsilon, sigma);

        if (!forcesInitialized) {
            // 1. Рассчитываем ускорения в текущих координатах.
            calculateForces(epsilon, sigma);
            forcesInitialized = true;
        }

        previousAccelerations.resize(particles.size());

        for (std::size_t i = 0; i < particles.size(); ++i) {
            previousAccelerations[i] = particles[i].acceleration;
        }

        // 2. Предсказываем новые позиции и отражаем их через границы.
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

        // 3. Рассчитываем новые ускорения.
        calculateForces(epsilon, sigma);

        // 4. Корректируем скорости.
        for (std::size_t i = 0; i < particles.size(); ++i) {
            particles[i].correctVelocity(
                previousAccelerations[i],
                dt
            );
            particles[i].velocity.x *= velocitySigns[i].x;
            particles[i].velocity.y *= velocitySigns[i].y;
        }
    }

    // Обработка границ симуляции.
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
