
#pragma once

#include <cmath>
#include "Vec2.hpp"

struct Particle {
    Vec2 position;
    Vec2 velocity;
    Vec2 acceleration;
    double mass;

    bool isValid() const {
        return std::isfinite(position.x) && std::isfinite(position.y) &&
               std::isfinite(velocity.x) && std::isfinite(velocity.y) &&
               std::isfinite(acceleration.x) && std::isfinite(acceleration.y) &&
               std::isfinite(mass) && mass > 0.0;
    }

    void applyForce(const Vec2& force) {
        acceleration.x += force.x / mass;
        acceleration.y += force.y / mass;
    }

    // Первый этап Velocity Verlet:
    // предсказание новой позиции.
    void predictPosition(double dt) {
        position.x += velocity.x * dt
                    + 0.5 * acceleration.x * dt * dt;

        position.y += velocity.y * dt
                    + 0.5 * acceleration.y * dt * dt;
    }

    // Второй этап Velocity Verlet:
    // коррекция скорости с использованием
    // старого и нового ускорения.
    void correctVelocity(const Vec2& previousAcceleration, double dt) {
        velocity.x += 0.5 * (previousAcceleration.x + acceleration.x) * dt;

        velocity.y += 0.5 * (previousAcceleration.y + acceleration.y) * dt;
    }
};
