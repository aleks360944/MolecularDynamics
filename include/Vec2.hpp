#pragma once
#include <cmath>

struct Vec2 {
    double x;
    double y;

    double magnitude() const {
        return std::sqrt(x * x + y * y);
    }

    Vec2 operator+(const Vec2& other) const {
        return {
            x + other.x,
            y + other.y
        };
    }

    Vec2 operator-(const Vec2& other) const {
        return {
            x - other.x,
            y - other.y
        };
    }
};
