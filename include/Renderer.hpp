#pragma once

#include <optional>
#include <SFML/Graphics.hpp>
#include "ParticleSystem.hpp"

class Renderer {
public:
    struct ParticleSpawnDrag {
        sf::Vector2i startPosition;
        sf::Vector2i endPosition;
    };

private:
    sf::RenderWindow window;
    sf::View cameraView;      // Добавлено: камера для зума
    float currentZoom = 1.0f; // Добавлено: текущий масштаб

    sf::Font font;
    bool fontLoaded = false;

    std::optional<ParticleSpawnDrag> pendingParticleDrag;
    sf::Vector2i dragStartPosition;
    sf::Vector2i dragCurrentPosition;
    bool isDragging = false;

public:
    Renderer(unsigned int width, unsigned int height, bool fullscreen = false);
    sf::Vector2f pixelToWorld(sf::Vector2i pixel) const;
    
    bool isOpen() const;
    void processEvents();

    std::optional<ParticleSpawnDrag> pollParticleCreation();

    void render(const ParticleSystem& system, double averageEnergy);
};
