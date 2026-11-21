#include "Renderer.hpp"

#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>

#include <algorithm>
#include <array>
#include <cmath>

Renderer::Renderer(unsigned int width, unsigned int height, bool fullscreen)
    : window(
        sf::VideoMode({width, height}),
        "Molecular Dynamics",
        fullscreen ? sf::State::Fullscreen : sf::State::Windowed
      )
{
    window.setVerticalSyncEnabled(true);
    cameraView = window.getDefaultView();

    fontLoaded = font.openFromFile("C:/Windows/Fonts/arial.ttf");

    if (!fontLoaded) {
        std::cerr << "Warning: could not load font from C:/Windows/Fonts/arial.ttf\n";
    }
}

bool Renderer::isOpen() const {
    return window.isOpen();
}

std::optional<Renderer::ParticleSpawnDrag> Renderer::pollParticleCreation() {
    auto drag = pendingParticleDrag;
    pendingParticleDrag.reset();

    return drag;
}

sf::Vector2f Renderer::pixelToWorld(sf::Vector2i pixel) const {
    return window.mapPixelToCoords(pixel, cameraView);
}

void Renderer::processEvents() {
    while (const std::optional event = window.pollEvent()) {
        if (event->is<sf::Event::Closed>()) {
            window.close();
        }
        
        //Отдаление камеры
        if (const auto* wheel = event->getIf<sf::Event::MouseWheelScrolled>()) {
            if (wheel->wheel == sf::Mouse::Wheel::Vertical) {
                const float minZoom = 0.1f;
                const float maxZoom = 5.0f;
                const float zoomStep = 0.1f;

                const sf::Vector2f beforeCoord = window.mapPixelToCoords(wheel->position, cameraView);

                float factor = (wheel->delta > 0) ? (1.f - zoomStep) : (1.f + zoomStep);
                float nextZoom = std::clamp(currentZoom * factor, minZoom, maxZoom);

                cameraView.zoom(nextZoom / currentZoom);
                currentZoom = nextZoom;

                const sf::Vector2f afterCoord = window.mapPixelToCoords(wheel->position, cameraView);
                cameraView.move(beforeCoord - afterCoord);
            }
        }

        // 1. Нажали ЛКМ — запомнили начальную позицию
        if (const auto* mousePressed = event->getIf<sf::Event::MouseButtonPressed>()) {
            if (mousePressed->button == sf::Mouse::Button::Left) {
                dragStartPosition = mousePressed->position;
                dragCurrentPosition = mousePressed->position;
                isDragging = true;
            }
        }

        // 2. Двигаем мышь — обновляем текущую позицию
        if (const auto* mouseMoved = event->getIf<sf::Event::MouseMoved>()) {
            if (isDragging) {
                dragCurrentPosition = mouseMoved->position;
            }
        }

        // 3. Отпустили ЛКМ — фиксируем конечную позицию для main.cpp
        if (const auto* mouseReleased = event->getIf<sf::Event::MouseButtonReleased>()) {
            if (mouseReleased->button == sf::Mouse::Button::Left && isDragging) {
                dragCurrentPosition = mouseReleased->position;
                isDragging = false;
                pendingParticleDrag = ParticleSpawnDrag{
                    dragStartPosition,
                    dragCurrentPosition
                };
            }
        }
    }
}

sf::Color speedToColor(double speed) {
    constexpr double MAX_SPEED = 10.0;

    constexpr std::array<sf::Color, 5> colors = {
        sf::Color(45, 100, 255),  // Blue
        sf::Color(40, 190, 230),  // Cyan
        sf::Color(70, 210, 130),  // Green
        sf::Color(245, 210, 65),  // Yellow
        sf::Color(245, 90, 45)    // Red
    };

    double normalized = std::clamp(speed / MAX_SPEED, 0.0, 1.0);

    double position = normalized * (colors.size() - 1);

    std::size_t index = static_cast<std::size_t>(position);

    if (index >= colors.size() - 1) {
        return colors.back();
    }

    double t = position - index;

    const sf::Color& a = colors[index];
    const sf::Color& b = colors[index + 1];

    auto interpolate = [t](std::uint8_t start, std::uint8_t end) {
        return static_cast<std::uint8_t>(
            start + (end - start) * t
        );
    };

    return sf::Color(
        interpolate(a.r, b.r),
        interpolate(a.g, b.g),
        interpolate(a.b, b.b)
    );
}

void Renderer::render(const ParticleSystem& system, double averageEnergy) {
    window.clear(sf::Color(15, 17, 25));

    window.setView(cameraView);

    // Границы симуляции
    sf::RectangleShape walls({
        static_cast<float>(SimulationWorld::Width * SimulationWorld::PixelsPerUnit),
        static_cast<float>(SimulationWorld::Height * SimulationWorld::PixelsPerUnit)
    });

    walls.setPosition({0.f, 0.f});
    walls.setFillColor(sf::Color::Transparent);
    walls.setOutlineThickness(3.f);
    walls.setOutlineColor(sf::Color(65, 170, 220));

    window.draw(walls);

    for (const Particle& particle : system.particles) {
        sf::CircleShape circle(15.f);

        double speed = std::hypot(
            particle.velocity.x,
            particle.velocity.y
        );

        circle.setFillColor(speedToColor(speed));
        circle.setOrigin({15.f, 15.f});

        circle.setPosition({
            static_cast<float>(particle.position.x * SimulationWorld::PixelsPerUnit),
            static_cast<float>(particle.position.y * SimulationWorld::PixelsPerUnit)
        });

        window.draw(circle);
    }

    // Визуализация будущей частицы и её траектории
    if (isDragging) {
        sf::Vector2f startWorld =
            window.mapPixelToCoords(dragStartPosition, cameraView);

        sf::Vector2f currentWorld =
            window.mapPixelToCoords(dragCurrentPosition, cameraView);

        sf::CircleShape previewCircle(15.f);
        previewCircle.setFillColor(sf::Color(120, 170, 255, 100));
        previewCircle.setOrigin({15.f, 15.f});
        previewCircle.setPosition(startWorld);

        window.draw(previewCircle);

        sf::Vertex line[] = {
            sf::Vertex{startWorld, sf::Color(120, 170, 255)},
            sf::Vertex{currentWorld, sf::Color(255, 100, 100)}
        };

        window.draw(line, 2, sf::PrimitiveType::Lines);
    }

    // Индикатор средней энергии в левом нижнем углу.
    window.setView(window.getDefaultView());

    if (fontLoaded) {
        std::ostringstream energyTextStream;

        if (system.size() == 0) {
            energyTextStream << "Eavg: N/A";
        } else {
            energyTextStream << "Eavg: "
                             << std::fixed << std::setprecision(4)
                             << averageEnergy;
        }

        const float panelHeight = 34.0f;
        const float panelWidth = 190.0f;
        const float margin = 10.0f;
        const float panelY =
            static_cast<float>(window.getSize().y) - panelHeight - margin;

        sf::RectangleShape energyPanel({panelWidth, panelHeight});
        energyPanel.setPosition({margin, panelY});
        energyPanel.setFillColor(sf::Color(30, 34, 44, 220));

        sf::Text energyText(font, energyTextStream.str(), 18);
        energyText.setFillColor(sf::Color(225, 230, 240));
        energyText.setPosition({margin + 9.0f, panelY + 5.0f});

        window.draw(energyPanel);
        window.draw(energyText);
    }

    window.display();
}
