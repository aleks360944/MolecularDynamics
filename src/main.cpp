#include <SFML/Graphics.hpp>
#include "ParticleSystem.hpp"
#include "Renderer.hpp"
#include <iostream>


int main() {
    ParticleSystem system;

    constexpr double epsilon = 1.0;
    constexpr double sigma = 1.0;

    Renderer renderer(1920, 1080, true);

    sf::Clock clock;
    sf::Clock energyUpdateClock;

    constexpr double fixedDt = 0.001;
    double accumulator = 0.0;

    double averageEnergy = system.calculateAverageEnergy(epsilon, sigma);

    sf::Clock profilerClock;

    int totalSteps = 0;
    int totalFrames = 0;
    double totalPhysicsTime = 0.0;

    while (renderer.isOpen()) {
        renderer.processEvents();

        // 4. Получаем обе позиции в main.cpp после отпускания ЛКМ
        if (auto drag = renderer.pollParticleCreation()) {
            // Начальная позиция (где нажали ЛКМ)
            sf::Vector2f startWorld =
            renderer.pixelToWorld(drag->startPosition);

            sf::Vector2f endWorld =
                renderer.pixelToWorld(drag->endPosition);

            Vec2 startPos{
                startWorld.x / SimulationWorld::PixelsPerUnit,
                startWorld.y / SimulationWorld::PixelsPerUnit
            };

            Vec2 endPos{
                endWorld.x / SimulationWorld::PixelsPerUnit,
                endWorld.y / SimulationWorld::PixelsPerUnit
            };


            constexpr double K = 2.0; // коэфицент чувствительности

            Vec2 velocity{
                (endPos.x - startPos.x) * K,
                (endPos.y - startPos.y) * K
            };

            // Пока передаем startPos в позицию частицы, а endPos готова для расчетов
            Particle newParticle{
                startPos,
                velocity,
                {0.0, 0.0},
                1.0
            };
            
            // Минимальное расстояние примерно равно равновесному
            // расстоянию потенциала Леннарда — Джонса: 2^(1/6) * sigma.
            constexpr double minimumSpawnDistance = 1.122 * sigma;
            system.tryAddParticle(newParticle, minimumSpawnDistance);
        }

        double frameTime = clock.restart().asSeconds();

        // Защита от слишком большого скачка времени.
        if (frameTime > 0.25) {
            frameTime = 0.25;
        }

        accumulator += frameTime;

        sf::Clock physicsClock;

        int stepsThisFrame = 0;

        while (accumulator >= fixedDt) {
            system.step(fixedDt, epsilon, sigma);
            accumulator -= fixedDt;

            stepsThisFrame++;
        }

        totalSteps += stepsThisFrame;
        totalPhysicsTime += physicsClock.getElapsedTime().asMicroseconds();

        if (profilerClock.getElapsedTime().asSeconds() >= 1.0f) {
            double averageStepTime =
            totalSteps > 0
                ? totalPhysicsTime / totalSteps / 1000.0
                : 0.0;
                
            std::cout
            << "Average physics step: "
            << averageStepTime << " ms\n";
            

            std::cout
            << "FPS: " << totalFrames << '\n'
            << "Physics steps/sec: " << totalSteps << '\n'
            << "Average steps/frame: "
            << (totalFrames > 0
                ? static_cast<double>(totalSteps) / totalFrames
                : 0.0) << '\n'
            << "Physics time/frame: "
            << (totalFrames > 0
                ? totalPhysicsTime / totalFrames / 1000.0
                : 0.0) << " ms\n"
            << "-------------------------\n";

            totalSteps = 0;
            totalPhysicsTime = 0.0;

            profilerClock.restart();
        }

        totalFrames++;

        // Энергия пересчитывается 10 раз в секунду,
        // чтобы не выполнять дополнительный расчёт пар на каждом кадре.
        if (energyUpdateClock.getElapsedTime().asSeconds() >= 0.1f) {
            averageEnergy = system.calculateAverageEnergy(epsilon, sigma);
            energyUpdateClock.restart();
        }

        renderer.render(system, averageEnergy);
    }

    return 0;
}   
