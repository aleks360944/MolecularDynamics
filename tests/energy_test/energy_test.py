import numpy as np
import matplotlib.pyplot as plt

EPSILON = 1.0
SIGMA = 1.0
MASS = 1.0

DT = 0.001
STEPS = 10_000

def lennard_jones_potential(r):
    inv_r6 = (1.0 / r) ** 6
    return 4.0 * (inv_r6 ** 2 - inv_r6)

def lennard_jones_force(r):
    inv_r6 = (1.0 / r) ** 6
    return (24.0 / r) * (inv_r6 - 2.0 * inv_r6 ** 2)


def calculate_total_energy(x1, x2, v1, v2, mass):
    kinetic_energy = 0.5 * mass * v1**2 + 0.5 * mass * v2**2

    distance = abs(x2 - x1)
    potential_energy = lennard_jones_potential(distance)

    return kinetic_energy + potential_energy


# Начальные условия
x1 = 0.0
x2 = 1.5

v1 = 0.5
v2 = -0.5

x1_initial = x1
x2_initial = x2
v1_initial = v1
v2_initial = v2

mass = 1.0

time = 0.0

times = []
total_energies = []

r = abs(x2 - x1)

force = lennard_jones_force(r)

a1 = force / mass
a2 = -force / mass

initial_energy = calculate_total_energy(
    x1, x2, v1, v2, mass
)

for step in range(STEPS):

    # 1. Обновляем координаты
    x1 += v1 * DT + 0.5 * a1 * DT**2
    x2 += v2 * DT + 0.5 * a2 * DT**2

    # 2. Вычисляем новые силы
    r = abs(x2 - x1)

    force = lennard_jones_force(r)

    new_a1 = force / mass
    new_a2 = -force / mass

    # 3. Обновляем скорости
    v1 += 0.5 * (a1 + new_a1) * DT
    v2 += 0.5 * (a2 + new_a2) * DT

    # 4. Сохраняем новые ускорения
    a1 = new_a1
    a2 = new_a2

    # 5. Время и энергия
    time += DT

    if x1 >= x2:
        print("WARNING: Particles crossed!")
        break
    
    energy = calculate_total_energy(
        x1, x2, v1, v2, mass
    )

    times.append(time)
    total_energies.append(energy)
    
print("Initial energy:", initial_energy)
print("Final energy:", total_energies[-1])
print("Energy difference:", total_energies[-1] - initial_energy)
print("Distance:", abs(x2 - x1))
print("Potential energy:", lennard_jones_potential(abs(x2 - x1)))

plt.figure(figsize=(10, 5))

plt.plot(times, total_energies, label="Total Energy")

plt.axhline(
    y=initial_energy,
    linestyle="--",
    label="Initial Energy"
)

plt.xlabel("Time")
plt.ylabel("Energy")
plt.title("Energy Conservation — Velocity Verlet")
plt.legend()
plt.grid()

plt.show()