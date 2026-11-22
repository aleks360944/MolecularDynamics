def lennard_jones_potential(r: float, epsilon: float, sigma: float) -> float:
    inv_r6 = (sigma / r) ** 6
    return 4.0 * epsilon * (inv_r6 ** 2 - inv_r6)


def lennard_jones_force_magnitude(r: float, epsilon: float, sigma: float) -> float:
    inv_r6 = (sigma / r) ** 6
    # Возвращает скаляр: > 0 при отталкивании, < 0 при притяжении
    return (24.0 * epsilon / r) * (2.0 * inv_r6 ** 2 - inv_r6)


