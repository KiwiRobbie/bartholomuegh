from math import sqrt, cos, sin

# Ray = (r, theta, phi, a, rho, Delta, Sigma, alpha, omega, omega_bar, p_r, p_theta, b, q, P, R, Theta)
# Metric = (r, theta, phi, a, rho, Delta, Sigma, alpha, omega, omega_bar)


def metric_values(r, theta, phi, a):
    rho = sqrt(r**2 + a**2 * cos(theta))
    Delta = r**2 - 2 * r + a**2
    Sigma = sqrt((r**2 + a**2) ** 2 - a**2 * Delta * sin(theta) ** 2)
    alpha = rho * sqrt(Delta) / Sigma
    omega = 2 * a * r / Sigma**2
    omega_bar = Sigma * sin(theta) / rho
    return (r, theta, phi, a, rho, Delta, Sigma, alpha, omega, omega_bar)


def ray_values(r, theta, phi, a, rho, Delta, Sigma, alpha, omega, omega_bar, b, q):
    P = r**2 + a**2 - a * b
    R = P**2 - Delta * ((b - a) ** 2 + q)
    Theta = q - cos(theta) ** 2 * (b**2 / sin(theta) ** 2 - a**2)
    return (P, R, Theta)


def derivatives(
    # Metric
    r,
    theta,
    phi,
    a,
    rho,
    Delta,
    Sigma,
    alpha,
    omega,
    omega_bar,
    # Ray
    p_r,
    p_theta,
    b,
    q,
    P,
    R,
    Theta,
):
    # Partial derivatives of relevant quantities
    drho_r = r / sqrt(a**2 * cos(theta) + r**2)
    drho_theta = -(a**2) * sin(theta) / (2 * sqrt(a**2 * cos(theta) + r**2))
    dDelta_r = 2 * r - 2
    dP_b = -a
    dR_b = -Delta * (-2 * a + 2 * b) + 2 * P * dP_b
    dTheta_b = -2 * b * cos(theta) ** 2 / sin(theta) ** 2
    dTheta_theta = 2 * b**2 * cos(theta) ** 3 / sin(theta) ** 3 + 2 * (
        -(a**2) + b**2 / sin(theta) ** 2
    ) * sin(theta) * cos(theta)

    # Derivatives
    d_r = Delta / rho**2 * p_r

    d_theta = 1 / rho**2 * p_theta

    d_phi = (Delta * dTheta_b + dR_b) / (2 * Delta * rho**2)
    d_p_r = (
        Theta * dDelta_r / (2 * Delta * rho**2)
        + p_r**2 * Delta * drho_r / rho**3
        - p_r**2 * dDelta_r / (2 * rho**2)
        + p_theta**2 * drho_r / rho**3
        - (R + Theta * Delta) * drho_r / (Delta * rho**3)
        - (R + Theta * Delta) * dDelta_r / (2 * Delta**2 * rho**2)
    )

    d_p_theta = (
        Delta * p_r**2 * drho_theta / rho**3
        + p_theta**2 * drho_theta / rho**3
        - (Delta * Theta + R) * drho_theta / (Delta * rho**3)
        + (Delta * dTheta_theta) / (2 * Delta * rho**2)
    )

    return (d_r, d_theta, d_phi, d_p_r, d_p_theta)
