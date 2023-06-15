import numpy as np
import matplotlib.pyplot as plt
import expressions
from math import sqrt, cos, sin, pi


class Observer:
    def __init__(self, r, theta, phi, a, B_r, B_theta, B_phi):
        self.metric_values = expressions.metric_values(r, theta, phi, a)
        (
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
        ) = self.metric_values

        self.B_r = B_r
        self.B_theta = B_theta
        self.B_phi = B_phi

        Omega = 1 / (a + r**1.5)
        self.beta = omega_bar / alpha * (Omega - omega)

        print(self.metric_values, (B_r, B_theta, B_phi), self.beta)

    def fido_ray(self, x, y, z):
        (
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
        ) = self.metric_values

        # Normalize
        L = sqrt(x * x + y * y + z * z)
        N_x = x / L
        N_y = y / L
        N_z = z / L

        # Cartesian FIDO ray
        _ = 1 - self.beta * N_y
        n_Fy = (self.beta - N_y) / _
        n_Fx = -N_x * sqrt(1 - self.beta**2) / _
        n_Fz = -N_z * sqrt(1 - self.beta**2) / _

        # Spherical FIDO ray
        kappa = sqrt(1 - self.B_theta**2)
        n_Fr = (
            self.B_phi / kappa * n_Fx
            + self.B_r * n_Fy
            + self.B_r * self.B_theta / kappa * n_Fz
        )
        n_Ftheta = self.B_theta * n_Fy - kappa * n_Fz
        n_Fphi = (
            -self.B_r / kappa * n_Fx
            + self.B_phi * n_Fy
            + self.B_theta * self.B_phi / kappa * n_Fz
        )

        # Determine conjugate momentum
        E_f = 1 / (alpha + omega * omega_bar * n_Fphi)
        p_t = -1
        p_r = E_f * rho / sqrt(Delta) * n_Fr
        p_theta = E_f * rho * n_Ftheta
        p_phi = E_f * omega_bar * n_Fphi

        # Constants of motion for photon
        b = p_phi
        q = p_theta**2 + cos(theta) ** 2 * (b**2 / sin(theta) ** 2 - a**2)

        return Ray(r, theta, phi, p_r, p_theta, b, q, a)


class Ray:
    def __init__(self, r, theta, phi, p_r, p_theta, b, q, a):
        self._calculate_metric_values(r, theta, phi, a)

        self.p_r = p_r
        self.p_theta = p_theta
        self.b = b
        self.q = q

        self._calculate__ray_values()

    def _calculate_metric_values(self, *args):
        (
            self.r,
            self.theta,
            self.phi,
            self.a,
            self.rho,
            self.Delta,
            self.Sigma,
            self.alpha,
            self.omega,
            self.omega_bar,
        ) = expressions.metric_values(*args)

    def _calculate__ray_values(self):
        (
            self.P,
            self.R,
            self.Theta,
        ) = expressions.ray_values(
            self.r,
            self.theta,
            self.phi,
            self.a,
            self.rho,
            self.Delta,
            self.Sigma,
            self.alpha,
            self.omega,
            self.omega_bar,
            self.b,
            self.q,
        )

    def euler_step(self, h):
        (
            d_r,
            d_theta,
            d_phi,
            d_p_r,
            d_p_theta,
        ) = expressions.derivatives(
            self.r,
            self.theta,
            self.phi,
            self.a,
            self.rho,
            self.Delta,
            self.Sigma,
            self.alpha,
            self.omega,
            self.omega_bar,
            self.p_r,
            self.p_theta,
            self.b,
            self.q,
            self.P,
            self.R,
            self.Theta,
        )

        self.r = self.r + h * d_r
        self.theta = self.theta + h * d_theta
        self.phi = self.phi + h * d_phi
        self.p_r = self.p_r + h * d_p_r
        self.p_theta = self.p_theta + h * d_p_theta

        if (
            d_r * d_r
            + d_theta * d_theta
            + d_phi * d_phi
            + d_p_r * d_p_r
            + d_p_theta * d_p_theta
        ) > 25.0**2:
            raise Exception("Tolerance Exceeded")

        self._update_values()

    def _update_values(self):
        self._calculate_metric_values(self.r, self.theta, self.phi, self.a)
        self._calculate__ray_values()


ax = plt.figure().add_subplot(projection="3d")


observer = Observer(100.0, pi / 2, 0, 0.0, 0.0, 0.0, 1.0)


for phi in np.linspace(-pi / 2, pi / 2, 2, endpoint=False):
    for theta in np.linspace(0, 2 * pi, 4, endpoint=False):
        print(theta, phi)
        x, y, z = (cos(theta) * cos(phi), sin(theta), cos(theta) * sin(phi))

        ray = observer.fido_ray(x, y, z)
        data = ([], [], [])
        for i in range(5000):
            try:
                ray.euler_step(0.025)
            except:
                break
            _r = ray.r
            _theta = ray.theta
            _phi = ray.phi

            data[0].append(sqrt(_r**2 + ray.a**2) * sin(_theta) * cos(_phi))
            data[1].append(sqrt(_r**2 + ray.a**2) * sin(_theta) * sin(_phi))
            data[2].append(_r * cos(_theta))

        ax.plot(
            *data,
            "r" if y > 0 else ("k" if y == 0 else "b"),
            label="parametric curve",
            linewidth=1,
            alpha=0.4
        )
# ax.plot(*data, label="parametric curve", linewidth=1, alpha=0.4)

R = 150
ax.axes.set_xlim3d(left=-R, right=R)
ax.axes.set_ylim3d(bottom=-R, top=R)
ax.axes.set_zlim3d(bottom=-R, top=R)

plt.show()
