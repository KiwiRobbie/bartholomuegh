import sympy as sp
from sympy.parsing.sympy_parser import parse_expr as expr
import numpy as np
import matplotlib.pyplot as plt


# Expressions
rho = expr("sqrt(r**2 + a**2*cos(theta))")
Delta = expr("r**2 -2*r+a**2")
Sigma = expr("sqrt((r**2+a**2)**2 -a**2*Delta*sin(theta)**2)")
alpha = expr("rho*sqrt(Delta)/Sigma")
omega = expr("2* a*r/Sigma**2")
omega_bar = expr("Sigma*sin(theta)/rho")

P = expr("r**2 + a**2 - a*b")
R = expr("P**2 - Delta*((b-a)**2+q)")
Theta = expr("q-cos(theta)**2* (b**2/sin(theta)**2-a**2)")

d_r = expr("Delta/rho**2*p_r")
d_theta = expr("1/rho**2*p_theta")

d_phi = expr("(Delta*Derivative(Theta(b), b) + Derivative(R(b), b))/(2*Delta*rho**2)")
d_p_r = expr(
    "Theta*Derivative(Delta(r), r)/(2*Delta*rho**2) + p_r**2*Delta*Derivative(rho(r), r)/rho**3 - p_r**2*Derivative(Delta(r), r)/(2*rho**2) + p_theta**2*Derivative(rho(r), r)/rho**3 - (R + Theta*Delta)*Derivative(rho(r), r)/(Delta*rho**3) - (R + Theta*Delta)*Derivative(Delta(r), r)/(2*Delta**2*rho**2)"
)
d_p_theta = expr(
    "Delta*p_r**2*Derivative(rho(theta), theta)/rho**3 + p_theta**2*Derivative(rho(theta), theta)/rho**3 - (Delta*Theta + R)*Derivative(rho(theta), theta)/(Delta*rho**3) + (Delta*Derivative(Theta(theta), theta) + Theta*Derivative(Delta(theta), theta) + Derivative(R(theta), theta))/(2*Delta*rho**2)"
)

# d_phi = expr("(R(b)+Delta*Theta(b))/(2*Delta*rho**2)").diff("b")
# d_p_r = sp.diff(
#     expr(
#         "-Delta(r)/(2*rho(r)**2)*p_r**2-1/(2*rho(r)**2)*p_theta**2+((R+Delta(r)*Theta)/(2*Delta(r)*rho(r)**2))"
#     ),
#     "r",
# )
# d_p_theta = sp.diff(
#     expr(
#         "-Delta/(2*rho(theta)**2)*p_r**2-1/(2*rho(theta)**2)*p_theta**2+((R(theta)+Delta(theta)*Theta(theta))/(2*Delta*rho(theta)**2))"
#     ),
#     "theta",
# )


def metric_values(r, theta, phi, a):
    rho_v = rho.evalf(subs={"r": r, "a": a, "theta": theta})
    Delta_v = Delta.evalf(subs={"r": r, "a": a})
    Sigma_v = Sigma.evalf(subs={"r": r, "a": a, "theta": theta, "Delta": Delta_v})
    alpha_v = alpha.evalf(subs={"rho": rho_v, "Delta": Delta_v, "Sigma": Sigma_v})
    omega_v = omega.evalf(subs={"a": a, "r": r, "Sigma": Sigma_v})
    omega_bar_v = omega_bar.evalf(
        subs={"a": a, "r": r, "Sigma": Sigma_v, "theta": theta, "rho": rho_v}
    )

    return {
        "r": r,
        "theta": theta,
        "phi": phi,
        "a": a,
        "rho": rho_v,
        "Delta": Delta_v,
        "Sigma": Sigma_v,
        "alpha": alpha_v,
        "omega": omega_v,
        "omega_bar": omega_bar_v,
    }


def sub_metric_values(expr, values):
    return expr.evalf(subs=values)


class Observer:
    def __init__(self, r, theta, phi, a, B_r, B_theta, B_phi):
        self.r = r
        self.theta = theta
        self.phi = phi
        self.a = a
        self.metric_values = metric_values(r, theta, phi, a)

        self.B_r = B_r
        self.B_theta = B_theta
        self.B_phi = B_phi

        Omega = expr("1/(a+r**1.5)")
        beta = omega_bar / alpha * (Omega - omega)
        self.beta = self.sub_values(beta)

    def sub_values(self, expr):
        return sub_metric_values(expr, self.metric_values)

    def fido_ray(self, x, y, z):
        # Normalize
        L = sp.sqrt(x * x + y * y + z * z)
        N_x = (x / L).evalf()
        N_y = (y / L).evalf()
        N_z = (z / L).evalf()

        # Cartesian FIDO ray
        _ = self.sub_values(1 - self.beta * N_y)
        n_Fy = self.sub_values((self.beta - N_y) / _)
        n_Fx = self.sub_values(-N_x * sp.sqrt(1 - self.beta**2) / _)
        n_Fz = self.sub_values(-N_z * sp.sqrt(1 - self.beta**2) / _)

        # Spherical FIDO ray
        kappa = sp.sqrt(1 - self.B_theta**2)
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
        E_f = self.sub_values(1 / (expr("alpha") + expr("omega*omega_bar") * n_Fphi))
        p_t = -1
        p_r = self.sub_values(E_f * rho / sp.sqrt(Delta) * n_Fr)
        p_theta = self.sub_values(E_f * rho * n_Ftheta)
        p_phi = self.sub_values(E_f * omega_bar * n_Fphi)

        # Constants of motion for photon
        b = p_phi
        q = self.sub_values(
            p_theta**2
            + sp.cos(self.theta) ** 2 * (b**2 / sp.sin(self.theta) ** 2 - self.a**2)
        )

        return Ray(self.r, self.theta, self.phi, p_r, p_theta, b, q, self.a)

        # # Partial derivative stuff:
        # def _partials(sym_b, sym_r=None, sym_theta=None):
        #     if sym_r is None and sym_theta is None:
        #         sym_r = sym_b
        #         sym_theta = sym_b

        #     return (sym_b.diff("b"), sym_r.diff("r"), sym_theta.diff("theta"))

        # partials_rho = _partials(rho)
        # partials_Delta = _partials(Delta)
        # partials_Sigma = _partials(
        #     Sigma,
        #     Sigma.subs("Delta", "Delta(r)"),
        #     Sigma,
        # )

        # partials_alpha = _partials(
        #     alpha,
        #     alpha.subs({"rho": "rho(r)", "Delta": "Delta(r)", "Sigma": "Sigma(r)"}),
        #     alpha.subs({"rho": "rho(theta)", "Sigma": "Sigma(theta)"}),
        # )

        # partials_P = _partials(P)
        # partials_R = _partials(
        #     R.subs({"P": "P(b)"}),
        #     R.subs({"Delta": "Delta(r)"}),
        #     R,
        # )
        # partials_Theta = _partials(Theta)


# partial_derivative_expressions = {
#     "rho": (
#         expr("0"),
#         expr("r / sqrt(a**2 * cos(theta) + r**2)"),
#         expr("-(a**2) * sin(theta) / (2 * sqrt(a**2 * cos(theta) + r**2))"),
#     ),
#     "Delta": (expr("0"), expr("2 * r - 2"), expr("0")),
#     "Sigma": (
#         expr("0"),
#         expr(
#             """
#             (
#                 -(a**2) * sin(theta) ** 2 * Derivative(Delta(r), r) / 2
#                 + 2 * r * (a**2 + r**2)
#             )
#             / sqrt(-(a**2) * Delta(r) * sin(theta) ** 2 + (a**2 + r**2) ** 2)
#             """
#         ),
#         expr(
#             """
#             -Delta
#             * a**2
#             * sin(theta)
#             * cos(theta)
#             / sqrt(-Delta * a**2 * sin(theta) ** 2 + (a**2 + r**2) ** 2)
#             """
#         ),
#     ),
#     "alpha": (
#         expr("0"),
#         expr(
#             """
#             sqrt(Delta(r)) * Derivative(rho(r), r) / Sigma(r)
#             - sqrt(Delta(r)) * rho(r) * Derivative(Sigma(r), r) / Sigma(r) ** 2
#             + rho(r) * Derivative(Delta(r), r) / (2 * sqrt(Delta(r)) * Sigma(r))
#             """
#         ),
#         expr(
#             """
#             sqrt(Delta) * Derivative(rho(theta), theta) / Sigma(theta)
#             - sqrt(Delta)
#             * rho(theta)
#             * Derivative(Sigma(theta), theta)
#             / Sigma(theta) ** 2
#             """
#         ),
#     ),
#     "P": (expr("-a"), expr("2 * r"), expr("0")),
#     "R": (
#         expr("-Delta * (-2 * a + 2 * b) + 2 * P(b) * Derivative(P(b), b)"),
#         expr("-(q + (-a + b) ** 2) * Derivative(Delta(r), r)"),
#         expr("0"),
#     ),
#     "Theta": (
#         expr("-2 * b * cos(theta) ** 2 / sin(theta) ** 2"),
#         expr("0"),
#         expr(
#             """
#             2 * b**2 * cos(theta) ** 3 / sin(theta) ** 3
#             + 2 * (-(a**2) + b**2 / sin(theta) ** 2) * sin(theta) * cos(theta)
#             """
#         ),
#     ),
# }


class Ray:
    def __init__(self, r, theta, phi, p_r, p_theta, b, q, a):
        self.values = metric_values(r, theta, phi, a)

        self.values["p_r"] = p_r
        self.values["p_theta"] = p_theta
        self.values["b"] = b
        self.values["q"] = q

        self._calculate_values()

    def _calculate_values(self):
        self.values["P"] = self.eval(P)
        self.values["R"] = self.eval(R)
        self.values["Theta"] = self.eval(Theta)

    def eval(self, expression):
        return expression.evalf(subs=self.values)

    def _calculate_partials(self):
        partials = {}
        partials["rho"] = (
            self.eval(expr("0")),
            self.eval(expr("r / sqrt(a**2 * cos(theta) + r**2)")),
            self.eval(
                expr("-(a**2) * sin(theta) / (2 * sqrt(a**2 * cos(theta) + r**2))")
            ),
        )

        partials["Delta"] = (
            self.eval(expr("0")),
            self.eval(expr("2 * r - 2")),
            self.eval(expr("0")),
        )
        partials["Sigma"] = (
            0.0,
            self.eval(
                expr(
                    """ ( -(a**2) * sin(theta) ** 2 * Derivative(Delta(r), r) / 2 + 2 * r * (a**2 + r**2) ) / sqrt(-(a**2) * Delta * sin(theta) ** 2 + (a**2 + r**2) ** 2) """
                ).subs(expr("Derivative(Delta(r), r)"), partials["Delta"][1])
            ),
            self.eval(
                expr(
                    """ -Delta * a**2 * sin(theta) * cos(theta) / sqrt(-Delta * a**2 * sin(theta) ** 2 + (a**2 + r**2) ** 2) """
                )
            ),
        )
        partials["alpha"] = (
            0.0,
            self.eval(
                expr(
                    """ sqrt(Delta) * Derivative(rho(r), r) / Sigma - sqrt(Delta) * rho * Derivative(Sigma(r), r) / Sigma ** 2 + rho * Derivative(Delta(r), r) / (2 * sqrt(Delta) * Sigma) """
                ).subs(
                    {
                        expr("Derivative(rho(r), r)"): partials["rho"][1],
                        expr("Derivative(Sigma(r), r)"): partials["Sigma"][1],
                        expr("Derivative(Delta(r), r)"): partials["Delta"][1],
                    }
                )
            ),
            self.eval(
                expr(
                    """ sqrt(Delta) * Derivative(rho(theta), theta) / Sigma - sqrt(Delta) * rho * Derivative(Sigma(theta), theta) / Sigma ** 2 """
                ).subs(
                    {
                        expr("Derivative(rho(theta), theta)"): partials["rho"][2],
                        expr("Derivative(Sigma(theta), theta)"): partials["Sigma"][2],
                    }
                )
            ),
        )
        partials["P"] = (-self.a, 2 * self.r, 0.0)
        partials["R"] = (
            self.eval(
                expr("-Delta * (-2 * a + 2 * b) + 2 * P * Derivative(P(b), b)").subs(
                    {
                        expr("Derivative(P(b), b)"): partials["P"][0],
                    }
                )
            ),
            self.eval(
                expr("-(q + (-a + b) ** 2) * Derivative(Delta(r), r)").subs(
                    {
                        expr("Derivative(Delta(r), r)"): partials["Delta"][1],
                    }
                )
            ),
            0.0,
        )
        partials["Theta"] = (
            self.eval(expr("-2 * b * cos(theta) ** 2 / sin(theta) ** 2")),
            0.0,
            self.eval(
                expr(
                    """ 2 * b**2 * cos(theta) ** 3 / sin(theta) ** 3 + 2 * (-(a**2) + b**2 / sin(theta) ** 2) * sin(theta) * cos(theta) """
                )
            ),
        )
        self.partial_values = partials
        self.partials_b = dict(
            (expr(f"Derivative({symbol}(b),b)"), expressions[0])
            for (symbol, expressions) in partials.items()
        )
        self.partials_r = dict(
            (expr(f"Derivative({symbol}(r),r)"), expressions[1])
            for (symbol, expressions) in partials.items()
        )
        self.partials_theta = dict(
            (expr(f"Derivative({symbol}(theta),theta)"), expressions[2])
            for (symbol, expressions) in partials.items()
        )

    def euler_step(self, h):
        self._calculate_partials()

        deltas = (
            h * self.eval(d_r),
            h * self.eval(d_theta),
            h * self.eval(d_phi.subs(self.partials_b)),
            h * self.eval(d_p_r.subs(self.partials_r)),
            h * self.eval(d_p_theta.subs(self.partials_theta)),
        )

        # print(deltas)

        # for a, b, c in zip(
        #     self.partials_b.items(),
        #     self.partials_r.items(),
        #     self.partials_theta.items(),
        # ):
        #     print(a)
        #     print(b)
        #     print(c)

        # print(self)
        # quit()
        self.r = (self.r + deltas[0]).evalf()
        self.theta = (self.theta + deltas[1]).evalf()
        self.phi = (self.phi + deltas[2]).evalf()
        self.p_r = (self.p_r + deltas[3]).evalf()
        self.p_theta = (self.p_theta + deltas[4]).evalf()
        self._update_values()

    def __str__(self):
        string = f"=== RAY {id(self)} ===\n"
        for key, val in self.values.items():
            string += f"\t{key}: {val}\n"
        return string

    def _update_values(self):
        self.values.update(metric_values(self.r, self.theta, self.phi, self.a))
        self._calculate_values()

    @property
    def r(self):
        return self.values["r"]

    @r.setter
    def r(self, value):
        self.values["r"] = value

    @property
    def theta(self):
        return self.values["theta"]

    @theta.setter
    def theta(self, value):
        self.values["theta"] = value

    @property
    def phi(self):
        return self.values["phi"]

    @phi.setter
    def phi(self, value):
        self.values["phi"] = value

    @property
    def p_r(self):
        return self.values["p_r"]

    @p_r.setter
    def p_r(self, value):
        self.values["p_r"] = value

    @property
    def p_theta(self):
        return self.values["p_theta"]

    @p_theta.setter
    def p_theta(self, value):
        self.values["p_theta"] = value

    @property
    def a(self):
        return self.values["a"]


observer = Observer(4, sp.pi / 2, 0, 0.999, 0.0, 0.0, 1.0)
ray = observer.fido_ray(1, 1, 1)

data = ([], [], [])

for i in range(100):
    print(i)
    ray.euler_step(0.01)
    _r = ray.r
    _theta = ray.theta
    _phi = ray.phi
    data[0].append(
        (sp.sqrt(_r**2 + ray.a**2) * sp.sin(_theta) * sp.cos(_phi)).evalf()
    )
    data[1].append(
        (sp.sqrt(_r**2 + ray.a**2) * sp.sin(_theta) * sp.sin(_phi)).evalf()
    )
    data[2].append((_r * sp.cos(_theta)).evalf())


# ax = plt.figure().add_subplot(projection="3d")

# x = data[0]
# y = data[1]
# z = data[2]

# ax.plot(x, y, z, label="parametric curve")
# eps = 1e-16
# ax.axes.set_xlim3d(left=-5.0 - eps, right=5 + eps)
# ax.axes.set_ylim3d(bottom=-5.0 - eps, top=5 + eps)
# ax.axes.set_zlim3d(bottom=-5.0 - eps, top=5 + eps)
# ax.legend()

# plt.show()
