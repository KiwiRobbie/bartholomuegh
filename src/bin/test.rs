use std::{borrow::Borrow, f32::consts::PI};

use bevy_math::Vec3;

#[derive(Clone, Copy, Debug)]
struct BhParameters {
    a: f32,
}

fn delta(r: f32, a: f32) -> f32 {
    r * r - 2.0 * r + a * a
}

fn P(r: f32, a: f32, b: f32) -> f32 {
    r * r + a * a - a * b
}

fn rho(r: f32, theta: f32, a: f32) -> f32 {
    (r * r + a * a * theta.cos().powi(2)).powf(0.5)
}
fn Theta(theta: f32, a: f32, b: f32, q: f32) -> f32 {
    q - theta.cos().powi(2) * (b * b / theta.sin().powi(2) - a * a)
}
fn R(a: f32, b: f32, q: f32, P: f32, Delta: f32) -> f32 {
    P * P - Delta * ((b - a).powi(2) + q)
}

fn omega(a: f32, r: f32, theta: f32) -> f32 {
    let sigma2 = (r * r + a * a).powi(2) - a * a * delta(r, a) * theta.sin().powi(2);
    2.0 * a * r / sigma2
}
fn omega_bar(r: f32, theta: f32, a: f32) -> f32 {
    let sigma = ((r * r + a * a).powi(2) - a * a * delta(r, a) * theta.sin().powi(2)).sqrt();
    sigma * theta.sin() / rho(r, theta, a)
}
fn alpha(r: f32, theta: f32, a: f32) -> f32 {
    let sigma = ((r * r + a * a).powi(2) - a * a * delta(r, a) * theta.sin().powi(2)).sqrt();
    rho(r, theta, a) * delta(r, a).sqrt() / sigma
}

#[derive(Debug)]
struct RayState {
    bh: BhParameters,

    r: f32,
    phi: f32,
    theta: f32,

    p_r: f32,
    p_theta: f32,

    rho: f32,
    P: f32,
    R: f32,
    Delta: f32,
    Theta: f32,

    b: f32, // Axial angular momentum p_phi
    q: f32, // Carter Constant
}
impl RayState {
    // Equations of motion
    fn dr(&self) -> f32 {
        self.Delta / (self.rho * self.rho) * self.p_r
    }
    fn dtheta(&self) -> f32 {
        self.p_theta / (self.rho * self.rho)
    }
    fn dphi(&self) -> f32 {
        (2.0 * (self.b - self.bh.a) - 2.0 * self.b * (self.theta.cos() / self.theta.sin()).powi(2))
            * 0.5
            / self.rho.powi(2)
    }
    fn dp_r(&self) -> f32 {
        let a = self.bh.a;
        let b = self.b;
        let r = self.r;
        let r2 = r * r;
        let cos2theta = self.theta.cos().powi(2);

        let c_1 = self.p_r.powi(2);
        let c_2 = self.p_theta.powi(2);

        let f = 2.0 * self.rho * self.rho;
        let f_prime = 4.0 * r;

        let term_1 = ((c_1 * (a * a + (r - 2.0) * r) - c_2) * f_prime) / f.powi(2)
            - (2.0 * c_1 * (r - 1.0)) / f;

        let c_3 = ((b - a).powi(2) + self.q);
        let c_4 = self.Theta;
        let q = self.P.powi(2) / self.Delta;
        let q_prime =
            2.0 * (self.P) * (a * a * (r + 1.0) + a * b * (r - 1.0) + (r - 3.0) * r * r) / (q * q);

        let term_2 = (f * q_prime - (-c_3 + c_4 + q) * f_prime) / f.powi(2);

        return term_1 + term_2;

        // -(
        //     (-self.bh.a.powi(2)+self.r.powi(2)+(self.r-1.0)*self.theta.cos().powi(2))/
        //     (self.bh.a.powi(2) + self.theta.cos().powi(2)+self.r.powi(2)).powi(2)
        // )*self.p_r.powi(2) // https://www.wolframalpha.com/input?i=Partial%5B%28r%5E2-2r%2Ba%5E2%29%2F%5B2%5Br%5E2%2Ba%5E2%2BCos%5B%CE%98%5D%5E2%5D%5D%2Cr%5D
        // +(
        //     self.r/(self.bh.a.powi(2)+self.theta.cos().powi(2) +self.r.powi(2)).powi(2)
        // ) * self.p_theta.powi(2) // https://www.wolframalpha.com/input?i=Partial%5B1%2F%5B2%5Br%5E2%2Ba%5E2%2BCos%5B%CE%98%5D%5E2%5D%5D%2Cr%5D
        // +(
        //     (a2 + a*b*r - a*b - r2) / ((a2 + r2 - 2.0*r).powi(2)*(a2 + cos2theta + r2)) - (r*(a2*self.Theta + a2 - a*b + r2*self.Theta + r2 - 2.0*r*self.Theta))/((a2 + r2 - 2.0*r)*(a2 + cos2theta + r2).powi(2))
        // ) //https://www.wolframalpha.com/input?i2d=true&i=Partial%5B%5C%2840%29Divide%5BPower%5Br%2C2%5D%2BPower%5Ba%2C2%5D-ab%2B%5C%2840%29Power%5Br%2C2%5D-2r%2BPower%5Ba%2C2%5D%5C%2841%29T%2C2%5C%2840%29Power%5Br%2C2%5D-2r%2BPower%5Ba%2C2%5D%5C%2841%29%5C%2840%29Power%5Br%2C2%5D%2BPower%5Ba%2C2%5D%2BPower%5BCos%5B%CE%B8%5D%2C2%5D%5C%2841%29%5D%5C%2841%29%2Cr%5D
    }
    fn dp_theta(&self) -> f32 {
        let denominator = 2.0 * self.Delta * self.rho.powi(2);

        let a = self.bh.a;
        let b = self.b;
        let b2 = b * b;
        let r = self.r;
        let cos_theta = self.theta.cos();
        let sin_theta = self.theta.sin();
        let cot_theta = cos_theta / sin_theta;

        let sin2_theta = sin_theta * sin_theta;
        let cos2_theta = cos_theta * cos_theta;

        let partial_reciprocal_rho2 = sin2_theta / (a * a + r * r + cos2_theta).powi(2);

        (
            -0.5*self.Delta*partial_reciprocal_rho2*self.p_r.powi(2)
        )//https://www.wolframalpha.com/input?i2d=true&i=Partial%5BDivide%5B1%2CPower%5Br%2C2%5D%2BPower%5Ba%2C2%5D%2BPower%5BCos%5B%CE%B8%5D%2C2%5D%5D%2C%CE%B8%5D
        +(
            -0.5*partial_reciprocal_rho2*self.p_theta.powi(2)
        )
        +(
            (2.0*self.Delta*(-a*a*sin_theta * cos_theta + b2*cot_theta.powi(3) + b2*cot_theta))/denominator
        ) // https://www.wolframalpha.com/input?i2d=true&i=Partial%5BDivide%5BR%2B%5C%2840%29T%5C%2841%29%5C%2840%29q-Power%5BCos%5B%CE%B8%5D%2C2%5D%5C%2840%29Divide%5BPower%5Bb%2C2%5D%2CPower%5BSin%5B%CE%B8%5D%2C2%5D%5D-Power%5Ba%2C2%5D%5C%2841%29%5C%2841%29%2CQ%5D%2C%CE%B8%5D
    }

    // Calculators
    fn omega(&self) -> f32 {
        omega(self.bh.a, self.r, self.theta)
    }

    fn omega_bar(&self) -> f32 {
        omega_bar(self.r, self.theta, self.bh.a)
    }

    fn alpha(&self) -> f32 {
        alpha(self.r, self.theta, self.bh.a)
    }

    pub fn new(
        bh: BhParameters,

        r: f32,
        theta: f32,
        phi: f32,

        p_r: f32,
        p_theta: f32,

        b: f32,
        q: f32,
    ) -> Self {
        let P = P(r, bh.a, b);
        let Delta = delta(r, bh.a);
        Self {
            bh,
            r,
            phi,
            theta,
            p_r,
            p_theta,
            b,
            q,
            rho: rho(r, theta, bh.a),
            Delta: Delta,
            Theta: Theta(theta, bh.a, b, q),
            P: P,
            R: R(bh.a, b, q, P, Delta),
        }
    }

    pub fn euler_step(&mut self, h: f32) {
        // Calculate derivatives
        let dr = self.dr();
        let dphi = self.dphi();
        let dtheta = self.dtheta();
        let dp_r = self.dp_r();
        let dp_theta = self.dp_theta();

        if dr > 100.0 || dphi > 100.0 || dtheta > 100.0 || dp_r > 100.0 || dp_theta > 100.0 {
            // dbg!(&self);
            // loop {}
        }

        // Update state values
        self.r += h * dr;
        self.phi += h * dphi;
        self.theta += h * dtheta;
        self.p_r += h * dp_r;
        self.p_theta += h * dp_theta;

        // Update calculated values
        self.rho = rho(self.r, self.theta, self.bh.a);
        self.P = P(self.r, self.bh.a, self.b);
        self.Delta = delta(self.r, self.bh.a);
        self.Theta = Theta(self.theta, self.bh.a, self.b, self.q);
        self.R = R(self.bh.a, self.b, self.q, self.P, self.Delta);
    }
}

struct BoyerLindquistObserver {
    r: f32,
    theta: f32,
    phi: f32,
    beta: f32,
    B_r: f32,
    B_theta: f32,
    B_phi: f32,
}
impl BoyerLindquistObserver {
    fn fido_ray(&self, n: Vec3, bh: BhParameters) -> RayState {
        let cartesian_fido = Vec3 {
            x: -n.x * (1.0 - self.beta.powi(2)).sqrt(),
            y: -n.y + self.beta,
            z: -n.z * (1.0 - self.beta.powi(2)).sqrt(),
        } / (1.0 - self.beta * n.y);

        let kappa = (1.0 - self.B_theta.powi(2)).sqrt();
        let spherical_fido = Vec3 {
            x: self.B_phi / kappa * cartesian_fido.x
                + self.B_r * cartesian_fido.y
                + self.B_r * self.B_theta / kappa * cartesian_fido.z, //r
            y: self.B_theta * cartesian_fido.y - kappa * cartesian_fido.z, //theta
            z: -self.B_r / kappa * cartesian_fido.x
                + self.B_phi * cartesian_fido.y
                + self.B_theta * self.B_theta / kappa * cartesian_fido.z, // phi
        };
        let mut ray = RayState::new(bh, self.r, self.theta, self.phi, 0.0, 0.0, 0.0, 0.0);

        let E_f = 1.0 / (ray.alpha() + ray.omega() * ray.omega_bar() * spherical_fido.z);

        ray.p_r = E_f * ray.rho / ray.Delta.sqrt() * spherical_fido.x;
        ray.p_theta = E_f * ray.rho * spherical_fido.y;
        ray.b = E_f * ray.omega_bar() * spherical_fido.z;

        ray.q = ray.p_theta * ray.p_theta
            + ray.theta.cos().powi(2)
                * (ray.b * ray.b / ray.theta.sin().powi(2) - ray.bh.a * ray.bh.a);
        return ray;
    }
}

const OUT_FILE_NAME: &'static str = "plots/3d-plot.svg";
use plotters::prelude::*;
fn main() {
    let bh = BhParameters { a: 0.0 };
    let r_c: f32 = 4.0;
    let theta_c: f32 = 0.5 * PI;
    let phi_c: f32 = 0.0;

    let Omega: f32 = 1.0 / (bh.a + r_c.powf(3.0 / 2.0));
    let observer = BoyerLindquistObserver {
        r: r_c,
        theta: theta_c,
        phi: phi_c,
        beta: omega_bar(r_c, theta_c, bh.a) / alpha(r_c, theta_c, bh.a)
            * (Omega - omega(bh.a, r_c, theta_c)),
        B_r: 0.0,
        B_theta: 0.0,
        B_phi: 1.0,
    };

    let area = SVGBackend::new(OUT_FILE_NAME, (1024, 760)).into_drawing_area();

    area.fill(&WHITE).unwrap();

    let x_axis = (-3.0..3.0).step(0.1);
    let z_axis = (-3.0..3.0).step(0.1);

    let mut chart = ChartBuilder::on(&area)
        .caption(format!("3D Plot Test"), ("sans", 20))
        .build_cartesian_3d(x_axis.clone(), -3.0..3.0, z_axis.clone())
        .unwrap();

    chart.with_projection(|mut pb| {
        pb.yaw = 0.5;
        pb.scale = 0.9;
        pb.into_matrix()
    });

    chart
        .configure_axes()
        .light_grid_style(BLACK.mix(0.15))
        .max_light_lines(3)
        .draw()
        .unwrap();

    for i in 0..10 {
        let n_theta = PI * (i as f32) / 5.0;

        let mut ray =
            observer.fido_ray(Vec3::new(n_theta.sin(), 0.0, n_theta.cos()).normalize(), bh);

        let mut data: Vec<(f64, f64, f64)> = vec![];
        for i in 0..10_000 {
            dbg!(&ray);
            ray.euler_step(0.001);

            if ray.r.is_nan()
                || ray.theta.is_nan()
                || ray.phi.is_nan()
                || ray.p_r.is_nan()
                || ray.p_theta.is_nan()
            {
                break;
            }

            let x = (ray.r * ray.r + ray.bh.a * ray.bh.a).sqrt() * ray.theta.sin() * ray.phi.cos();
            let y = (ray.r * ray.r + ray.bh.a * ray.bh.a).sqrt() * ray.theta.sin() * ray.phi.sin();
            let z = ray.r * ray.theta.cos();
            data.push((x as f64, y as f64, z as f64));
        }

        chart
            .draw_series(LineSeries::new(
                data.iter().map(|e| (e.0 as f64, e.1 as f64, e.2 as f64)),
                &BLACK,
            ))
            .unwrap();
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    area.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {}", OUT_FILE_NAME);
}
