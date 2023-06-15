#import bevy_core_pipeline::fullscreen_vertex_shader

struct MainPassUniforms {
    camera: mat4x4<f32>,
    camera_inverse: mat4x4<f32>,
    time: f32,
    surface_bool: u32,
    disk_mode: u32, 
    misc_bool: u32,
    step_count: i32,
    initial_step: f32,
    abs_error: f32,
    rel_error: f32,
    max_step: f32,
    method: u32,
    disk_start: f32,
    disk_end: f32,

    r: f32,
    theta: f32,
    phi: f32,
    a: f32,

    rho: f32,
    Delta: f32,
    Sigma: f32,
    alpha: f32,
    omega: f32,
    omega_bar: f32,

    B_r: f32,
    B_theta: f32,
    B_phi: f32,
    beta: f32,
};



@group(0) @binding(0)
var<uniform> uniforms: MainPassUniforms;

struct Ray {
    pos: vec3<f32>,
    dir: vec3<f32>,
};

fn get_camera(clip_space: vec2<f32>) -> Ray {
    let pos = uniforms.camera_inverse * vec4(clip_space.x, clip_space.y, 1.0, 1.0);
    let dir = uniforms.camera_inverse * vec4(clip_space.x, clip_space.y, 0.01, 1.0);
    let pos2 = pos.xyz / pos.w;
    let dir2 = normalize(dir.xyz / dir.w - pos2);
    return Ray(pos2, dir2);
}

#import "noise.wgsl"
#import "blackbody.wgsl"

fn fbm(v: vec3<f32>) -> f32 {
    return (0.5 * noise(1.0 * v) + 0.25 * noise(2.0 * v) + 0.125 * noise(4.0 * v) + 0.0625 * noise(8.0 * v));
}

fn skybox(dir: vec3<f32>) -> vec3<f32> {

    let color = BlackBodyRadiation(1000.0 + pow(10.0, 4.0 * noise(200.0 * dir)));
    let brightness = pow(
        noise(100.0 * dir) * noise(200.0 * dir) * noise(400.0 * dir),
        10.0
    );
    return 2.0 * log(color.a) * color.rgb * pow(1.0 * max(0.0, brightness - 0.001), 0.7);
    // return vec3(fbm(10.0 * dir)) * mix(vec3(0.75), vec3(0.25), checker(dir, 5.0));
}

fn checker(dir: vec3<f32>, frequency: f32) -> f32 {
    let r = fract(frequency * dir) - vec3(0.5);
    return sign(r.x * r.y * r.z) * 0.5 + 0.5;
}

fn surface(point: vec3<f32>) -> vec3<f32> {
    if uniforms.surface_bool == 1u {
        return vec3(fbm(10.0 * point)) * mix(vec3(0.8, 0.8, 0.8), vec3(0.2, 0.2, 0.8), checker(point, 1.0));
    } else {
        return vec3(0.0);
    }
}


fn disk_sample(point: vec3<f32>, t_offset: f32) -> vec3<f32> {
    let fract_t = fract(-uniforms.time * 0.05 + t_offset);
    let radius = dot(point, point);
    let cos_t = cos(1000.0 * fract_t / radius);
    let sin_t = sin(1000.0 * fract_t / radius);
    let rotation = mat3x3(cos_t, 0.0, - sin_t, 0.0, 0.0, 0.0, sin_t, 0.0, cos_t);
    let p = point * rotation ;


    let r = length(point) * 1.0 + 0.02 * radius * fbm(vec3(p + vec3(0.0, uniforms.time + fract_t / radius, 0.0)));
    let color = pow(1.0 - 1.0 / r, 2.0) * BlackBodyRadiation(1.0e4 * pow(sqrt(radius), -0.75)) / 255.0;
    let density = 10.0 * pow(fbm(vec3(r) + uniforms.time * 0.01 + 1000.0 * t_offset), 1.25) * fbm(p);
    return color.rgb * color.a * density * (1.0 - abs(1.0 - 2.0 * fract_t));
}

fn disk(point: vec3<f32>) -> vec3<f32> {
    if uniforms.disk_mode == 1u {
        return disk_sample(point, 0.0) + disk_sample(point, 1.0 / 3.0) + disk_sample(point, 2.0 / 3.0);
    } else if uniforms.disk_mode == 2u {

        return vec3(0.75) * mix(vec3(0.8, 0.2, 0.2), vec3(0.2, 0.2, 0.2), checker(point, 1.0));
    } else {
        return vec3(0.0);
    }
}

fn integrand(ray: Ray, h2: f32) -> Ray {
    return Ray(
        ray.dir,
        -1.5 * h2 * ray.pos / pow(dot(ray.pos, ray.pos), 2.5)
    );
}
fn rk4_step(ray: Ray, h: f32, h2: f32) -> Ray {
    let k1 = integrand(ray, h2);
    let k2 = integrand(Ray(ray.pos + 0.5 * h * k1.pos, ray.dir + 0.5 * h * k1.dir), h2);
    let k3 = integrand(Ray(ray.pos + 0.5 * h * k2.pos, ray.dir + 0.5 * h * k2.dir), h2);
    let k4 = integrand(Ray(ray.pos + h * k3.pos, ray.dir + h * k3.dir), h2);
    return Ray(
        ray.pos + (k1.pos + 2.0 * k2.pos + 2.0 * k3.pos + k4.pos) * h / 6.0,
        ray.dir + (k1.dir + 2.0 * k2.dir + 2.0 * k3.dir + k4.dir) * h / 6.0,
    );
}


fn metric_values(
    r: f32,
    theta: f32,
    phi: f32,
    a: f32,
    rho: ptr<function, f32>,
    Delta: ptr<function, f32>,
    Sigma: ptr<function, f32>,
    alpha: ptr<function, f32>,
    omega: ptr<function, f32>,
    omega_bar: ptr<function, f32>
) {
    *rho = sqrt(r * 2.0 + a * 2.0 * cos(theta));
    *Delta = r * r - 2.0 * r + a * a;
    *Sigma = sqrt(pow(r * r + a * a, 2.0) - a * a * *Delta * sin(theta) * sin(theta));
    *alpha = *rho * sqrt(*Delta) / *Sigma;
    *omega = 2.0 * a * r / (*Sigma * *Sigma);
    *omega_bar = *Sigma * sin(theta) / *rho;
}


fn ray_values(r: f32, theta: f32, phi: f32, a: f32, rho: f32, Delta: f32, Sigma: f32, alpha: f32, omega: f32, omega_bar: f32, b: f32, q: f32, P: ptr<function, f32>, R: ptr<function, f32>, Theta: ptr<function, f32>) {
    *P = r * r + a * a - a * b;
    *R = *P * *P - Delta * ((b - a) * (b - a) + q);
    *Theta = q - cos(theta) * cos(theta) * (b * b / (sin(theta) * sin(theta)) - a * a);
}


struct State {
    r: f32,
    theta: f32,
    phi: f32, 
    p_r: f32,
    p_theta: f32,
}

fn derivatives(
    r: f32,
    theta: f32,
    phi: f32,
    a: f32,
    rho: f32,
    Delta: f32,
    Sigma: f32,
    alpha: f32,
    omega: f32,
    omega_bar: f32,
    p_r: f32,
    p_theta: f32,
    b: f32,
    q: f32,
    P: f32,
    R: f32,
    Theta: f32,
) -> State {

    let x0_0: f32 = a * a;
    let x0_1: f32 = r * r;
    let o0: f32 = p_r * (-2.0 * r + x0_0 + x0_1) / (x0_0 * cos(theta) * cos(theta) + x0_1);

    let o1: f32 = p_theta / (a * a * cos(theta) * cos(theta) + r * r);

    let x2_0: f32 = a * a;
    let x2_1: f32 = 2.0 * r;
    let x2_2: f32 = r * r;
    let x2_3: f32 = cos(theta) * cos(theta);
    let x2_4: f32 = x2_0 * x2_3;
    let x2_5: f32 = a * x2_1;
    let o2: f32 = (-b * x2_1 + b * x2_2 + b * x2_4 - x2_3 * x2_5 + x2_5) / ((x2_2 + x2_4) * (x2_0 - x2_1 + x2_2) * sin(theta) * sin(theta));

    let x3_0: f32 = sin(theta) * sin(theta);
    let x3_1: f32 = 2.0 * r;
    let x3_2: f32 = r * r;
    let x3_3: f32 = a * a;
    let x3_4: f32 = x3_2 + x3_3;
    let x3_5: f32 = -x3_1 + x3_4;
    let x3_6: f32 = x3_5 * x3_5;
    let x3_7: f32 = cos(theta) * cos(theta);
    let x3_8: f32 = x3_2 + x3_3 * x3_7;
    let x3_9: f32 = p_r * p_r;
    let x3_10: f32 = x3_0 * x3_6;
    let x3_11: f32 = r - 1.0;
    let x3_12: f32 = -x3_11 * x3_8;
    let x3_13: f32 = q * x3_0 + x3_7 * (-b * b + x3_0 * x3_3);
    let x3_14: f32 = -a * b + x3_4;
    let x3_15: f32 = q + (a - b) * (a - b) ;
    let x3_16: f32 = x3_0 * (x3_14 * x3_14 - x3_15 * x3_5) + x3_13 * x3_5;
    let o3: f32 = (r * x3_10 * (p_theta * p_theta + x3_5 * x3_9) - r * x3_16 * x3_5 + x3_10 * x3_12 * x3_9 + x3_12 * x3_16 + x3_5 * x3_8 * (x3_0 * (x3_1 * x3_14 - x3_11 * x3_15) + x3_11 * x3_13)) / (x3_0 * x3_6 * x3_8 * x3_8);

    let x4_0: f32 = cos(theta);
    let x4_1: f32 = sin(theta);
    let x4_2: f32 = r * r;
    let x4_3: f32 = a * a;
    let x4_4: f32 = x4_2 + x4_3;
    let x4_5: f32 = -2.0 * r + x4_4;
    let x4_6: f32 = x4_0 * x4_0 ;
    let x4_7: f32 = x4_2 + x4_3 * x4_6;
    let x4_8: f32 = -b * b;
    let x4_9: f32 = pow(x4_1, 4.0) * x4_3;
    let x4_10: f32 = x4_1 * x4_1;
    let x4_11: f32 = x4_10 * x4_3;
    let o4: f32 = x4_0 * (x4_11 * (x4_10 * (-x4_5 * (q + (a - b) * (a - b)) + (-a * b + x4_4) * (-a * b + x4_4)) + x4_5 * (q * x4_10 + x4_6 * (x4_11 + x4_8))) + x4_5 * x4_7 * (-x4_8 - x4_9) + x4_5 * x4_9 * (-p_r * p_r * x4_5 - p_theta * p_theta)) / (x4_1 * x4_1 * x4_1 * x4_5 * x4_7 * x4_7);



    return State(
        o0,
        o1,
        o2,
        o3,
        o4
    );


    // // Prartial Derivatives
    // let drho_r = r / sqrt(a * a * cos(theta) + r * r);
    // let drho_theta = -(a * a) * sin(theta) / (2.0 * sqrt(a * a * cos(theta) + r * r));
    // let dDelta_r = 2.0 * r - 2.0;
    // let dP_b = -a;
    // let dR_b = -Delta * (-2.0 * a + 2.0 * b) + 2.0 * P * dP_b;
    // let dTheta_b = -2.0 * b * cos(theta) * cos(theta) / (sin(theta) * sin(theta));
    // let dTheta_theta = 2.0 * b * b * pow(cos(theta) / sin(theta), 3.0) + 2.0 * (-(a * a) + b * b / (sin(theta) * sin(theta))) * sin(theta) * cos(theta);

    // // State Derivatives
    // return State(
    //     Delta / rho * rho * p_r,
    //     1.0 / (rho * rho) * p_theta,
    //     (Delta * dTheta_b + dR_b) / (2.0 * Delta * rho * rho),
    //     (Theta * dDelta_r / (2.0 * Delta * rho * rho) + p_r * p_r * Delta * drho_r / (rho * rho * rho) - p_r * p_r * dDelta_r / (2.0 * rho * rho) + p_theta * p_theta * drho_r / (rho * rho * rho) - (R + Theta * Delta) * drho_r / (Delta * rho * rho * rho) - (R + Theta * Delta) * dDelta_r / (2.0 * Delta * Delta * rho * rho)),
    //     (Delta * p_r * p_r * drho_theta / (rho * rho * rho) + p_theta * p_theta * drho_theta / (rho * rho * rho) - (Delta * Theta + R) * drho_theta / (Delta * rho * rho * rho) + (Delta * dTheta_theta) / (2.0 * Delta * rho * rho)),
    // );
}


fn spherical_to_dir(theta: f32, phi: f32) -> vec3<f32> {
    return vec3(
        cos(phi) * cos(theta),
        sin(phi),
        cos(phi) * sin(theta)
    );
}

fn test(u: ptr<function, f32>) {
    *u = 1.0;
}


@fragment
fn fragment(in: FullscreenVertexOutput) -> @location(0) vec4<f32> {
    let clip_space = vec2(1.0, -1.0) * (in.uv * 2.0 - 1.0);
    var output_colour = vec3(0.0);

    var ray = get_camera(clip_space);
    




    // // Observer 
    var r: f32 = uniforms.r;
    var theta: f32 = uniforms.theta;
    var phi: f32 = uniforms.phi;
    var a: f32 = uniforms.disk_end / 100.0;
    var rho: f32 = uniforms.rho;
    var Delta: f32 = uniforms.Delta;
    var Sigma: f32 = uniforms.Sigma;
    var alpha: f32 = uniforms.alpha;
    var omega: f32 = uniforms.omega;
    var omega_bar: f32 = uniforms.omega_bar;
    let B_r: f32 = uniforms.B_r;
    let B_theta: f32 = uniforms.B_theta;
    let B_phi: f32 = uniforms.B_phi;
    let beta: f32 = uniforms.beta;


    // var r: f32 = 4.0;
    // var theta: f32 = 1.5707963267948966;
    // var phi: f32 = uniforms.time;
    // var a: f32 = uniforms.max_step;
    // var rho: f32 = 4.0;
    // var Delta: f32 = 8.0;
    // var Sigma: f32 = 16.0;
    // var alpha: f32 = 0.7071067811865476;
    // var omega: f32 = 0.0;
    // var omega_bar: f32 = 4.0;
    // let B_r: f32 = 0.0;
    // let B_theta: f32 = 0.0;
    // let B_phi: f32 = 1.0;
    // let beta: f32 = 0.7071067811865475;
    // (4.0, 1.5707963267948966, 0, 0.0, 4.0, 8.0, 16.0, 0.7071067811865476, 0.0, 4.0) (0.0, 0.0, 1.0) 0.7071067811865475



    // Camera ray
    let n: vec3<f32> = normalize(ray.dir);

    // Cartesian FIDO ray
    let quotient: f32 = 1.0 - beta * n.y;
    let nF: vec3<f32> = vec3(
        -n.x * sqrt(1.0 - beta * beta),
        (beta - n.y),
        -n.z * sqrt(1.0 - beta * beta)
    ) / quotient;


    // Spherical FIDO ray
    let kappa: f32 = sqrt(1.0 - B_theta * 2.0);
    let nF_r: f32 = (B_phi / kappa * nF.x + B_r * nF.y + B_r * B_theta / kappa * nF.z);
    let nF_theta: f32 = B_theta * nF.y - kappa * nF.z;
    let nF_phi: f32 = (-B_r / kappa * nF.x + B_phi * nF.y + B_theta * B_phi / kappa * nF.z);

    let E_f: f32 = 1.0 / (alpha + omega * omega_bar * nF_phi);
    let p_t: f32 = -1.0;
    var p_r: f32 = E_f * rho / sqrt(Delta) * nF_r;
    var p_theta: f32 = E_f * rho * nF_theta;
    let p_phi: f32 = E_f * omega_bar * nF_phi;

    // Constants of motion for photon
    let b = p_phi;
    let q = p_theta * p_theta + cos(theta) * cos(theta) * (b * b / (sin(theta) * sin(theta)) - a * a);

    // Calculate ray values 
    metric_values(
        r,
        theta,
        phi,
        a,
        &rho,
        &Delta,
        &Sigma,
        &alpha,
        &omega,
        &omega_bar,
    );

    var P: f32 = 0.0;
    var R: f32 = 0.0;
    var Theta: f32 = 0.0;
    ray_values(
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
        b,
        q,
        &P,
        &R,
        &Theta
    );
    let r_0 = r;
    // let h = -1.0;
    var h = uniforms.initial_step;
    let rel_tol = uniforms.rel_error;
    let abs_tol = uniforms.abs_error;
    let max_step = uniforms.max_step;


    let theta_0 = theta;

    var dstate: State;
    for (var i = 0; i < uniforms.step_count; i += 1) {
        let k1 = derivatives(
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
            p_r,
            p_theta,
            b,
            q,
            P,
            R,
            Theta
        );

        let coarse = State(
            r - h * k1.r,
            theta - h * k1.theta,
            phi - h * k1.phi,
            p_r - h * k1.p_r,
            p_theta - h * k1.p_theta,
        );


        r = r - 0.5 * h * k1.r;
        theta = theta - 0.5 * h * k1.theta;
        phi = phi - 0.5 * h * k1.phi;
        p_r = p_r - 0.5 * h * k1.p_r;
        p_theta = p_theta - 0.5 * h * k1.p_theta;

        metric_values(
            r,
            theta,
            phi,
            a,
            &rho,
            &Delta,
            &Sigma,
            &alpha,
            &omega,
            &omega_bar,
        );
        ray_values(
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
            b,
            q,
            &P,
            &R,
            &Theta
        );

        let k2 = derivatives(
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
            p_r,
            p_theta,
            b,
            q,
            P,
            R,
            Theta
        );

        r = r - 0.5 * h * k2.r;
        theta = theta - 0.5 * h * k2.theta;
        phi = phi - 0.5 * h * k2.phi;
        p_r = p_r - 0.5 * h * k2.p_r;
        p_theta = p_theta - 0.5 * h * k2.p_theta;

        let fine = State(
            r,
            theta,
            phi,
            p_r,
            p_theta,
        );

        metric_values(
            r,
            theta,
            phi,
            a,
            &rho,
            &Delta,
            &Sigma,
            &alpha,
            &omega,
            &omega_bar,
        );
        ray_values(
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
            b,
            q,
            &P,
            &R,
            &Theta
        );

        let error_state = State(
            coarse.r - fine.r,
            coarse.theta - fine.theta,
            coarse.phi - fine.phi,
            coarse.p_r - fine.p_r,
            coarse.p_r - fine.p_r,
        );

        let error = sqrt(
            error_state.r * error_state.r + error_state.theta * error_state.theta + error_state.phi * error_state.phi + error_state.p_r * error_state.p_r + error_state.p_r * error_state.p_r
        );

        let y = sqrt(dot(ray.pos, ray.pos) + dot(ray.dir, ray.dir));
        h = min(h * clamp(sqrt(max(abs_tol, abs(y) * rel_tol) / abs(error)), 0.3, 2.0), max_step);
    }

    // let h_vec = cross(ray.pos, ray.dir);
    // let h2 = dot(h_vec, h_vec);

 


    // var hit = 0.0;
    // var phi = 0.0;


    // let disk_start = uniforms.disk_start;
    // let disk_end = uniforms.disk_end;

    // var color = vec3(0.0);


    // for (var i = 0; i < uniforms.step_count; i += 1) {
    //     // Integrate using selected method
    //     var ray_coarse = Ray();
    //     var ray_fine = Ray();
    //     if uniforms.method == 0u {
    //         ray_coarse = rk4_step(ray, 2.0 * h, h2);
    //         let ray_mid = rk4_step(ray, h, h2);
    //         ray_fine = rk4_step(ray_mid, h, h2);
    //     } else {
    //         let k1 = integrand(ray, h2);
    //         ray_coarse = Ray(
    //             ray.pos + k1.pos * h,
    //             ray.dir + k1.dir * h,
    //         );
    //         let ray_mid = Ray(
    //             ray.pos + k1.pos * h * 0.5,
    //             ray.dir + k1.dir * h * 0.5,
    //         );
    //         let k2 = integrand(ray_mid, h2);
    //         ray_fine = Ray(
    //             ray_mid.pos + k2.pos * h * 0.5,
    //             ray_mid.dir + k2.dir * h * 0.5,
    //         );
    //     }

    //     // Determine error and adapt step size
    //     let error_ray = Ray(
    //         ray_coarse.pos - ray_fine.pos,
    //         ray_coarse.dir - ray_fine.dir
    //     );
    //     let error = sqrt(dot(error_ray.pos, error_ray.pos) + dot(error_ray.dir, error_ray.dir));
    //     let y = sqrt(dot(ray.pos, ray.pos) + dot(ray.dir, ray.dir));
    //     h = min(h * clamp(sqrt(max(abs_tol, abs(y) * rel_tol) / abs(error)), 0.3, 2.0), max_step);

    //     // Check for y-plane intersection
    //     if (ray_fine.pos.y * ray.pos.y) <= 0.0 {
    //         // Approximate hit location
    //         let t = -ray.pos.y / (ray_fine.pos.y - ray.pos.y);
    //         let hit_ray = Ray(
    //             mix(ray.pos, ray_fine.pos, t) * vec3(1.0, 0.0, 1.0),
    //             mix(ray.dir, ray_fine.dir, t)
    //         );

    //         // Check for disk intersection
    //         let hit_distance = dot(hit_ray.pos, hit_ray.pos);
    //         if hit_distance >= disk_start * disk_start && hit_distance < disk_end * disk_end {
    //             color += disk(hit_ray.pos);
    //         }
    //     }

    //     // If light ray passes below horizon
    //     if dot(ray_fine.pos, ray_fine.pos) <= 1.0 && dot(ray.pos, ray.pos) >= 1.0 {
    //         let l_i = sqrt(dot(ray.pos, ray.pos));
    //         let l_f = sqrt(dot(ray_fine.pos, ray_fine.pos));

    //         let t = - (l_i - 1.0) / (l_f - l_i);
    //         ray = Ray(
    //             mix(ray.pos, ray_fine.pos, t),
    //             mix(ray.dir, ray_fine.dir, t)
    //         );
    //         hit = 1.0; 
    //         break;
    //     }

    //     ray = ray_fine;
    // }

    // let r = ray.pos;
    // let v = ray.dir;

    // let r_hat = normalize(r);
    // let v_hat = normalize(v);


    // output_colour = max(mix(skybox(v_hat), surface(r_hat), hit), vec3(0.0));
    // output_colour += color;

    // output_colour = vec3(
    //     checker(n, 1.0),
    //     0.0,
    //     checker(nF, 1.0)
    // );


    // output_colour = vec3(
    //     checker(
    //         spherical_to_dir(
    //             p_theta,
    //             p_phi,
    //         ),
    //         5.0
    //     )
    // );

    // output_colour = skybox(
    //     spherical_to_dir(theta, phi)
    // );




    let _d = derivatives(
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
        p_r,
        p_theta,
        b,
        q,
        P,
        R,
        Theta
    );

    let r_hat = spherical_to_dir(theta, phi);
    let phi_hat = normalize(cross(vec3(0.0, 1.0, 0.0), r_hat)) ;
    let theta_hat = normalize(cross(phi_hat, r_hat)) ;



    let f = select(1.0, 0.0, length(r_hat) == 0.0);

    let dir = r_hat * _d.r + phi_hat * _d.phi / r + theta_hat * _d.theta / r;
    // output_colour = vec3(_d.r, _d.phi, _d.theta);
    output_colour = vec3(checker(
        // length(dir),
        spherical_to_dir(theta, phi),
        1.0
    ) * select(
        vec3(1.0, 0.0, 0.0),
        vec3(0.0, 0.0, 1.0),
        r > 2.0
    ));
    output_colour.y = vec3(checker(dir, 1.0)).y;
    return vec4<f32>(output_colour, 1.0);
}
