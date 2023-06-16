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


fn disk_sample(r: f32, phi: f32, t_offset: f32) -> vec3<f32> {
    let fract_t = fract(-uniforms.time * 0.05 + t_offset);
    let radius = r * 100.0 / uniforms.disk_end + 2.0 ;
    let rotation = phi + 1000.0 * fract_t / radius;
    let p = vec3(sin(rotation), 0.0, cos(rotation)) * r;

    let r = r * 1.0 + 0.02 * radius * fbm(vec3(p + vec3(0.0, uniforms.time + fract_t / radius, 0.0)));
    let color = pow(1.0 - 1.0 / r, 2.0) * BlackBodyRadiation(1.55e4 * pow(sqrt(radius), -0.75)) / 255.0;
    let density = 10.0 * pow(fbm(vec3(r) + uniforms.time * 0.01 + 1000.0 * t_offset), 1.25) * fbm(p);
    return color.rgb * color.a * density * (1.0 - abs(1.0 - 2.0 * fract_t)) ;
}

fn disk(r: f32, theta: f32) -> vec3<f32> {
    if uniforms.disk_mode == 1u {
        return disk_sample(r, theta, 0.0) + disk_sample(r, theta, 1.0 / 3.0) + disk_sample(r, theta, 2.0 / 3.0);
    } else if uniforms.disk_mode == 2u {

        return vec3(0.75) * mix(vec3(0.8, 0.2, 0.2), vec3(0.2, 0.2, 0.2), checker(vec3(r, theta, 0.0), 1.0));
    } else {
        return vec3(0.0);
    }
}

// fn disk(point: vec3<f32>) -> vec3<f32> {
//     if uniforms.disk_mode == 1u {
//         return disk_sample(point, 0.0) + disk_sample(point, 1.0 / 3.0) + disk_sample(point, 2.0 / 3.0);
//     } else if uniforms.disk_mode == 2u {

//         return vec3(0.75) * mix(vec3(0.8, 0.2, 0.2), vec3(0.2, 0.2, 0.2), checker(point, 1.0));
//     } else {
//         return vec3(0.0);
//     }
// }


fn rk4_step(ray: State, h: f32, constants: Constants) -> State {
    let k1 = integrand(ray, constants);
    let k2 = integrand(State(ray.x + 0.5 * h * k1.x, ray.p + 0.5 * h * k1.p), constants);
    let k3 = integrand(State(ray.x + 0.5 * h * k2.x, ray.p + 0.5 * h * k2.p), constants);
    let k4 = integrand(State(ray.x + h * k3.x, ray.p + h * k3.p), constants);
    return State(
        ray.x + (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x) * h / 6.0,
        ray.p + (k1.p + 2.0 * k2.p + 2.0 * k3.p + k4.p) * h / 6.0,
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
    x: vec3<f32>,
    p: vec2<f32>,
    // r: f32,
    // theta: f32,
    // phi: f32, 
    // p_r: f32,
    // p_theta: f32,
}

struct Constants { 
    a: f32,
    b: f32,
    q: f32, 
}


const PI = 3.141592653589793;

fn integrand(
    state: State,
    constants: Constants,
) -> State {
    let r = state.x.x;
    let theta = state.x.y;
    let phi = state.x.z;
    let p_r = state.p.x;
    let p_theta = state.p.y;
    let b = constants.b;
    let q = constants.q;
    let a = constants.a;

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
        -vec3(o0, o1, o2),
        -vec2(o3, o4)
    );
}


fn spherical_to_dir(theta: f32, phi: f32) -> vec3<f32> {
    return vec3(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
    );
}

fn test(u: ptr<function, f32>) {
    *u = 1.0;
}


fn rotationMatrix(axis: vec3<f32>, angle: f32) -> mat3x3<f32> {
    let axis = normalize(axis);
    let  s = sin(angle);
    let c = cos(angle);
    let oc = 1.0 - c;

    return mat3x3(oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s, oc * axis.z * axis.x + axis.y * s, oc * axis.x * axis.y + axis.z * s, oc * axis.y * axis.y + c, oc * axis.y * axis.z - axis.x * s, oc * axis.z * axis.x - axis.y * s, oc * axis.y * axis.z + axis.x * s, oc * axis.z * axis.z + c);
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
    var a: f32 = uniforms.a;
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

    // Camera ray
    var n: vec3<f32> = normalize(vec3(
        -ray.dir.z,
        -ray.dir.y,
        ray.dir.x,
    ));
    // n = n * rotationMatrix(cross(vec3(0.0, 1.0, 0.0), ray.pos), -theta);
    // n = n * rotationMatrix(vec3(0.0, 1.0, 0.0), phi);
    n = n * rotationMatrix(vec3(0.0, 1.0, 0.0), -phi);
    n = n * rotationMatrix(normalize(vec3(0.0, 0.0, 1.0)), theta + 0.5 * 3.141592653589793);
    
    
    
    // n = vec3(    
    //     n.x * cos(theta),
    //     n.y * sin(theta),
    //     n.z * cos(theta) 
    // )


    // Cartesian FIDO ray
    let quotient: f32 = 1.0 - beta * n.y;
    let nF: vec3<f32> = vec3(
        -n.x * sqrt(1.0 - beta * beta),
        (beta - n.z),
        -n.y * sqrt(1.0 - beta * beta)
    ) / quotient;


    // Spherical FIDO ray
    let kappa: f32 = sqrt(1.0 - B_theta * B_theta);
    let nF_r: f32 = (B_phi / kappa * nF.x + B_r * nF.y + B_r * B_theta / kappa * nF.z);
    let nF_theta: f32 = B_theta * nF.y - kappa * nF.z;
    let nF_phi: f32 = (-B_r / kappa * nF.x + B_phi * nF.y + B_theta * B_phi / kappa * nF.z);

    let E_f: f32 = 1.0 / (alpha + omega * omega_bar * nF_phi);
    let p_t: f32 = -1.0;
    var p_r: f32 = E_f * rho / sqrt(Delta) * nF_r;
    var p_theta: f32 = E_f * rho * nF_theta;
    let p_phi: f32 = E_f * omega_bar * nF_phi;

    var state = State(vec3(r, theta, phi), vec2(p_r, p_theta));


    // Constants of motion for photon
    let b = p_phi;
    let q = p_theta * p_theta + cos(theta) * cos(theta) * (b * b / (sin(theta) * sin(theta)) - a * a);
    let constants = Constants(a, b, q);


    var h = 100.0 * uniforms.initial_step;
    let rel_tol = uniforms.rel_error;
    let abs_tol = uniforms.abs_error;



    var color = vec3(0.0);

    var fine: State;
    var coarse: State;

    for (var i = 0; i < uniforms.step_count; i += 1) {
        let max_step = uniforms.max_step * state.x.x; // Hacky but works well

        if uniforms.method == 0u {
            coarse = rk4_step(state, 2.0 * h, constants);
            let mid = rk4_step(state, h, constants);
            fine = rk4_step(mid, h, constants);
        } else {
            let k1 = integrand(
                state,
                constants
            );

            coarse = State(
                state.x - h * k1.x,
                state.p - h * k1.p,
            );

            let mid = State(
                state.x - 0.5 * h * k1.x,
                state.p - 0.5 * h * k1.p,
            );

            let k2 = integrand(
                mid,
                constants
            );

            fine = State(
                mid.x - 0.5 * h * k2.x,
                mid.p - 0.5 * h * k2.p,
            );
        }



        let error_state = State(
            fine.x - coarse.x,
            fine.p - coarse.p
        );

        let error = sqrt(
            dot(error_state.x, error_state.x) + dot(error_state.p, error_state.p)
        );

        let y = sqrt(
            dot(state.x, state.x) + dot(state.p, state.p)
        );
        h = min(h * clamp(sqrt(max(abs_tol, abs(y) * rel_tol) / abs(error)), 0.3, 2.0), max_step);


        let ct_a = cos(state.x.y) ;
        let cb_b = cos(fine.x.y);

        // Check for y-plane intersection
        if (ct_a * cb_b) <= 0.0 {

            // Approximate hit location
            let t = -ct_a / (cb_b - ct_a);
            let hit_state = State(
                mix(state.x, fine.x, t),
                mix(state.p, fine.p, t)
            );
            // hit_state = state;
            let disk_start = uniforms.disk_start;
            let disk_end = uniforms.disk_end;


            // Check for disk intersection
            let hit_distance = hit_state.x.x;
            if hit_distance >= disk_start * disk_start && hit_distance < disk_end * disk_end {
                color += disk(hit_state.x.x, hit_state.x.z);
            }
        }

        state = fine;
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




    let _d = integrand(
        state,
        constants
    );

    r = state.x.x;
    theta = state.x.y;
    phi = state.x.z;

    let v1 = vec3(sin(theta) * cos(phi), sin(phi) * sin(theta), cos(theta)) * _d.x.x;
    let v2 = vec3(r * cos(phi) * cos(theta), r * sin(phi) * cos(theta), -r * sin(theta)) * _d.x.y;
    let v3 = vec3(-r * sin(phi) * sin(theta), r * sin(theta) * cos(phi), 0.0) * _d.x.z;
    let dir = normalize(v1 + v2 + v3);

    output_colour = vec3(checker(
        dir,
        5.0
    ) * select(
        vec3(1.0, 0.0, 0.0),
        vec3(0.0, 0.0, 1.0),
        r > 50.0
    ));

    // output_colour = dir;
    // output_colour.y = vec3(checker(dir, 5.0)).y;
    output_colour = skybox(dir);
    // output_colour = vec3(select(0.0, 1.0, state.x.z > 0.0));
    output_colour += color;
    return vec4<f32>(output_colour, 1.0);
}
