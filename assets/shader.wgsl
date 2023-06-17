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


struct State {
    x: vec3<f32>,
    p: vec2<f32>,
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

    let x0: f32 = 2.0 * r;
    let x1: f32 = pow(r, 2.0);
    let x2: f32 = pow(a, 2.0);
    let x3: f32 = x1 + x2;
    let x4: f32 = -x0 + x3;
    let x5: f32 = cos(theta);
    let x6: f32 = pow(x5, 2.0);
    let x7: f32 = x2 * x6;
    let x8: f32 = x1 + x7;
    let x9: f32 = 1.0 / x8;
    let x10: f32 = sin(theta);
    let x11: f32 = pow(x10, 2.0);
    let x12: f32 = 1.0 / x11;
    let x13: f32 = 1.0 / x4;
    let x14: f32 = a * x0;
    let x15: f32 = pow(x4, 2.0);
    let x16: f32 = pow(x8, -2.0);
    let x17: f32 = pow(p_theta, 2.0);
    let x18: f32 = pow(p_r, 2.0);
    let x19: f32 = x11 * x15;
    let x20: f32 = r - 1.0;
    let x21: f32 = -x20 * x8;
    let x22: f32 = -pow(b, 2.0);
    let x23: f32 = x11 * x2;
    let x24: f32 = q * x11 + x6 * (x22 + x23);
    let x25: f32 = -a * b + x3;
    let x26: f32 = q + pow(a - b, 2.0);
    let x27: f32 = x4 * x8;
    let x28: f32 = x11 * (pow(x25, 2.0) - x26 * x4) + x24 * x4;
    let x29: f32 = pow(x10, 4.0) * x2;

    return State(
        -vec3(
            p_r * x4 * x9,
            p_theta * x9,
            x12 * x13 * x9 * (-b * x0 + b * x1 + b * x7 - x14 * x6 + x14),
        ),
        -vec2(
            x12 * x16 * (r * x19 * (x17 + x18 * x4) - r * x28 * x4 + x18 * x19 * x21 + x21 * x28 + x27 * (x11 * (x0 * x25 - x20 * x26) + x20 * x24)) / x15,
            x13 * x16 * x5 * (x23 * x28 + x27 * (-x22 - x29) + x29 * x4 * (-x17 - x18 * x4)) / pow(x10, 3.0),
        )
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
    var disk_color = vec3(0.0);

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

    // Camera ray (Nobody likes people who use different coordinates, this is honestly worse that spherical)
    var n: vec3<f32> = normalize(vec3(
        -ray.dir.z,
        -ray.dir.y,
        ray.dir.x,
    )) * rotationMatrix(vec3(0.0, 1.0, 0.0), -phi) * rotationMatrix(normalize(vec3(0.0, 0.0, 1.0)), theta + 0.5 * 3.141592653589793);

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


    var fine: State;
    var coarse: State;

    for (var i = 0; i < uniforms.step_count; i += 1) {
        let max_step = uniforms.max_step; // * state.x.x; // Hacky but works well

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


            // Check for disk intersection
            let hit_distance = hit_state.x.x;
            if hit_distance >= uniforms.disk_start * uniforms.disk_start && hit_distance < uniforms.disk_end * uniforms.disk_end {
                disk_color += disk(hit_state.x.x, hit_state.x.z);
            }
        }

        state = fine;
    }


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
    output_colour += disk_color;
    return vec4<f32>(output_colour, 1.0);
}
