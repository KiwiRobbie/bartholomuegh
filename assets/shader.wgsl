#import bevy_core_pipeline::fullscreen_vertex_shader

struct MainPassUniforms {
    camera: mat4x4<f32>,
    camera_inverse: mat4x4<f32>,
    time: f32,
    show_ray_steps: u32,
    indirect_lighting: u32,
    shadows: u32,
    misc_bool: u32,
    step_count: i32,
    initial_step: f32,
    abs_error: f32,
    rel_error: f32,
};

@group(0) @binding(0)
var<uniform> uniforms: MainPassUniforms;

struct Ray {
    pos: vec3<f32>,
    dir: vec3<f32>,
};

// returns the closest intersection and the furthest intersection
fn ray_box_dist(r: Ray, vmin: vec3<f32>, vmax: vec3<f32>) -> vec2<f32> {
    let v1 = (vmin.x - r.pos.x) / r.dir.x;
    let v2 = (vmax.x - r.pos.x) / r.dir.x;
    let v3 = (vmin.y - r.pos.y) / r.dir.y;
    let v4 = (vmax.y - r.pos.y) / r.dir.y;
    let v5 = (vmin.z - r.pos.z) / r.dir.z;
    let v6 = (vmax.z - r.pos.z) / r.dir.z;
    let v7 = max(max(min(v1, v2), min(v3, v4)), min(v5, v6));
    let v8 = min(min(max(v1, v2), max(v3, v4)), max(v5, v6));
    if v8 < 0.0 || v7 > v8 {
        return vec2(0.0);
    }

    return vec2(v7, v8);
}

fn get_camera(clip_space: vec2<f32>) -> Ray {
    let pos = uniforms.camera_inverse * vec4(clip_space.x, clip_space.y, 1.0, 1.0);
    let dir = uniforms.camera_inverse * vec4(clip_space.x, clip_space.y, 0.01, 1.0);
    let pos2 = pos.xyz / pos.w;
    let dir2 = normalize(dir.xyz / dir.w - pos2);
    return Ray(pos2, dir2);
}

#import "noise.wgsl"

fn fbm(v: vec3<f32>) -> f32 {
    return (0.5 * noise(1.0 * v) + 0.25 * noise(2.0 * v) + 0.125 * noise(4.0 * v));
}

fn skybox(dir: vec3<f32>) -> vec3<f32> {
    return vec3(fbm(10.0 * dir)) * mix(vec3(0.75), vec3(0.25), checker(dir, 5.0));
}

fn checker(dir: vec3<f32>, frequency: f32) -> f32 {
    let r = fract(frequency * dir) - vec3(0.5);
    return sign(r.x * r.y * r.z) * 0.5 + 0.5;
}

fn surface(point: vec3<f32>) -> vec3<f32> {
    return vec3(fbm(5.0 * point));
    // return vec3(fbm(10.0 * point)) * mix(vec3(0.8, 0.2, 0.2), vec3(0.2, 0.2, 0.8), checker(point, 5.0));
}

fn disk(point: vec3<f32>) -> vec3<f32> {
    return vec3(fbm(1.0 * point)) * mix(vec3(0.8, 0.2, 0.2), vec3(0.2, 0.2, 0.8), checker(point, 1.0));
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

@fragment
fn fragment(in: FullscreenVertexOutput) -> @location(0) vec4<f32> {
    let clip_space = vec2(1.0, -1.0) * (in.uv * 2.0 - 1.0);
    var output_colour = vec3(0.0);

    var ray = get_camera(clip_space);

    let h_vec = cross(ray.pos, ray.dir);
    let h2 = dot(h_vec, h_vec);

    var h = uniforms.initial_step;
    let rel_tol = uniforms.rel_error;
    let abs_tol = uniforms.abs_error;
    let max_step = 0.25;

    var hit = 0.0;
    var phi = 0.0;


    let disk_start = 2.0;
    let disk_end = 10.0;

    var disk_hit = 0.0;

    for (var i = 0; i < uniforms.step_count; i += 1) {
        let k1 = integrand(ray, h2);

        let ray_coarse = Ray(
            ray.pos + k1.pos * h,
            ray.dir + k1.dir * h,
        );
        let ray_mid = Ray(
            ray.pos + k1.pos * h * 0.5,
            ray.dir + k1.dir * h * 0.5,
        );
        let k2 = integrand(ray_mid, h2);

        let ray_fine = Ray(
            ray_mid.pos + k2.pos * h * 0.5,
            ray_mid.dir + k2.dir * h * 0.5,
        );


        // let ray_coarse = rk4_step(ray, 2.0 * h, h2);
        // let ray_fine = rk4_step(
        //     rk4_step(ray, h, h2),
        //     h,
        //     h2
        // );
        // ray = ray_coarse;
        // ray = ray_fine;

        let error_ray = Ray(
            ray_coarse.pos - ray_fine.pos,
            ray_coarse.dir - ray_fine.dir
        );
        let error = sqrt(dot(error_ray.pos, error_ray.pos) + dot(error_ray.dir, error_ray.dir));

        let y = sqrt(dot(ray.pos, ray.pos) + dot(ray.dir, ray.dir));
        h = min(h * clamp(sqrt(max(abs_tol, abs(y) * rel_tol) / abs(error)), 0.3, 2.0), max_step);

        if dot(ray_fine.pos, ray_fine.pos) < 1.0 && dot(ray.pos, ray.pos) > 1.0 {
            let l_i = sqrt(dot(ray.pos, ray.pos));
            let l_f = sqrt(dot(ray_fine.pos, ray_fine.pos));

            let t = -l_i / (l_f - l_i);
            ray = Ray(
                mix(ray.pos, ray_fine.pos, t),
                mix(ray.dir, ray_fine.dir, t)
            );
            hit = 1.0; 
            break;
        }
        if (ray_fine.pos.y * ray.pos.y) <= 0.0 {
            let total_dist = ray_fine.pos.y - ray.pos.y;
            let t = -ray.pos.y / total_dist;

            let hit_ray = Ray(
                mix(ray.pos, ray_fine.pos, t) * vec3(1.0, 0.0, 1.0),
                mix(ray.dir, ray_fine.dir, t)
            );

            let hit_distance = dot(hit_ray.pos, hit_ray.pos);
            if hit_distance > disk_start * disk_start && hit_distance < disk_end * disk_end {
                ray = hit_ray;
                disk_hit = 1.0;
            break;
            }
        }

        ray = ray_fine;
    }

    let r = ray.pos;
    let v = ray.dir;

    let r_hat = normalize(r);
    let v_hat = normalize(v);

    let radius = length(r);
    let early = max(0.0, (radius - 1.0) * (3.0 - radius));

    let warning_color = vec3(2.0, 0.5, 2.0);
    output_colour = max(mix(skybox(v_hat), surface(r_hat), hit), vec3(0.0));
    // output_colour = mix(output_colour, warning_color, early);
    output_colour = mix(output_colour, disk(r), disk_hit);
    return vec4<f32>(output_colour, 1.0);
}