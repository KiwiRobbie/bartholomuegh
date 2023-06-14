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
    let max_step = uniforms.max_step;


    var hit = 0.0;
    var phi = 0.0;


    let disk_start = uniforms.disk_start;
    let disk_end = uniforms.disk_end;

    var color = vec3(0.0);


    for (var i = 0; i < uniforms.step_count; i += 1) {
        // Integrate using selected method
        var ray_coarse = Ray();
        var ray_fine = Ray();
        if uniforms.method == 0u {
            ray_coarse = rk4_step(ray, 2.0 * h, h2);
            let ray_mid = rk4_step(ray, h, h2);
            ray_fine = rk4_step(ray_mid, h, h2);
        } else {
            let k1 = integrand(ray, h2);
            ray_coarse = Ray(
                ray.pos + k1.pos * h,
                ray.dir + k1.dir * h,
            );
            let ray_mid = Ray(
                ray.pos + k1.pos * h * 0.5,
                ray.dir + k1.dir * h * 0.5,
            );
            let k2 = integrand(ray_mid, h2);
            ray_fine = Ray(
                ray_mid.pos + k2.pos * h * 0.5,
                ray_mid.dir + k2.dir * h * 0.5,
            );
        }

        // Determine error and adapt step size
        let error_ray = Ray(
            ray_coarse.pos - ray_fine.pos,
            ray_coarse.dir - ray_fine.dir
        );
        let error = sqrt(dot(error_ray.pos, error_ray.pos) + dot(error_ray.dir, error_ray.dir));
        let y = sqrt(dot(ray.pos, ray.pos) + dot(ray.dir, ray.dir));
        h = min(h * clamp(sqrt(max(abs_tol, abs(y) * rel_tol) / abs(error)), 0.3, 2.0), max_step);

        // Check for y-plane intersection
        if (ray_fine.pos.y * ray.pos.y) <= 0.0 {
            // Approximate hit location
            let t = -ray.pos.y / (ray_fine.pos.y - ray.pos.y);
            let hit_ray = Ray(
                mix(ray.pos, ray_fine.pos, t) * vec3(1.0, 0.0, 1.0),
                mix(ray.dir, ray_fine.dir, t)
            );

            // Check for disk intersection
            let hit_distance = dot(hit_ray.pos, hit_ray.pos);
            if hit_distance >= disk_start * disk_start && hit_distance < disk_end * disk_end {
                color += disk(hit_ray.pos);
            }
        }

        // If light ray passes below horizon
        if dot(ray_fine.pos, ray_fine.pos) <= 1.0 && dot(ray.pos, ray.pos) >= 1.0 {
            let l_i = sqrt(dot(ray.pos, ray.pos));
            let l_f = sqrt(dot(ray_fine.pos, ray_fine.pos));

            let t = - (l_i - 1.0) / (l_f - l_i);
            ray = Ray(
                mix(ray.pos, ray_fine.pos, t),
                mix(ray.dir, ray_fine.dir, t)
            );
            hit = 1.0; 
            break;
        }

        ray = ray_fine;
    }

    let r = ray.pos;
    let v = ray.dir;

    let r_hat = normalize(r);
    let v_hat = normalize(v);


    output_colour = max(mix(skybox(v_hat), surface(r_hat), hit), vec3(0.0));
    output_colour += color;
    return vec4<f32>(output_colour, 1.0);
}