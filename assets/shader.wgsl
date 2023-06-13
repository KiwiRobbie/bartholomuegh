#import bevy_core_pipeline::fullscreen_vertex_shader

struct MainPassUniforms {
    camera: mat4x4<f32>,
    camera_inverse: mat4x4<f32>,
    time: f32,
    show_ray_steps: u32,
    indirect_lighting: u32,
    shadows: u32,
    misc_bool: u32,
    misc_float: f32,
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
    if (v8 < 0.0 || v7 > v8) {
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

@fragment
fn fragment(in: FullscreenVertexOutput) -> @location(0) vec4<f32> {
    let clip_space = vec2(1.0, -1.0) * (in.uv * 2.0 - 1.0);
    var output_colour = vec3(0.0);

    var ray = get_camera(clip_space);

    var r = ray.pos;
    var v = ray.dir;

    let r_hat = normalize(r);
    let phi_hat = normalize(cross(cross(r_hat, v), r_hat));

    // let r_dash = dot(r_hat, v);
    let phi_dash = dot(phi_hat, v);

    var u = 1.0 / length(r);
    var u_dash = u - 1.0 / (length(r) + dot(v, r_hat));

    let step_size = uniforms.misc_float;
    let steps = 1000;
    for (var i = 0; i < steps; i += 1) {
        let u_double_dash = -u * (1.0 - 1.5 * u * u);
        u_dash += u_double_dash * step_size;
        u += u_dash * step_size;

        let phi_hat = normalize(cross(cross(r_hat, v), r_hat));
        r += phi_hat * phi_dash * step_size;

    //     v = r_dash * r_hat + phi_dash * phi_hat;
    //     r += v * step_size;
    }

    // let dir = 

    r = normalize(r) * (1.0 / u);
    output_colour = r;

    // output_colour = r + ray_box_dist(ray, vec3(-0.5), vec3(0.5)).x;

    // output_colour = vec3(r_dash, phi_dash, 0.0);
    // output_colour = r_dash * r_hat;

    output_colour = max(output_colour, vec3(0.0));
    return vec4<f32>(output_colour, 1.0);
}