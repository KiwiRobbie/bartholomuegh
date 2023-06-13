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

#import "noise.wgsl"

fn fbm(v: vec3<f32>) -> f32 {
    return (
           0.5*noise(1.0*v)
         +0.25*noise(2.0*v)
        +0.125*noise(4.0*v)
     
    );
}


fn skybox(dir: vec3<f32>) -> vec3<f32> {



    return vec3(fbm(10.0 * dir));
    // return round(10.0*dir)/10.0;
    
    // let v = fract(10.0*normalize(dir))*2.0 - 1.0;
    // return vec3(sign(v.x*v.y*v.z));
    // return vec3(pow(noise(100.0 * dir), 3.0));
}

fn surface(point: vec3<f32>) -> vec3<f32> {
    return vec3(fbm(point* 10.0), 0.0, fbm(vec3(1000.0) + point* 10.0));
    // return round(10.0*dir)/10.0;
    
    // let v = fract(10.0*normalize(dir))*2.0 - 1.0;
    // return vec3(sign(v.x*v.y*v.z));
    // return vec3(pow(noise(100.0 * dir), 3.0));
}


@fragment
fn fragment(in: FullscreenVertexOutput) -> @location(0) vec4<f32> {
    let clip_space = vec2(1.0, -1.0) * (in.uv * 2.0 - 1.0);
    var output_colour = vec3(0.0);

    var ray = get_camera(clip_space);

    var r = ray.pos;
    var v = ray.dir;


    let h_vec = cross(r, v);
    let h2 = dot(h_vec,h_vec);



    // var r_hat = normalize(r);
    // let normal = normalize(cross(r_hat, v));
    // var phi_hat = normalize(cross(normal, r_hat));

    // let phi_dash = dot(phi_hat, v);
    // var u = 1.0 / length(r);
    // var d_u = -1.0/sqrt(length(r))*dot(r_hat, v);
    // var d_u = 1.0 / (length(r) + dot(v, r_hat)) - u;

    let step_size =  uniforms.misc_float;
    let steps = 2000;

    var hit = 0.0; 
    var phi = 0.0; 
    for (var i = 0; i < steps; i += 1) {
        r += v * step_size;
        let dv = -1.5*h2*r / pow(dot(r,r), 2.5);
        v+=dv*step_size;

        if (dot(r,r) < 2.0) { 
            hit = 1.0; 
            break; 
        }
    }

    // r = normalize(r) * (1.0 / u);

    let r_hat = normalize(r);
    let v_hat = normalize(v);


    output_colour = max( mix( skybox( v_hat), surface(r_hat), hit ), vec3(0.0));
    return vec4<f32>(output_colour, 1.0);
}