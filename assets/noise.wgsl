fn mod289_(x: f32) -> f32 {
    var x1: f32;

    x1 = x;
    let e2: f32 = x1;
    let e3: f32 = x1;
    let e8: f32 = x1;
    return (e2 - (floor((e8 * (1.0 / 289.0))) * 289.0));
}

fn mod289_1(x2: vec4<f32>) -> vec4<f32> {
    var x3: vec4<f32>;

    x3 = x2;
    let e2: vec4<f32> = x3;
    let e3: vec4<f32> = x3;
    let e8: vec4<f32> = x3;
    return (e2 - (floor((e8 * (1.0 / 289.0))) * 289.0));
}

fn perm(x4: vec4<f32>) -> vec4<f32> {
    var x5: vec4<f32>;

    x5 = x4;
    let e2: vec4<f32> = x5;
    let e8: vec4<f32> = x5;
    let e10: vec4<f32> = x5;
    let e16: vec4<f32> = x5;
    let e18: vec4<f32> = mod289_1((((e10 * 34.0) + vec4<f32>(1.0)) * e16));
    return e18;
}

fn noise(p: vec3<f32>) -> f32 {
    var p1: vec3<f32>;
    var a: vec3<f32>;
    var d: vec3<f32>;
    var b: vec4<f32>;
    var k1_: vec4<f32>;
    var k2_: vec4<f32>;
    var c: vec4<f32>;
    var k3_: vec4<f32>;
    var k4_: vec4<f32>;
    var o1_: vec4<f32>;
    var o2_: vec4<f32>;
    var o3_: vec4<f32>;
    var o4_: vec2<f32>;

    p1 = p;
    let e3: vec3<f32> = p1;
    a = floor(e3);
    let e6: vec3<f32> = p1;
    let e7: vec3<f32> = a;
    d = (e6 - e7);
    let e10: vec3<f32> = d;
    let e11: vec3<f32> = d;
    let e15: vec3<f32> = d;
    d = ((e10 * e11) * (vec3<f32>(3.0) - (2.0 * e15)));
    let e20: vec3<f32> = a;
    b = (e20.xxyy + vec4<f32>(0.0, 1.0, 0.0, 1.0));
    let e29: vec4<f32> = b;
    let e31: vec4<f32> = b;
    let e33: vec4<f32> = perm(e31.xyxy);
    k1_ = e33;
    let e35: vec4<f32> = k1_;
    let e37: vec4<f32> = b;
    let e40: vec4<f32> = k1_;
    let e42: vec4<f32> = b;
    let e45: vec4<f32> = perm((e40.xyxy + e42.zzww));
    k2_ = e45;
    let e47: vec4<f32> = k2_;
    let e48: vec3<f32> = a;
    c = (e47 + e48.zzzz);
    let e53: vec4<f32> = c;
    let e54: vec4<f32> = perm(e53);
    k3_ = e54;
    let e56: vec4<f32> = c;
    let e60: vec4<f32> = c;
    let e64: vec4<f32> = perm((e60 + vec4<f32>(1.0)));
    k4_ = e64;
    let e66: vec4<f32> = k3_;
    let e71: vec4<f32> = k3_;
    o1_ = fract((e71 * (1.0 / 41.0)));
    let e78: vec4<f32> = k4_;
    let e83: vec4<f32> = k4_;
    o2_ = fract((e83 * (1.0 / 41.0)));
    let e90: vec4<f32> = o2_;
    let e91: vec3<f32> = d;
    let e94: vec4<f32> = o1_;
    let e96: vec3<f32> = d;
    o3_ = ((e90 * e91.z) + (e94 * (1.0 - e96.z)));
    let e102: vec4<f32> = o3_;
    let e104: vec3<f32> = d;
    let e107: vec4<f32> = o3_;
    let e110: vec3<f32> = d;
    o4_ = ((e102.yw * e104.x) + (e107.xz * (1.0 - e110.x)));
    let e116: vec2<f32> = o4_;
    let e118: vec3<f32> = d;
    let e121: vec2<f32> = o4_;
    let e124: vec3<f32> = d;
    return ((e116.y * e118.y) + (e121.x * (1.0 - e124.y)));
}