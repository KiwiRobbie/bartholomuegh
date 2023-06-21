use bevy::{
    input::mouse::{MouseMotion, MouseWheel},
    prelude::*,
};

const SENSITIVITY: f32 = 0.004;

#[derive(Component)]
pub struct CharacterEntity {
    pub speed: f32,
    pub velocity: Vec3,
    pub grounded: bool,
    pub look_at: Vec3,
    pub up: Vec3,
}

impl Default for CharacterEntity {
    fn default() -> Self {
        Self {
            speed: 10.0,
            velocity: Vec3::ZERO,
            grounded: false,
            look_at: Vec3::Z,
            up: Vec3::Y,
        }
    }
}

pub struct CharacterPlugin;

impl Plugin for CharacterPlugin {
    fn build(&self, app: &mut App) {
        app.add_system(update_character);
    }
}

fn update_character(
    mut character: Query<(&mut Transform, &mut CharacterEntity)>,
    keys: Res<Input<KeyCode>>,
    mouse: Res<Input<MouseButton>>,
    mut mouse_motion_events: EventReader<MouseMotion>,
    time: Res<Time>,
    mut mouse_wheel_events: EventReader<MouseWheel>,
) {
    let (mut transform, mut character) = character.single_mut();
    // speed
    for event in mouse_wheel_events.iter() {
        character.speed *= 1.0 + event.y.min(50.0) / 100.0;
    }

    // rotation
    if mouse.pressed(MouseButton::Left) {
        let mut mouse_delta = Vec2::new(0.0, 0.0);
        for event in mouse_motion_events.iter() {
            mouse_delta += event.delta;
        }
        if mouse_delta != Vec2::ZERO {
            let angle = character.look_at.dot(character.up).acos();
            let max_angle = 0.01;

            // Order is important to prevent unintended roll
            character.look_at = Quat::from_axis_angle(Vec3::Y, -mouse_delta.x * SENSITIVITY)
                * Quat::from_axis_angle(
                    transform.local_x(),
                    (-mouse_delta.y * SENSITIVITY)
                        .min(angle - max_angle)
                        .max(angle + max_angle - std::f32::consts::PI),
                )
                * character.look_at;
        }

        let pos = transform.translation;
        transform.look_at(pos + character.look_at, character.up);
    }

    // movement
    let mut input = Vec3::new(
        (keys.pressed(KeyCode::D) as i32 - keys.pressed(KeyCode::A) as i32) as f32,
        (keys.pressed(KeyCode::Space) as i32 - keys.pressed(KeyCode::LShift) as i32) as f32,
        (keys.pressed(KeyCode::S) as i32 - keys.pressed(KeyCode::W) as i32) as f32,
    );
    if input != Vec3::ZERO {
        input = input.normalize();
    }
    input *= character.speed;

    let target_velocity = input.z * transform.local_z()
        + input.x * transform.local_x()
        + input.y * transform.local_y();

    let acceleration = if character.grounded { 0.03 } else { 0.01 };

    character.velocity = lerp(
        character.velocity,
        target_velocity,
        acceleration,
        time.delta_seconds(),
    );

    character.up = slerp(
        character.up.normalize(),
        Vec3::Y,
        0.04,
        time.delta_seconds(),
    );

    transform.translation += character.velocity * time.delta_seconds();
}

fn lerp(i: Vec3, f: Vec3, s: f32, dt: f32) -> Vec3 {
    let s = (1.0 - s).powf(dt * 120.0);
    i * s + f * (1.0 - s)
}

// https://youtu.be/ibkT5ao8kGY
fn slerp(i: Vec3, f: Vec3, s: f32, dt: f32) -> Vec3 {
    let s = (1.0 - s).powf(dt * 120.0);
    let theta = i.dot(f).acos();
    if theta.sin() == 0.0 {
        return i + Vec3::splat(0.00000001);
    }
    ((s * theta).sin() / theta.sin()) * i + (((1.0 - s) * theta).sin() / theta.sin()) * f
}
