use bevy::{prelude::*, render::camera::CameraRenderGraph};
use character::CharacterEntity;
use render_pipeline::MainPassSettings;
use test;

mod character;
mod render_pipeline;
mod ui;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(AssetPlugin {
            watch_for_changes: true,
            ..default()
        }))
        .add_plugin(render_pipeline::RenderPlugin)
        .add_plugin(character::CharacterPlugin)
        .add_plugin(ui::UiPlugin)
        .add_startup_system(setup)
        .run();
}

fn setup(mut commands: Commands) {
    let character_transform = Transform::from_xyz(2.0, 2.0, -1.0).looking_at(Vec3::ZERO, Vec3::Y);
    commands.spawn((
        Camera3dBundle {
            transform: character_transform,
            camera_render_graph: CameraRenderGraph::new("voxel"),
            camera: Camera {
                hdr: true,
                ..default()
            },
            projection: Projection::Perspective(PerspectiveProjection {
                fov: 1.57,
                ..default()
            }),
            ..default()
        },
        MainPassSettings::default(),
        CharacterEntity {
            velocity: Vec3::ZERO,
            speed: 10.0,
            grounded: false,
            in_spectator: true,
            look_at: -character_transform.local_z(),
            up: Vec3::Y,
        },
    ));
}
