use bevy::diagnostic::FrameTimeDiagnosticsPlugin;
use bevy::{core_pipeline::bloom::BloomSettings, prelude::*, render::camera::CameraRenderGraph};
use character::CharacterEntity;
use render_pipeline::{KerrPassSettings, MainPassSettings, SchwarzschildPassSettings};

mod character;
mod render_pipeline;
mod ui;

fn main() {
    #[cfg(target_arch = "wasm32")]
    console_error_panic_hook::set_once();

    let mut app = App::new();
    app.add_plugins(
        DefaultPlugins
            .set(AssetPlugin {
                watch_for_changes: true,
                ..default()
            })
            .set(WindowPlugin {
                primary_window: Some(Window {
                    fit_canvas_to_parent: true,
                    ..default()
                }),
                ..default()
            }),
    )
    .add_plugin(FrameTimeDiagnosticsPlugin::default())
    .add_plugin(render_pipeline::RenderPlugin)
    .add_plugin(character::CharacterPlugin)
    .add_plugin(ui::UiPlugin)
    .add_startup_system(setup);

    #[cfg(not(target_arch = "wasm32"))]
    {
        let settings = bevy_mod_debugdump::render_graph::Settings::default();
        let dot = bevy_mod_debugdump::render_graph_dot(&mut app, &settings);
        std::fs::write("render-graph.dot", dot).expect("Failed to write render-graph.dot");
    }

    app.run();
}

fn setup(mut commands: Commands) {
    let character_transform = Transform::from_xyz(10.0, 1.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y);
    commands.spawn((
        Camera3dBundle {
            transform: character_transform,
            camera_render_graph: CameraRenderGraph::new("main_render"),
            camera: Camera {
                hdr: true,
                ..default()
            },
            projection: Projection::Perspective(PerspectiveProjection {
                fov: 1.57,
                near: 0.001,
                far: 100.0,
                ..default()
            }),
            ..default()
        },
        KerrPassSettings::default(),
        SchwarzschildPassSettings::default(),
        BloomSettings::default(),
        MainPassSettings::default(),
        CharacterEntity {
            velocity: Vec3::ZERO,
            speed: 10.0,
            grounded: false,
            look_at: -character_transform.local_z(),
            up: Vec3::Y,
        },
    ));
}
