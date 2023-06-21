use bevy::{
    core_pipeline::{bloom::BloomSettings, fxaa::Fxaa, tonemapping::Tonemapping},
    diagnostic::FrameTimeDiagnosticsPlugin,
    prelude::*,
    render::{
        camera::{CameraRenderGraph, RenderTarget},
        render_resource::{Extent3d, TextureDimension, TextureFormat, TextureUsages},
    },
    window::{PrimaryWindow, WindowResized, WindowScaleFactorChanged},
};
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
    .add_startup_system(setup)
    .add_system(update_render_texture);

    #[cfg(not(target_arch = "wasm32"))]
    {
        let settings = bevy_mod_debugdump::render_graph::Settings::default();
        let dot = bevy_mod_debugdump::render_graph_dot(&mut app, &settings);
        std::fs::write("render-graph.dot", dot).expect("Failed to write render-graph.dot");
    }

    app.run();
}

#[derive(Resource)]
struct CameraData {
    image: Handle<Image>,
    sprite: Entity,
}

fn setup(mut commands: Commands, mut images: ResMut<Assets<Image>>) {
    let mut image = Image::new_fill(
        Extent3d {
            width: 1,
            height: 1,
            depth_or_array_layers: 1,
        },
        TextureDimension::D2,
        &[0, 255, 0, 255],
        TextureFormat::Rgba8UnormSrgb,
    );
    image.texture_descriptor.usage =
        TextureUsages::TEXTURE_BINDING | TextureUsages::RENDER_ATTACHMENT;
    let image = images.add(image);

    let character_transform = Transform::from_xyz(10.0, 1.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y);
    commands.spawn((
        Camera3dBundle {
            transform: character_transform,
            camera_render_graph: CameraRenderGraph::new("main_render"),
            camera: Camera {
                hdr: true,
                target: RenderTarget::Image(image.clone()),
                order: -10,
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
        MainPassSettings::default(),
        SchwarzschildPassSettings::default(),
        KerrPassSettings::default(),
        CharacterEntity {
            velocity: Vec3::ZERO,
            speed: 10.0,
            grounded: false,
            look_at: -character_transform.local_z(),
            up: Vec3::Y,
        },
    ));

    let sprite = commands
        .spawn(SpriteBundle {
            texture: image.clone(),
            ..default()
        })
        .id();
    commands.spawn((
        Camera2dBundle {
            camera: Camera {
                hdr: true,
                ..default()
            },
            tonemapping: Tonemapping::TonyMcMapface,
            ..default()
        },
        BloomSettings::default(),
        Fxaa::default(),
    ));

    commands.insert_resource(CameraData { image, sprite });
}

fn update_render_texture(
    mut resize_reader: EventReader<WindowResized>,
    mut scale_factor_reader: EventReader<WindowScaleFactorChanged>,
    mut images: ResMut<Assets<Image>>,
    mut sprites: Query<&mut Sprite>,
    render_image: Res<CameraData>,
    windows: Query<&Window, With<PrimaryWindow>>,
) {
    let window = windows.single();

    let mut update = |width: f32, height: f32| {
        let new_size = Extent3d {
            width: width as u32,
            height: height as u32,
            depth_or_array_layers: 1,
        };

        let image = images.get_mut(&render_image.image).unwrap();
        image.resize(new_size);

        let mut sprite = sprites.get_mut(render_image.sprite).unwrap();
        sprite.custom_size = Some(Vec2::new(width, height));
    };

    for _ in resize_reader.iter() {
        update(window.width(), window.height());
    }

    for _ in scale_factor_reader.iter() {
        update(window.width(), window.height());
    }
}
