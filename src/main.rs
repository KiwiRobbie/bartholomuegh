use bevy::{
    asset::ChangeWatcher,
    core_pipeline::{bloom::BloomSettings, fxaa::Fxaa, tonemapping::Tonemapping},
    prelude::*,
    render::{
        camera::{CameraRenderGraph, RenderTarget},
        render_resource::{Extent3d, TextureDimension, TextureFormat, TextureUsages},
    },
    window::{PrimaryWindow, WindowResized, WindowScaleFactorChanged},
};
use character::CharacterEntity;
use render_pipeline::MainPassSettings;
use std::time::Duration;

mod character;
mod render_pipeline;
// mod ui;

fn main() {
    App::new()
        .add_plugins(
            DefaultPlugins
                .set(AssetPlugin {
                    watch_for_changes: ChangeWatcher::with_delay(Duration::from_millis(200)),
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
        .add_plugin(render_pipeline::RenderPlugin)
        .add_plugin(character::CharacterPlugin)
        // .add_plugin(ui::UiPlugin)
        .add_systems(Startup, setup)
        .add_systems(PreUpdate, update_render_texture)
        .run();
}

#[derive(Resource)]
struct CameraData {
    _camera: Entity,
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
        TextureFormat::Rgba16Float,
    );
    image.texture_descriptor.usage =
        TextureUsages::TEXTURE_BINDING | TextureUsages::RENDER_ATTACHMENT;
    let image = images.add(image);

    let character_transform = Transform::from_xyz(2.0, 2.0, -1.0).looking_at(Vec3::ZERO, Vec3::Y);
    let camera = commands
        .spawn((
            Camera3dBundle {
                transform: character_transform,
                camera_render_graph: CameraRenderGraph::new("voxel"),
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
            CharacterEntity {
                velocity: Vec3::ZERO,
                speed: 10.0,
                grounded: false,
                in_spectator: true,
                look_at: -character_transform.local_z(),
                up: Vec3::Y,
            },
        ))
        .id();

    let sprite = commands
        .spawn(SpriteBundle {
            sprite: Sprite {
                custom_size: Some(Vec2::new(512.0, 512.0)),
                ..default()
            },
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

    commands.insert_resource(CameraData {
        _camera: camera,
        image,
        sprite,
    });
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
