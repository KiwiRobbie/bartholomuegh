use crate::{
    character::{self, CharacterEntity},
    render_pipeline::MainPassSettings,
};
use bevy::{prelude::*, render::camera::Camera, render::camera::CameraRenderGraph};

use bevy::{
    core_pipeline::{bloom::BloomSettings, fxaa::Fxaa, tonemapping::Tonemapping},
    prelude::*,
    reflect::TypeRegistryInternal,
    window::PrimaryWindow,
};
use bevy_inspector_egui::{
    bevy_egui::{EguiContexts, EguiPlugin},
    egui::{self, DragValue, Grid, Slider},
    reflect_inspector::ui_for_value,
};

pub struct UiPlugin;

impl Plugin for UiPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugin(EguiPlugin).add_system(ui_system);
    }
}

fn ui_system(
    mut contexts: EguiContexts,
    mut camera_settings_query: Query<(
        &mut MainPassSettings,
        Option<&mut BloomSettings>,
        Option<&mut Tonemapping>,
        Option<&mut Fxaa>,
        Option<&mut Projection>,
        Option<&mut Transform>,
    )>,
    window: Query<Entity, With<PrimaryWindow>>,
) {
    egui::Window::new("Settings")
        .anchor(egui::Align2::RIGHT_TOP, [-5.0, 5.0])
        .show(contexts.ctx_for_window_mut(window.single()), |ui| {
            for (
                i,
                (mut trace_settings, bloom_settings, tonemapping, fxaa, projection, mut transform),
            ) in camera_settings_query.iter_mut().enumerate()
            {
                ui.collapsing(format!("Camera Settings {}", i), |ui| {
                    ui.heading("Integration");
                    ui.add(
                        Slider::new(&mut trace_settings.step_size, 0.0001..=0.1)
                            .text("Step Size")
                            .logarithmic(true),
                    );
                    ui.add(
                        Slider::new(&mut trace_settings.step_count, 1..=10000)
                            .text("Step Count")
                            .logarithmic(true),
                    );
                    ui.heading("Misc");

                    ui.checkbox(&mut trace_settings.show_ray_steps, "Show ray steps");
                    ui.checkbox(&mut trace_settings.indirect_lighting, "Indirect lighting");
                    ui.checkbox(&mut trace_settings.shadows, "Shadows");
                    ui.checkbox(&mut trace_settings.misc_bool, "Misc");

                    if let Some(bloom_settings) = bloom_settings {
                        ui.add(
                            Slider::new(&mut bloom_settings.into_inner().intensity, 0.0..=1.0)
                                .text("Bloom"),
                        );
                    }

                    let registry = &TypeRegistryInternal::default();
                    if let Some(tonemapping) = tonemapping {
                        ui_for_value(tonemapping.into_inner(), ui, registry);
                    }
                    // if let Some(bloom_settings) = bloom_settings {
                    //     ui_for_value(bloom_settings.into_inner(), ui, registry);
                    // }
                    if let Some(fxaa) = fxaa {
                        ui_for_value(fxaa.into_inner(), ui, registry);
                    }

                    if let Some(projection) = projection {
                        match projection.into_inner() {
                            Projection::Orthographic(orthographic_projection) => {
                                ui.add(
                                    Slider::new(&mut orthographic_projection.scale, 0.0..=1000.0)
                                        .text("Orthographic scale"),
                                );
                            }
                            Projection::Perspective(perspective_projection) => {
                                ui.add(
                                    Slider::new(&mut perspective_projection.fov, 0.0..=3.1415)
                                        .logarithmic(true)
                                        .text("Perspective fov"),
                                );
                            }
                        }
                    }
                    ui.heading("Transform");
                    if let Some(mut transform) = transform {
                        Grid::new("Position").show(ui, |ui| {
                            ui.add(DragValue::new(&mut transform.translation.x).speed(0.1));
                            ui.add(DragValue::new(&mut transform.translation.y).speed(0.1));
                            ui.add(DragValue::new(&mut transform.translation.z).speed(0.1));
                        });
                        // ui.add(Label::new(&transform.translation.x.to_string()));
                        // ui.add(Label::new(&transform.translation.y.to_string()));
                        // ui.add(Label::new(&transform.translation.z.to_string()));
                    }
                });
            }

            // let test = character.single();
            // if let CharacterEntity(character) = character.single() {}
        });
}
