use crate::{
    character::{self, CharacterEntity},
    render_pipeline::MainPassSettings,
};
use bevy::{prelude::*, render::camera::CameraRenderGraph};

use bevy::{
    core_pipeline::{bloom::BloomSettings, fxaa::Fxaa, tonemapping::Tonemapping},
    prelude::*,
    reflect::TypeRegistryInternal,
    window::PrimaryWindow,
};
use bevy_inspector_egui::{
    bevy_egui::{EguiContexts, EguiPlugin},
    egui::{self, DragValue, Label, Slider},
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
    )>,
    window: Query<Entity, With<PrimaryWindow>>,
    camera: Query<&Transform>,
) {
    egui::Window::new("Settings")
        .anchor(egui::Align2::RIGHT_TOP, [-5.0, 5.0])
        .show(contexts.ctx_for_window_mut(window.single()), |ui| {
            for (i, (mut trace_settings, bloom_settings, tonemapping, fxaa, projection)) in
                camera_settings_query.iter_mut().enumerate()
            {
                ui.collapsing(format!("Camera Settings {}", i), |ui| {
                    ui.checkbox(&mut trace_settings.show_ray_steps, "Show ray steps");
                    ui.checkbox(&mut trace_settings.indirect_lighting, "Indirect lighting");
                    ui.checkbox(&mut trace_settings.shadows, "Shadows");
                    ui.checkbox(&mut trace_settings.misc_bool, "Misc");
                    ui.add(
                        Slider::new(&mut trace_settings.misc_float, 0.0001..=0.1)
                            .text("Misc")
                            .logarithmic(true),
                    );
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
                });
            }
            for transform in &camera {
                ui.add(Label::new(&transform.translation.x.to_string()));
                ui.add(Label::new(&transform.translation.y.to_string()));
                ui.add(Label::new(&transform.translation.z.to_string()));
            }
            // let test = character.single();
            // if let CharacterEntity(character) = character.single() {}
        });
}
