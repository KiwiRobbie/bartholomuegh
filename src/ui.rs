use crate::render_pipeline::{
    KerrSettings, MainPassSettings, MetricSettings, SchwarzschildSettings,
};
use bevy::{
    core_pipeline::{bloom::BloomSettings, fxaa::Fxaa, tonemapping::Tonemapping},
    diagnostic::{Diagnostics, FrameTimeDiagnosticsPlugin},
    prelude::*,
    window::PrimaryWindow,
};
use bevy_inspector_egui::{
    bevy_egui::{EguiContexts, EguiPlugin},
    egui::{self, CollapsingHeader, Color32, ComboBox, DragValue, RichText, Slider},
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
        &mut Transform,
        Option<&mut Projection>,
    )>,
    mut two_d_camera_query: Query<(&mut BloomSettings, &mut Tonemapping, &mut Fxaa)>,
    window: Query<Entity, With<PrimaryWindow>>,
    diagnostics: Res<Diagnostics>,
    type_registry: ResMut<AppTypeRegistry>,
) {
    egui::Window::new("Settings")
        .anchor(egui::Align2::RIGHT_TOP, [-5.0, 5.0])
        .show(contexts.ctx_for_window_mut(window.single()), |ui| {
            if let Some(fps) = diagnostics.get_measurement(FrameTimeDiagnosticsPlugin::FPS) {
                ui.label(format!("{:.3}FPS", fps.value));
            }

            for (i, (mut main_pass_settings, mut transform, projection)) in
                camera_settings_query.iter_mut().enumerate()
            {
                CollapsingHeader::new(format!("Camera Settings {}", i))
                    .default_open(true)
                    .show(ui, |ui| {
                        ComboBox::from_label("Metric")
                            .selected_text(format!("{:?}", main_pass_settings.metric))
                            .show_ui(ui, |ui| {
                                ui.selectable_value(
                                    &mut main_pass_settings.metric,
                                    MetricSettings::Schwarzschild(SchwarzschildSettings::default()),
                                    "Schwarzschild",
                                );
                                ui.selectable_value(
                                    &mut main_pass_settings.metric,
                                    MetricSettings::Kerr(KerrSettings::default()),
                                    "Kerr",
                                );
                            });

                        CollapsingHeader::new(format!("Metric Settings"))
                            .default_open(true)
                            .show(ui, |ui| match main_pass_settings.metric {
                                MetricSettings::Kerr(ref mut kerr_settings) => {
                                    ui_for_value(kerr_settings, ui, &type_registry.read());
                                }
                                MetricSettings::Schwarzschild(ref mut schwarzschild_settings) => {
                                    ui_for_value(schwarzschild_settings, ui, &type_registry.read());
                                }
                            });

                        CollapsingHeader::new("General")
                            .default_open(true)
                            .show(ui, |ui| {
                                if let Some(projection) = projection {
                                    match projection.into_inner() {
                                        Projection::Orthographic(orthographic_projection) => {
                                            ui.add(
                                                Slider::new(
                                                    &mut orthographic_projection.scale,
                                                    0.0..=1000.0,
                                                )
                                                .text("Orthographic scale"),
                                            );
                                        }
                                        Projection::Perspective(perspective_projection) => {
                                            ui.add(
                                                Slider::new(
                                                    &mut perspective_projection.fov,
                                                    0.0..=3.1415,
                                                )
                                                .logarithmic(true)
                                                .text("Perspective fov"),
                                            );
                                        }
                                    }
                                }

                                ui.horizontal(|ui| {
                                    ui.label(
                                        RichText::new("X: ").color(Color32::from_rgb(204, 25, 38)),
                                    );
                                    ui.add(DragValue::new(&mut transform.translation.x).speed(0.1));
                                    ui.label(
                                        RichText::new("Y: ").color(Color32::from_rgb(51, 178, 51)),
                                    );
                                    ui.add(DragValue::new(&mut transform.translation.y).speed(0.1));
                                    ui.label(
                                        RichText::new("Z: ").color(Color32::from_rgb(25, 63, 204)),
                                    );
                                    ui.add(DragValue::new(&mut transform.translation.z).speed(0.1));
                                });
                            });
                    });
            }

            CollapsingHeader::new("Secondary Camera")
                .default_open(true)
                .show(ui, |ui| {
                    let (bloom_settings, tonemapping, fxaa) = two_d_camera_query.single_mut();

                    ui.add(
                        Slider::new(&mut bloom_settings.into_inner().intensity, 0.0..=1.0)
                            .text("Bloom"),
                    );

                    ui_for_value(tonemapping.into_inner(), ui, &type_registry.read());
                    ui_for_value(fxaa.into_inner(), ui, &type_registry.read());
                });
        });
}
