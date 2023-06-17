use crate::main;
use crate::render_pipeline::IntegrationMethod;
use crate::render_pipeline::KerrPassSettings;
use crate::render_pipeline::MainPassSettings;
use crate::render_pipeline::MainPasses;
use crate::render_pipeline::SchwarzschildPassSettings;
use bevy::diagnostic::Diagnostics;
use bevy::diagnostic::FrameTimeDiagnosticsPlugin;
use bevy::render::render_graph;

use bevy::{
    core_pipeline::{bloom::BloomSettings, fxaa::Fxaa, tonemapping::Tonemapping},
    prelude::*,
    reflect::TypeRegistryInternal,
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
        &mut KerrPassSettings,
        &mut SchwarzschildPassSettings,
        &mut Transform,
        Option<&mut BloomSettings>,
        Option<&mut Tonemapping>,
        Option<&mut Fxaa>,
        Option<&mut Projection>,
    )>,
    window: Query<Entity, With<PrimaryWindow>>,
    diagnostics: Res<Diagnostics>,
) {
    // voxel_graph.add_node("kerr_pass", kerr_node);
    // voxel_graph.add_slot_edge(input_node_id, "view_entity", "kerr_pass", "view");

    // RenderGraph::remove_slot_edge(&mut self, output_node, input_node);

    // //

    egui::Window::new("Settings")
        .anchor(egui::Align2::RIGHT_TOP, [-5.0, 5.0])
        .show(contexts.ctx_for_window_mut(window.single()), |ui| {
            if let Some(fps) = diagnostics.get_measurement(FrameTimeDiagnosticsPlugin::FPS) {
                ui.label(format!("{:.3}FPS", fps.value));
            }

            for (
                i,
                (
                    mut main_pass_settings,
                    mut kerr_pass_settings,
                    mut schwarzschild_pass_settings,
                    mut transform,
                    bloom_settings,
                    tonemapping,
                    fxaa,
                    projection,
                ),
            ) in camera_settings_query.iter_mut().enumerate()
            {
                if ui.button("Switch Metric").clicked() {
                    main_pass_settings.updated = true;
                    main_pass_settings.pass = match main_pass_settings.pass {
                        MainPasses::Kerr => MainPasses::Schwarzschild,
                        MainPasses::Schwarzschild => MainPasses::Kerr,
                    }
                } else {
                    main_pass_settings.updated = false;
                }

                match main_pass_settings.pass {
                    MainPasses::Kerr => {
                        CollapsingHeader::new(format!("Camera Settings {}", i))
                            .default_open(true)
                            .show(ui, |ui| {
                                CollapsingHeader::new("Integration")
                                    .default_open(true)
                                    .show(ui, |ui| {
                                        ui.add(
                                            Slider::new(
                                                &mut kerr_pass_settings.step_count,
                                                1..=10000,
                                            )
                                            .text("Step Count")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(
                                                &mut kerr_pass_settings.rel_error,
                                                0.00000001..=0.1,
                                            )
                                            .text("Relative Error Tolerance")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(
                                                &mut kerr_pass_settings.abs_error,
                                                0.00000001..=0.1,
                                            )
                                            .text("Absolute Error Tolerance")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(
                                                &mut kerr_pass_settings.initial_step,
                                                0.0001..=0.1,
                                            )
                                            .text("Initial Step")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(
                                                &mut kerr_pass_settings.max_step,
                                                0.0001..=1.0,
                                            )
                                            .text("Max Step")
                                            .logarithmic(true),
                                        );
                                        ComboBox::from_label("Integration Method")
                                            .selected_text(match kerr_pass_settings.method {
                                                IntegrationMethod::Rk4 => "RK4",
                                                IntegrationMethod::Euler => "Euler",
                                            })
                                            .show_ui(ui, |ui| {
                                                ui.selectable_value(
                                                    &mut kerr_pass_settings.method,
                                                    IntegrationMethod::Euler,
                                                    "Euler",
                                                );
                                                ui.selectable_value(
                                                    &mut kerr_pass_settings.method,
                                                    IntegrationMethod::Rk4,
                                                    "RK4",
                                                );
                                            });
                                    });
                                CollapsingHeader::new("Blackhole").default_open(true).show(
                                    ui,
                                    |ui| {
                                        ui.checkbox(
                                            &mut kerr_pass_settings.surface_bool,
                                            "Debug Surface",
                                        );
                                        ui.checkbox(
                                            &mut kerr_pass_settings.disk_bool,
                                            "Debug Disk",
                                        );
                                        ui.checkbox(&mut kerr_pass_settings.disk_hide, "Hide Disk");

                                        ui.add(
                                            Slider::new(
                                                &mut kerr_pass_settings.disk_start,
                                                1.0..=100.0,
                                            )
                                            .text("Inner Radius")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(
                                                &mut kerr_pass_settings.disk_end,
                                                1.0..=100.0,
                                            )
                                            .text("Outer Radius")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(&mut kerr_pass_settings.spin, 0.0..=1.0)
                                                .text("Spin")
                                                .logarithmic(true),
                                        );
                                    },
                                );

                                CollapsingHeader::new("Misc")
                                    .default_open(true)
                                    .show(ui, |ui| {
                                        if let Some(bloom_settings) = bloom_settings {
                                            ui.add(
                                                Slider::new(
                                                    &mut bloom_settings.into_inner().intensity,
                                                    0.0..=1.0,
                                                )
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
                                                Projection::Orthographic(
                                                    orthographic_projection,
                                                ) => {
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

                                        ui.checkbox(&mut kerr_pass_settings.misc_bool, "Misc");
                                    });

                                CollapsingHeader::new("Transform").default_open(true).show(
                                    ui,
                                    |ui| {
                                        ui.horizontal(|ui| {
                                            ui.label(
                                                RichText::new("X: ")
                                                    .color(Color32::from_rgb(204, 25, 38)),
                                            );
                                            ui.add(
                                                DragValue::new(&mut transform.translation.x)
                                                    .speed(0.1),
                                            );
                                            ui.label(
                                                RichText::new("Y: ")
                                                    .color(Color32::from_rgb(51, 178, 51)),
                                            );
                                            ui.add(
                                                DragValue::new(&mut transform.translation.y)
                                                    .speed(0.1),
                                            );
                                            ui.label(
                                                RichText::new("Z: ")
                                                    .color(Color32::from_rgb(25, 63, 204)),
                                            );
                                            ui.add(
                                                DragValue::new(&mut transform.translation.z)
                                                    .speed(0.1),
                                            );
                                        });
                                    },
                                );
                            });
                    }

                    MainPasses::Schwarzschild => {
                        CollapsingHeader::new(format!("Camera Settings {}", i))
                            .default_open(true)
                            .show(ui, |ui| {
                                CollapsingHeader::new("Integration")
                                    .default_open(true)
                                    .show(ui, |ui| {
                                        ui.add(
                                            Slider::new(
                                                &mut schwarzschild_pass_settings.step_count,
                                                1..=10000,
                                            )
                                            .text("Step Count")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(
                                                &mut schwarzschild_pass_settings.rel_error,
                                                0.00000001..=0.1,
                                            )
                                            .text("Relative Error Tolerance")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(
                                                &mut schwarzschild_pass_settings.abs_error,
                                                0.00000001..=0.1,
                                            )
                                            .text("Absolute Error Tolerance")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(
                                                &mut schwarzschild_pass_settings.initial_step,
                                                0.0001..=0.1,
                                            )
                                            .text("Initial Step")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(
                                                &mut schwarzschild_pass_settings.max_step,
                                                0.0001..=1.0,
                                            )
                                            .text("Max Step")
                                            .logarithmic(true),
                                        );
                                        ComboBox::from_label("Integration Method")
                                            .selected_text(
                                                match schwarzschild_pass_settings.method {
                                                    IntegrationMethod::Rk4 => "RK4",
                                                    IntegrationMethod::Euler => "Euler",
                                                },
                                            )
                                            .show_ui(ui, |ui| {
                                                ui.selectable_value(
                                                    &mut schwarzschild_pass_settings.method,
                                                    IntegrationMethod::Euler,
                                                    "Euler",
                                                );
                                                ui.selectable_value(
                                                    &mut schwarzschild_pass_settings.method,
                                                    IntegrationMethod::Rk4,
                                                    "RK4",
                                                );
                                            });
                                    });
                                CollapsingHeader::new("Blackhole").default_open(true).show(
                                    ui,
                                    |ui| {
                                        ui.checkbox(
                                            &mut schwarzschild_pass_settings.surface_bool,
                                            "Debug Surface",
                                        );
                                        ui.checkbox(
                                            &mut schwarzschild_pass_settings.disk_bool,
                                            "Debug Disk",
                                        );
                                        ui.checkbox(
                                            &mut schwarzschild_pass_settings.disk_hide,
                                            "Hide Disk",
                                        );

                                        ui.add(
                                            Slider::new(
                                                &mut schwarzschild_pass_settings.disk_start,
                                                1.0..=100.0,
                                            )
                                            .text("Inner Radius")
                                            .logarithmic(true),
                                        );
                                        ui.add(
                                            Slider::new(
                                                &mut schwarzschild_pass_settings.disk_end,
                                                1.0..=100.0,
                                            )
                                            .text("Outer Radius")
                                            .logarithmic(true),
                                        );
                                    },
                                );

                                CollapsingHeader::new("Misc")
                                    .default_open(true)
                                    .show(ui, |ui| {
                                        if let Some(bloom_settings) = bloom_settings {
                                            ui.add(
                                                Slider::new(
                                                    &mut bloom_settings.into_inner().intensity,
                                                    0.0..=1.0,
                                                )
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
                                                Projection::Orthographic(
                                                    orthographic_projection,
                                                ) => {
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

                                        ui.add(
                                            Slider::new(
                                                &mut schwarzschild_pass_settings.reletavistic_scale,
                                                0.0..=1.0,
                                            )
                                            .text("Scale"),
                                        );

                                        ui.checkbox(
                                            &mut schwarzschild_pass_settings.misc_bool,
                                            "Misc",
                                        );
                                    });

                                CollapsingHeader::new("Transform").default_open(true).show(
                                    ui,
                                    |ui| {
                                        ui.horizontal(|ui| {
                                            ui.label(
                                                RichText::new("X: ")
                                                    .color(Color32::from_rgb(204, 25, 38)),
                                            );
                                            ui.add(
                                                DragValue::new(&mut transform.translation.x)
                                                    .speed(0.1),
                                            );
                                            ui.label(
                                                RichText::new("Y: ")
                                                    .color(Color32::from_rgb(51, 178, 51)),
                                            );
                                            ui.add(
                                                DragValue::new(&mut transform.translation.y)
                                                    .speed(0.1),
                                            );
                                            ui.label(
                                                RichText::new("Z: ")
                                                    .color(Color32::from_rgb(25, 63, 204)),
                                            );
                                            ui.add(
                                                DragValue::new(&mut transform.translation.z)
                                                    .speed(0.1),
                                            );
                                        });
                                    },
                                );
                            });
                    }
                }
            }
        });
}
