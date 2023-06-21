use super::{kerr_pass::KerrSettings, schwarzschild_pass::SchwarzschildSettings};
use bevy::{
    prelude::*,
    render::{extract_component::ExtractComponent, render_resource::*},
};

#[derive(Component, Clone, ExtractComponent, Default)]
pub struct MainPassSettings {
    pub metric: MetricSettings,
}

#[derive(Component, Clone, PartialEq)]
pub enum MetricSettings {
    Kerr(KerrSettings),
    Schwarzschild(SchwarzschildSettings),
}

#[derive(Clone, Copy, PartialEq, Default, Reflect)]
pub enum IntegrationMethod {
    Rk4,
    #[default]
    Euler,
}

impl Default for MetricSettings {
    fn default() -> Self {
        Self::Schwarzschild(SchwarzschildSettings::default())
    }
}

#[derive(Clone, ShaderType)]
pub struct TraceUniforms {
    pub camera: Mat4,
    pub camera_inverse: Mat4,
    pub velocity: Vec3,
    pub time: f32,
    pub surface_bool: u32,
    pub disk_mode: u32,
    pub step_count: i32,
    pub rel_error: f32,
    pub abs_error: f32,
    pub initial_step: f32,
    pub max_step: f32,
    pub method: u32,
    pub disk_start: f32,
    pub disk_end: f32,
    pub misc_float: f32,
    pub misc_bool: u32,
}

// fn prepare_uniforms(
//     mut commands: Commands,
//     query: Query<(Entity, &MainPassSettings, &ExtractedView)>,
//     time: Res<Time>,
//     render_device: Res<RenderDevice>,
//     render_queue: Res<RenderQueue>,
// ) {
//     let elapsed = time.elapsed_seconds_f64();

//     for (entity, settings, view) in query.iter() {
//         let projection = view.projection;
//         let inverse_projection = projection.inverse();
//         let view = view.transform.compute_matrix();
//         let inverse_view = view.inverse();

//         let camera = projection * inverse_view;
//         let camera_inverse = view * inverse_projection;
//     }
// }

// fn queue_uniforms(
//     mut commands: Commands,
//     query: Query<(Entity, &MainPassSettings, &ExtractedView)>,
//     time: Res<Time>,
//     render_device: Res<RenderDevice>,
//     render_queue: Res<RenderQueue>,
// ) {
//     for (entity, settings, view) in query.iter() {
//         let mut uniform_buffer = UniformBuffer::from(uniforms);
//         uniform_buffer.write_buffer(&render_device, &render_queue);

//         commands
//             .entity(entity)
//             .insert(ViewSchwarzschildPassUniformBuffer(uniform_buffer));
//     }
// }

impl std::fmt::Debug for MetricSettings {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MetricSettings::Kerr(_) => write!(f, "Kerr"),
            MetricSettings::Schwarzschild(_) => write!(f, "Schwarzschild"),
        }
    }
}