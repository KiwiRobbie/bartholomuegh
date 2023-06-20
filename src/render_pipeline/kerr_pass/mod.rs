pub use crate::render_pipeline::IntegrationMethod;
use bevy::{
    core_pipeline::fullscreen_vertex_shader::fullscreen_shader_vertex_state,
    prelude::*,
    render::{
        extract_component::{ExtractComponent, ExtractComponentPlugin},
        render_resource::*,
        renderer::{RenderDevice, RenderQueue},
        view::{ExtractedView, ViewTarget},
        RenderApp, RenderSet,
    },
};
pub use node::KerrPassNode;

mod node;

pub struct KerrPassPlugin;

impl Plugin for KerrPassPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugin(ExtractComponentPlugin::<KerrPassSettings>::default());

        // setup custom render pipeline
        app.sub_app_mut(RenderApp)
            .init_resource::<KerrPassPipelineData>()
            .add_system(prepare_uniforms.in_set(RenderSet::Prepare));
    }
}

#[derive(Resource)]
struct KerrPassPipelineData {
    pipeline_id: CachedRenderPipelineId,
    bind_group_layout: BindGroupLayout,
}

#[derive(Component, Clone, ExtractComponent)]
pub struct KerrPassSettings {
    pub surface_bool: bool,
    pub disk_bool: bool,
    pub disk_hide: bool,
    pub misc_bool: bool,
    pub step_count: i32,
    pub rel_error: f32,
    pub abs_error: f32,
    pub initial_step: f32,
    pub max_step: f32,
    pub method: IntegrationMethod,
    pub disk_start: f32,
    pub disk_end: f32,
    pub spin: f32,
}

impl Default for KerrPassSettings {
    fn default() -> Self {
        Self {
            surface_bool: false,
            disk_bool: false,
            disk_hide: false,
            misc_bool: false,
            step_count: 100,
            initial_step: 0.00050,
            rel_error: 0.0000010,
            abs_error: 0.0000010,
            max_step: 1.0,
            method: IntegrationMethod::Rk4,
            disk_start: 1.0,
            disk_end: 12.0,
            spin: 1.0,
        }
    }
}

#[derive(Clone, ShaderType)]
pub struct TraceUniforms {
    pub camera: Mat4,
    pub camera_inverse: Mat4,

    pub time: f32,

    pub surface_bool: u32,
    pub disk_mode: u32,
    pub misc_bool: u32,

    pub step_count: i32,
    pub initial_step: f32,

    pub abs_error: f32,
    pub rel_error: f32,

    pub max_step: f32,

    pub method: u32,

    pub disk_start: f32,
    pub disk_end: f32,

    pub r: f32,
    pub theta: f32,
    pub phi: f32,
    pub a: f32,

    pub rho: f32,
    pub delta: f32,
    pub sigma: f32,
    pub alpha: f32,
    pub omega: f32,
    pub omega_bar: f32,

    pub camera_velocity_r: f32,
    pub camera_velocity_theta: f32,
    pub camera_velocity_phi: f32,
    pub camera_beta: f32,
}

#[derive(Component, Deref, DerefMut)]
struct ViewKerrPassUniformBuffer(UniformBuffer<TraceUniforms>);

fn metric_values(
    r: f32,
    theta: f32,
    phi: f32,
    a: f32,
) -> (f32, f32, f32, f32, f32, f32, f32, f32, f32, f32) {
    let rho = (r.powi(2) + a.powi(2) * theta.cos().powi(2)).sqrt();
    let delta = r.powi(2) - 2.0 * r + a.powi(2);
    let sigma = ((r.powi(2) + a.powi(2)).powi(2) - a.powi(2) * delta * theta.sin().powi(2)).sqrt();
    let alpha = rho * delta.sqrt() / sigma;
    let omega = 2.0 * a * r / sigma.powi(2);
    let omega_bar = sigma * theta.sin() / rho;
    return (r, theta, phi, a, rho, delta, sigma, alpha, omega, omega_bar);
}

fn prepare_uniforms(
    mut commands: Commands,
    query: Query<(Entity, &KerrPassSettings, &ExtractedView)>,
    time: Res<Time>,
    render_device: Res<RenderDevice>,
    render_queue: Res<RenderQueue>,
) {
    let elapsed = time.elapsed_seconds_f64();

    for (entity, settings, view) in query.iter() {
        let projection = view.projection;
        let inverse_projection = projection.inverse();
        let view = view.transform.compute_matrix();
        let inverse_view = view.inverse();

        let pos = Transform::from_matrix(view).translation;

        let r = pos.length(); //Vec3::new(transform.translation.x, 0.0, transform.translation.z).length();
        let theta = f32::acos(pos.y / pos.length());
        let phi: f32 = pos.x.atan2(pos.z); // f32::atan2(transform.translation.x, transform.translation.z);

        let view = view;

        // let camera = projection * inverse_view;
        // let camera_inverse = view * inverse_projection;

        let camera = projection * inverse_view;
        let camera_inverse = view * inverse_projection;

        let a = settings.spin;

        let (r, theta, phi, a, rho, delta, sigma, alpha, omega, omega_bar) =
            metric_values(r, theta, phi, a);

        let camera_velocity = Vec3::new(0.0, 0.0, 1.0);

        let big_omega: f32 = 1.0 / (a + r.powf(3.0 / 2.0));
        let beta = omega_bar / alpha * (big_omega - omega);

        let uniforms = TraceUniforms {
            camera,
            camera_inverse,

            time: elapsed as f32,

            surface_bool: settings.surface_bool as u32,
            disk_mode: (!settings.disk_hide) as u32
                + (!settings.disk_hide && settings.disk_bool) as u32,
            misc_bool: settings.misc_bool as u32,

            step_count: settings.step_count,

            initial_step: settings.initial_step,

            rel_error: settings.rel_error,
            abs_error: settings.abs_error,

            max_step: settings.max_step,

            method: match settings.method {
                IntegrationMethod::Rk4 => 0,
                IntegrationMethod::Euler => 1,
            },

            disk_start: settings.disk_start,
            disk_end: settings.disk_end,

            r,
            theta,
            phi,
            a,

            rho,
            delta,
            sigma,
            alpha,
            omega,
            omega_bar,

            camera_velocity_r: camera_velocity.x,
            camera_velocity_theta: camera_velocity.y,
            camera_velocity_phi: camera_velocity.z,
            camera_beta: beta,
        };

        let mut uniform_buffer = UniformBuffer::from(uniforms);
        uniform_buffer.write_buffer(&render_device, &render_queue);

        commands
            .entity(entity)
            .insert(ViewKerrPassUniformBuffer(uniform_buffer));
    }
}

impl FromWorld for KerrPassPipelineData {
    fn from_world(render_world: &mut World) -> Self {
        let asset_server = render_world.get_resource::<AssetServer>().unwrap();

        let bind_group_layout = render_world
            .resource::<RenderDevice>()
            .create_bind_group_layout(&BindGroupLayoutDescriptor {
                label: Some("kerr bind group layout"),
                entries: &[BindGroupLayoutEntry {
                    binding: 0,
                    visibility: ShaderStages::FRAGMENT,
                    ty: BindingType::Buffer {
                        ty: BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: BufferSize::new(TraceUniforms::SHADER_SIZE.into()),
                    },
                    count: None,
                }],
            });

        let trace_shader = asset_server.load("kerr.wgsl");

        let trace_pipeline_descriptor = RenderPipelineDescriptor {
            label: Some("kerr pipeline".into()),
            layout: vec![bind_group_layout.clone()],
            vertex: fullscreen_shader_vertex_state(),
            fragment: Some(FragmentState {
                shader: trace_shader,
                shader_defs: Vec::new(),
                entry_point: "fragment".into(),
                targets: vec![Some(ColorTargetState {
                    format: ViewTarget::TEXTURE_FORMAT_HDR,
                    blend: None,
                    write_mask: ColorWrites::ALL,
                })],
            }),
            primitive: PrimitiveState::default(),
            depth_stencil: None,
            multisample: MultisampleState::default(),
            push_constant_ranges: Vec::new(),
        };

        let cache = render_world.resource::<PipelineCache>();
        let pipeline_id = cache.queue_render_pipeline(trace_pipeline_descriptor);

        KerrPassPipelineData {
            pipeline_id,
            bind_group_layout,
        }
    }
}
