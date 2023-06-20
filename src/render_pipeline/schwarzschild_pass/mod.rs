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
pub use node::SchwarzschildPassNode;
mod node;

pub struct SchwarzschildPassPlugin;

impl Plugin for SchwarzschildPassPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugin(ExtractComponentPlugin::<SchwarzschildPassSettings>::default());

        // setup custom render pipeline
        app.sub_app_mut(RenderApp)
            .init_resource::<SchwarzschildPassPipelineData>()
            .add_system(prepare_uniforms.in_set(RenderSet::Prepare));
    }
}

#[derive(Resource)]
struct SchwarzschildPassPipelineData {
    pipeline_id: CachedRenderPipelineId,
    bind_group_layout: BindGroupLayout,
}

#[derive(Component, Clone, ExtractComponent)]
pub struct SchwarzschildPassSettings {
    pub surface_bool: bool,
    pub disk_bool: bool,
    pub disk_hide: bool,
    pub step_count: i32,
    pub rel_error: f32,
    pub abs_error: f32,
    pub initial_step: f32,
    pub max_step: f32,
    pub method: IntegrationMethod,
    pub disk_start: f32,
    pub disk_end: f32,
    pub velocity: Vec3,
    pub reletavistic_scale: f32,
    pub misc_float: f32,
    pub misc_bool: bool,
}

impl Default for SchwarzschildPassSettings {
    fn default() -> Self {
        Self {
            surface_bool: false,
            disk_bool: false,
            disk_hide: false,
            step_count: 128,
            rel_error: 1.0E-5,
            abs_error: 1.0E-5,
            initial_step: 0.00050,
            max_step: 0.25,
            method: IntegrationMethod::Rk4,
            disk_start: 1.0,
            disk_end: 100.0,
            velocity: Vec3::ZERO,
            reletavistic_scale: 0.15,
            misc_float: 0.0,
            misc_bool: false,
        }
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

#[derive(Component, Deref, DerefMut)]
pub struct ViewSchwarzschildPassUniformBuffer(UniformBuffer<TraceUniforms>);

fn prepare_uniforms(
    mut commands: Commands,
    query: Query<(Entity, &SchwarzschildPassSettings, &ExtractedView)>,
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

        let camera = projection * inverse_view;
        let camera_inverse = view * inverse_projection;

        let uniforms = TraceUniforms {
            camera,
            camera_inverse,
            velocity: settings.velocity,
            time: elapsed as f32,
            surface_bool: settings.surface_bool as u32,
            disk_mode: (!settings.disk_hide) as u32
                + (!settings.disk_hide && settings.disk_bool) as u32,
            step_count: settings.step_count,
            rel_error: settings.rel_error,
            abs_error: settings.abs_error,
            initial_step: settings.initial_step,
            max_step: settings.max_step,
            method: match settings.method {
                IntegrationMethod::Rk4 => 0,
                IntegrationMethod::Euler => 1,
            },
            disk_start: settings.disk_start,
            disk_end: settings.disk_end,
            misc_float: settings.misc_float,
            misc_bool: settings.misc_bool as u32,
        };

        let mut uniform_buffer = UniformBuffer::from(uniforms);
        uniform_buffer.write_buffer(&render_device, &render_queue);

        commands
            .entity(entity)
            .insert(ViewSchwarzschildPassUniformBuffer(uniform_buffer));
    }
}

impl FromWorld for SchwarzschildPassPipelineData {
    fn from_world(render_world: &mut World) -> Self {
        let asset_server = render_world.get_resource::<AssetServer>().unwrap();

        let bind_group_layout = render_world
            .resource::<RenderDevice>()
            .create_bind_group_layout(&BindGroupLayoutDescriptor {
                label: Some("schwarzschild bind group layout"),
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

        let trace_shader = asset_server.load("schwarzschild.wgsl");

        let trace_pipeline_descriptor = RenderPipelineDescriptor {
            label: Some("schwarzschild pipeline".into()),
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

        SchwarzschildPassPipelineData {
            pipeline_id,
            bind_group_layout,
        }
    }
}
