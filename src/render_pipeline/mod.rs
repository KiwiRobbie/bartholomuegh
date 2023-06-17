use self::main_pass::{MainPassNode, MainPassPlugin};
use bevy::{
    core_pipeline::upscaling::UpscalingNode,
    // core_pipeline::upscaling::UpscalingNode,
    prelude::*,
    render::{
        extract_resource::{ExtractResource, ExtractResourcePlugin},
        render_graph::{RenderGraphApp, ViewNodeRunner},
        RenderApp,
    },
};
pub use main_pass::{IntegrationMethod, MainPassSettings};

mod main_pass;

pub struct RenderPlugin;

impl Plugin for RenderPlugin {
    fn build(&self, app: &mut App) {
        app.insert_resource(RenderGraphSettings::default())
            .add_plugin(ExtractResourcePlugin::<RenderGraphSettings>::default())
            .add_plugin(MainPassPlugin);

        let render_app = match app.get_sub_app_mut(RenderApp) {
            Ok(render_app) => render_app,
            Err(_) => return,
        };

        render_app
            .add_render_sub_graph("voxel")
            .add_render_graph_node::<ViewNodeRunner<MainPassNode>>("voxel", "main_pass")
            .add_render_graph_node::<ViewNodeRunner<UpscalingNode>>("voxel", "upscaling")
            .add_render_graph_edges("voxel", &["main_pass", "upscaling"]);
    }
}

#[derive(Resource, Clone, ExtractResource)]
pub struct RenderGraphSettings {
    pub clear: bool,
    pub automata: bool,
    pub animation: bool,
    pub voxelization: bool,
    pub rebuild: bool,
    pub physics: bool,
    pub trace: bool,
    pub denoise: bool,
}

impl Default for RenderGraphSettings {
    fn default() -> Self {
        Self {
            clear: true,
            automata: true,
            animation: true,
            voxelization: true,
            rebuild: true,
            physics: true,
            trace: true,
            denoise: false,
        }
    }
}
