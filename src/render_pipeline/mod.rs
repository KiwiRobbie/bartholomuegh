use self::kerr_pass::{KerrPassNode, KerrPassPlugin};
use self::schwarzschild_pass::{SchwarzschildPassNode, SchwarzschildPassPlugin};
use bevy::{
    core_pipeline::{
        bloom::BloomNode, fxaa::FxaaNode, tonemapping::TonemappingNode, upscaling::UpscalingNode,
    },
    prelude::*,
    render::{
        extract_component::{ExtractComponent, ExtractComponentPlugin},
        extract_resource::{ExtractResource, ExtractResourcePlugin},
        render_graph::{RenderGraph, SlotInfo, SlotType},
        RenderApp,
    },
    ui::UiPassNode,
};
pub use kerr_pass::KerrPassSettings;
pub use schwarzschild_pass::SchwarzschildPassSettings;

mod kerr_pass;
mod schwarzschild_pass;

#[derive(Component, Clone, ExtractComponent, Default, Debug)]
pub struct MainPassSettings {
    pub updated: bool,
    pub pass: MainPasses,
}

#[derive(Component, Clone, Default, Debug, PartialEq, Eq)]
pub enum MainPasses {
    Kerr,
    #[default]
    Schwarzschild,
}

#[derive(Clone, Copy, PartialEq, Default)]
pub enum IntegrationMethod {
    Rk4,
    #[default]
    Euler,
}

pub struct RenderPlugin;

impl Plugin for RenderPlugin {
    fn build(&self, app: &mut App) {
        app.insert_resource(RenderGraphSettings::default())
            .add_plugin(ExtractResourcePlugin::<RenderGraphSettings>::default())
            .add_plugin(KerrPassPlugin)
            .add_plugin(SchwarzschildPassPlugin);

        app.add_plugin(ExtractComponentPlugin::<MainPassSettings>::default());
        let render_app = match app.get_sub_app_mut(RenderApp) {
            Ok(render_app) => render_app,
            Err(_) => return,
        };

        // build main render graph
        let mut render_graph = RenderGraph::default();
        let input_node_id =
            render_graph.set_input(vec![SlotInfo::new("view_entity", SlotType::Entity)]);

        let kerr_node = KerrPassNode::new(&mut render_app.world);
        let schwarzschild_node = SchwarzschildPassNode::new(&mut render_app.world);

        let bloom = BloomNode::new(&mut render_app.world);
        let tonemapping = TonemappingNode::new(&mut render_app.world);
        let fxaa = FxaaNode::new(&mut render_app.world);
        let ui = UiPassNode::new(&mut render_app.world);
        let upscaling = UpscalingNode::new(&mut render_app.world);

        render_graph.add_node("schwarzschild_pass", kerr_node);
        render_graph.add_node("kerr_pass", schwarzschild_node);
        render_graph.add_node("bloom", bloom);
        render_graph.add_node("tonemapping", tonemapping);
        render_graph.add_node("fxaa", fxaa);
        render_graph.add_node("ui", ui);
        render_graph.add_node("upscaling", upscaling);
        render_graph.add_slot_edge(input_node_id, "view_entity", "schwarzschild_pass", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "kerr_pass", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "bloom", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "tonemapping", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "fxaa", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "ui", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "upscaling", "view");
        render_graph.add_node_edge("schwarzschild_pass", "bloom");
        render_graph.add_node_edge("kerr_pass", "bloom");
        render_graph.add_node_edge("bloom", "tonemapping");
        render_graph.add_node_edge("tonemapping", "fxaa");
        render_graph.add_node_edge("fxaa", "ui");
        render_graph.add_node_edge("ui", "upscaling");

        let mut graph = render_app.world.resource_mut::<RenderGraph>();
        graph.add_sub_graph("main_render", render_graph);
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
