use self::kerr_pass::{KerrPassNode, KerrPassPlugin};
use self::schwarzschild_pass::{SchwarzschildPassNode, SchwarzschildPassPlugin};

use bevy::render::extract_component::{ExtractComponent, ExtractComponentPlugin};
use bevy::render::RenderSet;
use bevy::{
    core_pipeline::{
        bloom::BloomNode, fxaa::FxaaNode, tonemapping::TonemappingNode, upscaling::UpscalingNode,
    },
    prelude::*,
    render::{
        extract_resource::{ExtractResource, ExtractResourcePlugin},
        render_graph::{RenderGraph, SlotInfo, SlotType},
        RenderApp,
    },
    ui::UiPassNode,
};
pub use common::IntegrationMethod;
pub use kerr_pass::KerrPassSettings;
pub use schwarzschild_pass::SchwarzschildPassSettings;

mod common;
mod kerr_pass;
mod schwarzschild_pass;

fn change_state(
    query: Query<(
        &MainPassSettings,
        Changed<MainPassSettings>, // Why does this not work????
    )>,
    mut graph: ResMut<RenderGraph>,
) {
    for settings in query.iter() {
        if let (settings, true) = settings {
            if settings.updated {
                let render_graph = graph
                    .get_sub_graph_mut("main_render")
                    .expect("Failed switching passes: Unable to obtain 'main_render' subgraph!");

                let input_node_id = render_graph.get_input_node().unwrap().id;
                let (old_pass, new_pass) = match settings.pass {
                    MainPasses::Kerr => ("schwarzschild_pass", "kerr_pass"),
                    MainPasses::Schwarzschild => ("kerr_pass", "schwarzschild_pass"),
                };

                render_graph
                    .remove_slot_edge(input_node_id, "view_entity", old_pass, "view")
                    .expect("Failed to disconnect old pass from subgraph input!");
                render_graph
                    .remove_node_edge(old_pass, "bloom")
                    .expect("Failed to disconnect old pass from bloom pass!");

                render_graph.add_slot_edge(input_node_id, "view_entity", new_pass, "view");
                render_graph.add_node_edge(new_pass, "bloom");
            }
        }
    }
}

#[derive(Component, Clone, ExtractComponent, Default, Debug)]
pub struct MainPassSettings {
    pub updated: bool,
    pub pass: MainPasses,
}

#[derive(Component, Clone, Default, Debug)]
pub enum MainPasses {
    #[default]
    Kerr,
    Schwarzschild,
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

        render_app.add_system(change_state.in_set(RenderSet::Prepare));

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

        render_graph.add_node("kerr_pass", kerr_node);
        render_graph.add_node("schwarzschild_pass", schwarzschild_node);

        render_graph.add_node("bloom", bloom);
        render_graph.add_node("tonemapping", tonemapping);
        render_graph.add_node("fxaa", fxaa);
        render_graph.add_node("ui", ui);
        render_graph.add_node("upscaling", upscaling);
        render_graph.add_slot_edge(input_node_id, "view_entity", "kerr_pass", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "bloom", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "tonemapping", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "fxaa", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "ui", "view");
        render_graph.add_slot_edge(input_node_id, "view_entity", "upscaling", "view");
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
