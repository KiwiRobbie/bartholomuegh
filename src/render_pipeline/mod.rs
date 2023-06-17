use std::option::Iter;
use std::sync::mpsc;
use std::time::Duration;

use self::kerr_pass::{KerrPassNode, KerrPassPlugin};
use self::schwarzschild_pass::{SchwarzschildPassNode, SchwarzschildPassPlugin};
use bevy::app::{AppLabel, AppLabelId, SubApp};
use bevy::core_pipeline::bloom::BloomSettings;
use bevy::ecs::schedule;
use bevy::render::extract_component::{ExtractComponent, ExtractComponentPlugin};
use bevy::render::RenderSet;
use bevy::utils::tracing::instrument::WithSubscriber;
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

fn setup_ker(graph: ResMut<RenderGraph>) {
    dbg!(graph);
    // dbg!(graph);
}

fn setup_schwarzchild(graph: ResMut<RenderGraph>) {
    dbg!(graph);
    // dbg!(graph);
}

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
                let voxel_graph = graph
                    .get_sub_graph_mut("voxel")
                    .expect("Could not get subgraph to switch!");
                // let input_node_id = voxel_graph
                //     .get_input_node()
                //     .expect("Could not obtain input node!")
                //     .id;

                let input_node_id = voxel_graph.get_input_node().unwrap().id;
                match settings.pass {
                    MainPasses::Kerr => {
                        dbg!("Switching to ker");
                        voxel_graph.remove_slot_edge(
                            input_node_id,
                            "view_entity",
                            "schwarzschild_pass",
                            "view",
                        );
                        voxel_graph.remove_node_edge("schwarzschild_pass", "bloom");

                        voxel_graph.add_slot_edge(
                            input_node_id,
                            "view_entity",
                            "kerr_pass",
                            "view",
                        );
                        voxel_graph.add_node_edge("kerr_pass", "bloom");
                    }
                    MainPasses::Schwarzschild => {
                        dbg!("Switching to schwarzchild");
                        voxel_graph.remove_slot_edge(
                            input_node_id,
                            "view_entity",
                            "kerr_pass",
                            "view",
                        );
                        voxel_graph.remove_node_edge("kerr_pass", "bloom");

                        voxel_graph.add_slot_edge(
                            input_node_id,
                            "view_entity",
                            "schwarzschild_pass",
                            "view",
                        );
                        voxel_graph.add_node_edge("schwarzschild_pass", "bloom");
                    }
                }
            }
        }
    }

    // dbg!(world.components());
    // // app.world.resource_mut::<RenderGraph>();
    // dbg!(world.id());
    // // dbg!(App::from_world(world).world.id());

    // let render_app = match app.get_sub_app_mut(RenderApp) {
    //     Ok(render_app) => render_app,
    //     Err(_) => return,
    // };

    // dbg!(&render_app);
    // let sub_app = app.sub_app_mut(RenderApp);
    // let sub_world: &mut World = &mut app.world;

    // NextState::<RenderPassState>::from_world(sub_world).set(RenderPassState::Schwarzschild);
    // schedule::apply_state_transition::<RenderPassState>(sub_world);
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

        // build voxel render graph
        let mut voxel_graph = RenderGraph::default();
        let input_node_id =
            voxel_graph.set_input(vec![SlotInfo::new("view_entity", SlotType::Entity)]);

        // render graph
        let kerr_node = KerrPassNode::new(&mut render_app.world);
        let schwarzschild_node = SchwarzschildPassNode::new(&mut render_app.world);

        let bloom = BloomNode::new(&mut render_app.world);
        let tonemapping = TonemappingNode::new(&mut render_app.world);
        let fxaa = FxaaNode::new(&mut render_app.world);
        let ui = UiPassNode::new(&mut render_app.world);
        let upscaling = UpscalingNode::new(&mut render_app.world);

        voxel_graph.add_node("kerr_pass", kerr_node);
        voxel_graph.add_node("schwarzschild_pass", schwarzschild_node);

        voxel_graph.add_node("bloom", bloom);
        voxel_graph.add_node("tonemapping", tonemapping);
        voxel_graph.add_node("fxaa", fxaa);
        voxel_graph.add_node("ui", ui);
        voxel_graph.add_node("upscaling", upscaling);
        voxel_graph.add_slot_edge(input_node_id, "view_entity", "kerr_pass", "view");
        voxel_graph.add_slot_edge(input_node_id, "view_entity", "bloom", "view");
        voxel_graph.add_slot_edge(input_node_id, "view_entity", "tonemapping", "view");
        voxel_graph.add_slot_edge(input_node_id, "view_entity", "fxaa", "view");
        voxel_graph.add_slot_edge(input_node_id, "view_entity", "ui", "view");
        voxel_graph.add_slot_edge(input_node_id, "view_entity", "upscaling", "view");
        voxel_graph.add_node_edge("kerr_pass", "bloom");
        voxel_graph.add_node_edge("bloom", "tonemapping");
        voxel_graph.add_node_edge("tonemapping", "fxaa");
        voxel_graph.add_node_edge("fxaa", "ui");
        voxel_graph.add_node_edge("ui", "upscaling");

        // insert the voxel graph into the main render graph
        let mut graph = render_app.world.resource_mut::<RenderGraph>();
        graph.add_sub_graph("voxel", voxel_graph);

        // app.world.add_schedule(Timer::from_seconds(10.0, TimerMode::Once)., label)
        // .add_system()

        // let s = OnEnter::<RenderPassState>(RenderPassState::Kerr).with_subscriber(setup_ker);

        // .add_systems(SystemSet::on_enter(RenderPassState::Kerr).with_system(setup_ker))
        // .add_systems(
        // SystemSet::on_enter(RenderPassState::Schwarzschild).with_system(setup_schwarzchild),
        // );

        // render_app
        //     .world
        //     .insert_resource(Events::<SwitchPassEvent>::default());

        // // Create a schedule to store our systems
        // let mut schedule = Schedule::default();

        // // Events need to be updated in every frame in order to clear our buffers.
        // // This update should happen before we use the events.
        // // Here, we use system sets to control the ordering.
        // #[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
        // pub struct FlushEvents;

        // schedule.add_system(Events::<SwitchPassEvent>::update_system.in_set(FlushEvents));
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
