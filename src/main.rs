use std::cmp::max;
use clap::{Parser, Subcommand};
use graphbench::editgraph::EditGraph;
use graphbench::graph::{Graph, MutableGraph, Vertex};
use graphbench::io::LoadFromFile;
use peak_alloc::PeakAlloc;
use crate::adm_graph::AdmGraph;
use crate::utils::{load_graph, save_ordering_to_file};

mod utils;
mod adm_data;
mod adm_graph;
mod vias;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
struct Args {
    /// network file name
    network: String,

    /// start p value
    p: i32,

    /// Path to network
    network_path: String,

    #[clap(short, long, default_value_t = false)]
    /// Whether to track memory consumption
    track_memory: bool,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Whether to save ordering to file
    Save {
        /// The path to save ordering to
        #[arg(default_value= "results")]
        path: String,
    },
}

fn compute_ordering(p: usize, graph: &EditGraph, save_order: bool) -> Option<Vec<Vertex>> {
    let mut adm_graph = AdmGraph::new(graph, p);

    adm_graph.initialise_candidates();

    let mut next_vertex = adm_graph.get_next_v_in_ordering();
    let mut order = Vec::default();
    while next_vertex.is_some() && !adm_graph.is_all_vertices_in_r_or_candidates() {
        let v = next_vertex.unwrap();
        if save_order {
            order.push(v);
        }
        next_vertex = adm_graph.get_next_v_in_ordering();
    }
    if save_order {
        order.extend(next_vertex.iter()); // Adds vertex if not None
    }

    let found_order = adm_graph.is_all_vertices_in_r_or_candidates();

    if found_order {
        if save_order {
            order.extend(adm_graph.candidates.iter());
            assert_eq!(order.len(), graph.num_vertices());
        }
        Some(order)
    } else {
        None
    }
}

fn next_p_value(p: i32, is_p: bool, lowest_p: i32, highest_not_p: i32) -> i32 {
    //Stop where the lowest p is p or the highest p + 1 is p
    if (p - highest_not_p <= 1 && is_p) || (p - lowest_p).abs() == 1 {
        return -1;
    }
    //Continue to double the p value we check if we haven't found a value where G is p,2 admissible
    if lowest_p == -1 && !is_p {
        return (p * 2) as i32;
    }
    //Once we found a p value keep halving the search between the lowest p and the highest not p
    let x = max(p, lowest_p);
    return (x + highest_not_p) / 2;
}

fn main() {
    let args = Args::parse();

    let network_path = args.network_path;
    let network = args.network;
    let mut p = args.p;

    let track_memory = args.track_memory;

    let save_path = match args.command {
        None => None,
        Some(Commands::Save { path }) => Some(path)
    };

    let mut lowest_p: i32 = -1;
    let mut highest_not_p: i32 = -1;
    let mut best_order = None;

    let mut peak_mem : f32;

    let mut graph = load_graph(network_path, &network);

    graph.remove_loops();

    if track_memory{
        peak_mem = PEAK_ALLOC.peak_usage_as_kb();
        println!("Max memory used after graph loading in kb is {}", peak_mem);
    }

    loop {
        println!("checking for p {}", p);
        let result = compute_ordering(p as usize, &graph, save_path.is_some());
        let mut found_better = false;
        if let Some(order) = result {
            assert!(lowest_p == -1 || p < lowest_p);
            lowest_p = p;
            best_order = Some(order);
            found_better = true;
        } else {
            assert!(p > highest_not_p);
            highest_not_p = p;
        }

        let next_p = next_p_value(p, found_better, lowest_p, highest_not_p);
        if next_p == -1 {
            if !found_better {
                p = lowest_p;
            }
            break;
        }
        p = next_p;
    }

    println!("p is {}", p);

    if track_memory {
        peak_mem = PEAK_ALLOC.peak_usage_as_kb();
        println!("Max memory used in total kb is {}", peak_mem);
    }

    match save_path {
        None => {}
        Some(path) => {
            if let Some(order) = best_order {
                save_ordering_to_file(path, network, order);
            }
        }
    }
}
