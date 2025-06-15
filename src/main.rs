use clap::{Parser, Subcommand};
use graphbench::io::LoadFromFile;
use peak_alloc::PeakAlloc;

mod admData;
mod admGraph;
mod utils;

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
}