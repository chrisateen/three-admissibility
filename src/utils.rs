use flate2::Compression;
use flate2::write::GzEncoder;
use graphbench::editgraph::EditGraph;
use graphbench::graph::Vertex;
use graphbench::io::LoadFromFile;
use std::io::{BufRead, Write};
use std::path::PathBuf;

pub fn load_graph(network_path: String, network: &String) -> EditGraph {
    let file_dir = format!("{}/{}.txt.gz", network_path, network);
    EditGraph::from_gzipped(&file_dir)
        .unwrap_or_else(|_| panic!("Error occurred loading graph {}", network))
}

pub fn save_ordering_to_file(path: String, network: String, order: Vec<Vertex>) {
    let folder = PathBuf::from(path);
    std::fs::create_dir_all(&folder).unwrap();
    let file_path = folder.join(network.as_str().to_owned() + ".txt.gz");

    let file = std::fs::File::create(file_path).unwrap();
    let mut gz = GzEncoder::new(file, Compression::default());

    for v in order {
        writeln!(gz, "{}", v).unwrap();
    }

    gz.finish().unwrap();
}
