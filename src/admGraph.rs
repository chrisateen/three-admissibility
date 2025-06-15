use graphbench::editgraph::EditGraph;
use graphbench::graph::{Graph, VertexMap, VertexSet, Vertex};
use crate::admData::AdmData;

pub struct AdmGraph {
    l: VertexSet,
    checks: VertexSet,
    pub candidates: VertexSet,
    adm_data: VertexMap<AdmData>,
}

impl AdmGraph {
    pub fn new(graph: &EditGraph) -> Self {
        let mut adm_data = VertexMap::default();
        let l = graph.vertices().copied().collect();
        for u in graph.vertices() {
            let adm_vertex = AdmData::new(*u, graph.neighbours(u).copied().collect());
            adm_data.insert(*u, adm_vertex);
        }
        AdmGraph {
            l,
            checks: VertexSet::default(),
            candidates: VertexSet::default(),
            adm_data,
        }
    }

    pub fn initialise_candidates(&mut self, p: usize) {
        for (u, adm_data) in &self.adm_data {
            if adm_data.t1.len() <= p {
                self.candidates.insert(*u);
            }
        }
    }
    
    pub fn get_next_v_in_ordering(&mut self, p: usize) -> Option<Vertex>{
        let v = self.candidates.iter().next();
        match v { 
            Some(&v) => {
                self.candidates.remove(&v);
                Some(v)
            }
            None => None
        }
        
    }
}