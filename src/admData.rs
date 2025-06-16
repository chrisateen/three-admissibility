use graphbench::graph::{Vertex, VertexMap, VertexSet};

pub struct AdmData {
    pub id: Vertex,
    pub t1: VertexSet,
    pub p1: VertexSet,
    pub p2: VertexSet,
    pub packing: VertexMap<Vec<Vertex>>,
}

impl AdmData {
    pub fn new(v: Vertex, v_neighbours: VertexSet) -> Self {
        AdmData {
            id: v,
            t1: v_neighbours,
            p1: VertexSet::default(),
            p2: VertexSet::default(),
            packing: VertexMap::default(),
        }
    }

    pub fn is_size_of_packing_p(&self, p: usize) -> bool{
        self.t1.len() & self.packing.len() <= p
    }
    
    pub fn delete_packing(&mut self) {
        self.p1.clear();
        self.p2.clear();
        self.packing.clear();
    }
}