use graphbench::graph::{Vertex, VertexMap, VertexSet};

pub struct AdmData {
    pub id: Vertex,
    pub t1: VertexSet,
    pub n_in_r: VertexSet,
    pub p1: VertexSet,
    pub p2: VertexSet,
    pub vias: VertexSet,
    pub packing: VertexMap<Vec<Vertex>>,
}

impl AdmData {
    pub fn new(v: Vertex, v_neighbours: VertexSet) -> Self {
        AdmData {
            id: v,
            n_in_r: VertexSet::default(),
            t1: v_neighbours,
            p1: VertexSet::default(),
            p2: VertexSet::default(),
            vias: VertexSet::default(),
            packing: VertexMap::default(),
        }
    }

    pub fn size_of_packing(&self) -> usize{
        self.t1.len() + self.packing.len()
    }
    
    pub fn remove_v_from_packing(&mut self, v: &Vertex) ->  Vec<Vertex> {
        let p = self.packing.remove(v).unwrap();
        self.p1.remove(&p[0]);
        if p.len() == 2 {
            self.t1.remove(&p[1]);
        }
        p
    }
    
    pub fn is_an_endpoint_in_pack(&self, v: &Vertex) -> bool {
        self.t1.contains(v) | self.packing.contains_key(v)
    }
    
    pub fn add_t2_to_packing(&mut self, t2: &Vertex, p1: &Vertex) {
        self.packing.insert(*t2, vec![*p1]);
        self.p1.insert(*p1);
    }
    
    pub fn add_t3_to_packing(&mut self, t3: &Vertex, p1: &Vertex,  p2: &Vertex) {
        self.packing.insert(*t3, vec![*p1, *p2]);
        self.p1.insert(*p1);
        self.p2.insert(*p2);
    }
    
    pub fn is_vertex_in_p(&self, v: &Vertex) -> bool {
        self.p1.contains(v) | self.p2.contains(v)
    }

    pub fn get_all_t2_vertices(&self) -> VertexSet {
        self.packing.iter()
            .filter(|(_, p)| p.len() == 1)
            .map(|(v, _)| *v)
            .collect()
    }
    
    pub fn delete_packing(&mut self) {
        self.p1.clear();
        self.p2.clear();
        self.packing.clear();
    }
}