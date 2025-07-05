use graphbench::graph::{Vertex, VertexMap, VertexSet};

pub struct AdmData {
    pub id: Vertex,
    pub n_in_r: VertexSet,
    pub t1: VertexSet,
    pub t2: VertexSet,
    pub t3: VertexSet,
    pub vias: VertexSet,
    pub packing: VertexMap<Vec<Vertex>>,
}

impl AdmData {
    pub fn new(v: Vertex, v_neighbours: VertexSet) -> Self {
        AdmData {
            id: v,
            n_in_r: VertexSet::default(),
            t1: v_neighbours,
            t2: VertexSet::default(),
            t3: VertexSet::default(),
            vias: VertexSet::default(),
            packing: VertexMap::default(),
        }
    }

    pub fn can_add_t2_path_to_pack(&mut self, t2: &Vertex, t1: &Vertex) -> bool {
        !self.packing.contains_key(t2) & !self.t1.contains(t1)
    }

    pub fn can_add_t3_path_to_pack(&mut self, t3: &Vertex, t2: &Vertex, t1: &Vertex) -> bool {
        !self.packing.contains_key(t3) & !self.t1.contains(t1) & !self.t2.contains(t2)
    }
    pub fn add_t2_to_packing(&mut self, t2: &Vertex, t1: &Vertex) {
        self.packing.insert(*t2, vec![*t1]);
        self.t1.insert(*t1);
        self.t2.insert(*t2);
    }

    pub fn add_t3_to_packing(&mut self, t3: &Vertex, t1: &Vertex,  t2: &Vertex) {
        self.packing.insert(*t3, vec![*t1, *t2]);
        self.t1.insert(*t1);
        self.t2.insert(*t2);
        self.t3.insert(*t3);
    }

    pub fn delete_packing(&mut self) {
        self.packing.clear();
        self.t2.clear();
        self.t3.clear();
    }

    pub fn is_an_endpoint_in_pack(&self, v: &Vertex) -> bool {
        (self.t1.contains(v) & !self.n_in_r.contains(v)) | self.packing.contains_key(v)
    }

    pub fn remove_v_from_packing(&mut self, v: &Vertex) ->  Vec<Vertex> {
        if (self.t1.contains(v) & !self.n_in_r.contains(v)){
            self.t1.remove(v);
            self.n_in_r.insert(*v);
            return vec![];
        }

        let p = self.packing.remove(v).unwrap();
        self.t1.remove(&p[0]);
        if p.len() == 1 {
            self.t2.remove(v);
        }else {
            self.t2.remove(&p[1]);
            self.t3.remove(v);
        }
        p
    }

    pub fn size_of_packing(&self) -> usize{
        let num_t1_l = self.t1.difference(&self.n_in_r).count();
        num_t1_l + self.packing.len()
    }
}