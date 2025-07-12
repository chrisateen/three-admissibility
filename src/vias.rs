use graphbench::graph::{Vertex, VertexMap, VertexSet};

pub struct Vias {
    //Key vias[v_in_r][v_in_t2_l_of_v]
    pub vias: VertexMap<VertexMap<VertexSet>>,
    p: usize
}

impl Vias {
    pub fn new(p : usize) -> Self {
        Vias {
            vias: VertexMap::default(),
            p
        }
    }

    pub fn add_a_via(&mut self, v:  Vertex, t2_v: Vertex, via: Vertex) -> bool {
        let v_entries = self.vias.entry(v).or_insert(VertexMap::default());
        let t2_v_entries = v_entries.entry(t2_v).or_insert(VertexSet::default());
        
        let max_vias = (2*self.p + 1);
        
        if t2_v_entries.len() > max_vias{
            panic!("Number of vias for {} is too large", t2_v);
        }
        if t2_v_entries.len() == max_vias{
            return false;
        }
        t2_v_entries.insert(via);
        true
    }

    //TODO to delete vias for a vertex v moving to R, we need to identify all the t2 of v
    // and for each of them delete v from vias
    pub fn remove_vias(&mut self, v: Vertex) {

    }

    pub fn get_vias(&self, v:  Vertex, t2_v: Vertex) -> Option<&VertexSet> {
        let v_entries = self.vias.get(&v);
        match v_entries {
            Some(v_entries) => {
                let vias = v_entries.get(&t2_v);
                match vias {
                    Some(vias) => Some(vias),
                    None => None
                }
            }
            None => None
        }
    }
}