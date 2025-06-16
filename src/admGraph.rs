use std::collections::HashSet;
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
    
    pub fn get_t2_vertices(&self, v: &AdmData) -> VertexSet {
        let mut t2_vertices : VertexSet = VertexSet::default();
        for u in &v.p1 {
            let u_t1_vertices = self.adm_data.get(&u).unwrap().t1.clone();
            for x  in u_t1_vertices.difference(&v.t1) {
                if x != &v.id{
                    t2_vertices.insert(*x);
                }
            }
        }
        t2_vertices
    }
    
    pub fn get_t3_vertices(&self, v: &AdmData) -> VertexSet {
        let mut t3_vertices = VertexSet::default();
        for u in &v.p1 {
            let u_adm_data = self.adm_data.get(&u).unwrap();
            let t2_from_u = self.get_t2_vertices(u_adm_data);
            for x  in t2_from_u.difference(&v.t1) {
                if x != &v.id{
                    t3_vertices.insert(*x);
                }
            }
        }
        t3_vertices
    }

    pub fn store_maximal_2_packing(&self, v: &mut AdmData)  {
        let p1 = v.p1.clone();
        v.delete_packing();
        
        for u in p1{
            let u_adm_data = self.adm_data.get(&u).unwrap();
            let should_add_to_pack = u_adm_data.t1.difference(&v.t1).count() >  0;
            if should_add_to_pack {
                //TODO think about how 2 packing will be stored
                v.p1.insert(u);
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