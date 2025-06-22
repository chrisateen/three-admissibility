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

    pub fn store_maximal_2_packing(&self, v: &mut AdmData)  {
        let p1 = v.p1.clone();
        v.delete_packing();

        for u in p1{
            let u_adm_data = self.adm_data.get(&u).unwrap();
            let should_add_to_pack = u_adm_data.t1.difference(&v.t1).count() >  0;
            if should_add_to_pack {
                v.p1.insert(u);
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
    
    pub fn update_t1_of_v(&mut self, v:&AdmData, u: &mut AdmData) {
        u.t1.remove(&v.id);
        u.n_in_r.insert(v.id);
        
        for x in &v.t1{
            if !u.is_an_endpoint_in_pack(x) {
                u.add_t2_to_packing(x, &v.id);
                return
            }
        }
        
        for x in &v.p1{
            if !u.is_vertex_in_p(x) {
                let t1_x = self.adm_data.get(&x).unwrap().t1.clone();
                for y  in &t1_x {
                    if !u.is_an_endpoint_in_pack(y) {
                        u.add_t3_to_packing(y, x, &v.id);
                        return
                    }
                }
            }
        }
    }
    
    pub fn update_t2_of_v(&mut self, v:&Vertex, u: &mut AdmData) {
        let w = u.remove_v_from_packing(v)[0];
        let w_adm_data = self.adm_data.get(&w).unwrap();
        
        for x in  &w_adm_data.t1 {
            if !u.is_an_endpoint_in_pack(x) {
                u.add_t2_to_packing(x, &w);
                return
            }
        }
        
        for x in &w_adm_data.p1 {
            let x_adm_data = self.adm_data.get(&x).unwrap();
            if !u.is_vertex_in_p(x) {
                for y  in &x_adm_data.t1 {
                    if !u.is_an_endpoint_in_pack(y) {
                        u.add_t3_to_packing(y, &w, x);
                        return
                    }
                }
            }
        }
    }
    
    pub fn update_t3_of_v(&mut self, v:&Vertex, u: &mut AdmData) {
        let p = u.remove_v_from_packing(&v);
        let s = p[0];//p1
        let w = p[1];//p2
        let w_adm_data = self.adm_data.get(&w).unwrap();
        
        for x in  &w_adm_data.t1 {
            if !u.is_an_endpoint_in_pack(x) {
                u.add_t3_to_packing(x, &s, &w);
                return
            }
        }
        
        let s_adm_data = self.adm_data.get(&s).unwrap();
        for x in  &s_adm_data.p1 {
            let x_adm_data = self.adm_data.get(&x).unwrap();
            if !u.is_vertex_in_p(x) {
                for y  in &x_adm_data.t1 {
                    if !u.is_an_endpoint_in_pack(y) {
                        u.add_t3_to_packing(y, &s, x);
                        return
                    }
                }
            }
        }
    }
    
    pub fn update_v(&mut self, v:&mut AdmData) {
        self.candidates.remove(&v.id);
        self.l.remove(&v.id);
        
        let t2_v = self.get_t2_vertices(v);
        let t3_v = self.get_t3_vertices(v);
        
        self.store_maximal_2_packing(v);
        
        for u in &v.t1 {
            let mut u_adm_data = self.adm_data.remove(u).unwrap();
            self.update_t1_of_v(v, &mut u_adm_data);
            self.checks.insert(*u);
            self.adm_data.insert(*u, u_adm_data);
        }
        
        for u in &t2_v{
            let mut u_adm_data = self.adm_data.remove(u).unwrap();
            self.update_t2_of_v(&v.id, &mut u_adm_data);
            self.checks.insert(*u);
            self.adm_data.insert(*u, u_adm_data);
        }
        
        for  u in &t3_v{
            let mut u_adm_data = self.adm_data.remove(u).unwrap();
            self.update_t3_of_v(&v.id, &mut u_adm_data);
            self.checks.insert(*u);
            self.adm_data.insert(*u, u_adm_data);
        }
    }
    
    pub fn get_next_v_in_ordering(&mut self, p: usize) -> Option<Vertex>{
        let v = self.candidates.iter().next();
        match v { 
            Some(&v) => {
                let mut v_adm_data = self.adm_data.remove(&v).unwrap();
                self.update_v(&mut v_adm_data);
                self.adm_data.insert(v, v_adm_data);
                Some(v)
            }
            None => None
        }
        
    }
}