use std::collections::HashSet;
use graphbench::editgraph::EditGraph;
use graphbench::graph::{Graph, Vertex, VertexMap, VertexSet};
use crate::admData::AdmData;

pub struct AdmGraph {
    l: VertexSet,
    r: VertexSet,
    checks: VertexSet,
    pub candidates: VertexSet,
    adm_data: VertexMap<AdmData>,
}

impl AdmGraph {
    fn new(graph: &EditGraph) -> Self {
        let mut adm_data = VertexMap::default();
        let l = graph.vertices().copied().collect();
        for u in graph.vertices() {
            let adm_vertex = AdmData::new(*u, graph.neighbours(u).copied().collect());
            adm_data.insert(*u, adm_vertex);
        }
        AdmGraph {
            l,
            r: VertexSet::default(),
            checks: VertexSet::default(),
            candidates: VertexSet::default(),
            adm_data,
        }
    }

    /*
    Initialises candidates with vertices with degree <= p
    */
    pub fn initialise_candidates(&mut self, p: usize) {
        for (u, adm_data) in &self.adm_data {
            if adm_data.t1.len() <= p {
                self.candidates.insert(*u);
            }
        }
    }

    /*
    Computes vias and a maximal 2-packing for a vertex being moved to R
    */
    fn compute_vias(&mut self, p: usize, v: &mut AdmData){
        let mut counter : VertexMap<i32> = VertexMap::default();
        let mut t2_u = VertexSet::default();
        let t1_v: HashSet<Vertex> = v.t1.intersection(&self.l).into_iter().copied().collect();
        let n_in_r: Vec<Vertex> = v.n_in_r.iter().copied().collect();

        v.delete_packing(); //clear 3-packing of v as we now want to store a 2-packing for v

        for u in n_in_r {
            let mut u_adm_data = self.adm_data.remove(&u).unwrap();

            for w in u_adm_data.t1.intersection(&self.l) {
                if *w == v.id{
                    continue;
                }

                let num_vias_w = counter.entry(*w).or_insert(0);
                if !t2_u.contains(w){
                    t2_u.insert(*w);
                    v.vias.insert(u);
                    *num_vias_w += 1;
                    //If this is the first time we encounter a vertex in t2_l of v
                    //need to check if we can add w,u as a path in the 2-packing for v
                    if !t1_v.contains(w)  & v.can_add_t2_path_to_pack(w, &u) {
                        v.add_t2_to_packing(w, &u);
                    }
                }else {
                    if num_vias_w < &mut ((2 * p as i32) + 1) {
                        v.vias.insert(u);
                        *num_vias_w += 1;
                    }
                }
            }
            self.adm_data.insert(u, u_adm_data);
        }
    }

    /*
        Fetches all the vertices in L that is in T1, T2 and T3 of v
    */
    fn identify(&self, v: &AdmData)-> VertexSet{
        let mut t: VertexSet = v.t1.intersection(&self.l).cloned().collect();

        for u in v.t1.intersection(&self.r){
            let u_adm_data = self.adm_data.get(&u).unwrap();

            //identify t2_u
            for w in u_adm_data.t1.intersection(&self.l) {
                t.insert(*w);
            }

            //identify t3
            for w in u_adm_data.t1.intersection(&self.r){
                if w != &v.id {
                    let w_adm_data = self.adm_data.get(&w).unwrap();
                    for y in w_adm_data.t1.intersection(&self.l){
                        t.insert(*y);
                    }
                }
            }
        }
        t
    }

    /*
        Simple update of a packing
    */
    fn simple_update(&mut self, u: &mut AdmData, v: Vertex){
        if !u.is_an_endpoint_in_pack(&v){
            return;
        }

        let path: Vec<Vertex> = u.remove_v_from_packing(&v);
        let w = if path.len() > 0 { path[0] } else { v };
        let w_adm_data = self.adm_data.get(&w).unwrap();

        //check if we can add a path of length 2
        for x in w_adm_data.t1.intersection(&self.l){
            if !u.is_an_endpoint_in_pack(x){
                u.add_t2_to_packing(x, &w);
                return;
            }
        }

        //check if we can add a path of length 3
        //TODO is there an error in the pseudocode?
        for x in w_adm_data.t1.intersection(&self.r){
            if u.t2.contains(x) | u.t3.contains(x){
                continue;
            }
            let x_adm_data = self.adm_data.get(&x).unwrap();
            for y in x_adm_data.t1.intersection(&self.l){
                if !u.is_an_endpoint_in_pack(&y){
                    u.add_t3_to_packing(&y, x, &w);
                }
            }
        }
    }

    fn maximal_update(&mut self, u: &mut AdmData){
        for w in u.n_in_r.difference(&u.t1){
            let w_adm_data = self.adm_data.get(&w).unwrap();
            for x in &w_adm_data.vias{
                //TODO
            }
        }
    }

    //fn maximum_update(&mut self, v: &mut AdmData){}

    fn update(&mut self, p: usize, v:&Vertex){
        let mut v_adm_data = self.adm_data.remove(&v).unwrap();
        let t = self.identify(&v_adm_data);
        self.compute_vias(p, &mut v_adm_data);
        self.adm_data.insert(*v, v_adm_data);

        for u in t{
            let mut u_adm_data = self.adm_data.remove(&u).unwrap();
            self.simple_update(&mut u_adm_data, *v);
            if u_adm_data.size_of_packing() <= p{
                self.maximal_update(&mut u_adm_data);
            }
            self.adm_data.insert(u, u_adm_data);
        }
    }

    pub fn get_next_v_in_ordering(&mut self, p: usize) -> Option<Vertex>{
        let v = self.candidates.iter().next();

        match v {
            Some(&v) => {
                self.candidates.remove(&v);
                self.l.remove(&v);
                self.r.insert(v);
                //TODO updates
                Some(v)
            }
            None => None
        }

    }
}