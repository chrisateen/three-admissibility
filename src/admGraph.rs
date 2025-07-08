use std::collections::HashSet;
use graphbench::editgraph::EditGraph;
use graphbench::graph::{Graph, Vertex, VertexMap, VertexSet};
use crate::admData::AdmData;

pub struct AdmGraph {
    l: VertexSet,
    r: VertexSet,
    pub candidates: VertexSet,
    adm_data: VertexMap<AdmData>,
    num_of_vertices: usize,
    p: usize,
}

impl AdmGraph {
    pub(crate) fn new(graph: &EditGraph, p: usize) -> Self {
        let mut adm_data = VertexMap::default();
        let l = graph.vertices().copied().collect();
        for u in graph.vertices() {
            let adm_vertex = AdmData::new(*u, graph.neighbours(u).copied().collect());
            adm_data.insert(*u, adm_vertex);
        }
        AdmGraph {
            l,
            r: VertexSet::default(),
            candidates: VertexSet::default(),
            adm_data,
            num_of_vertices: graph.num_vertices(),
            p,
        }
    }

    /*
    Initialises candidates with vertices with degree <= p
    */
    pub fn initialise_candidates(&mut self) {
        for (u, adm_data) in &self.adm_data {
            if adm_data.t1.len() <= self.p {
                self.candidates.insert(*u);
            }
        }
    }

    /*
    Computes vias and a maximal 2-packing for a vertex being moved to R
    */
    fn compute_vias(&mut self, v: &mut AdmData){
        let mut counter : VertexMap<i32> = VertexMap::default();
        let mut t2_u = VertexSet::default();
        let t1_v: HashSet<Vertex> = v.t1.intersection(&self.l).into_iter().copied().collect();
        let n_in_r: Vec<Vertex> = v.n_in_r.iter().copied().collect();

        v.delete_packing(); //clear 3-packing of v as we now want to store a 2-packing for v

        for u in n_in_r {
            let u_adm_data = self.adm_data.remove(&u).unwrap();

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
                    if num_vias_w < &mut ((2 * self.p as i32) + 1) {
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
        for x in w_adm_data.t1.intersection(&self.r){
            if u.t2.contains(x) | u.t3.contains(x){
                continue;
            }
            let x_adm_data = self.adm_data.get(&x).unwrap();
            for y in x_adm_data.t1.intersection(&self.l){
                if !u.is_an_endpoint_in_pack(&y){
                    //Check if there is a shorter path to y
                    //TODO use editgraph to check N(u)
                    if u.n_in_r.contains(&x) {
                        //as there is a shorter path to y
                        //check if there is another path of length 3 through w
                        u.add_t2_to_packing(&y,x);
                        break;
                    }
                    u.add_t3_to_packing(&y, x, &w);
                }
            }
        }
    }

    /*
        Try to see if a disjoint path can be added to the packing of u
     */
    fn stage_1_update(&mut self, u: &mut AdmData){
        let n_r_u_not_in_pack: VertexSet = u.n_in_r.difference(&u.t1).cloned().collect();
        for w in n_r_u_not_in_pack{
            let w_adm_data = self.adm_data.get(&w).unwrap();
            //After w was moved to R some vertices in t1_l of u may have moved to r
            let t1_to_check: VertexSet = w_adm_data.t1.intersection(&self.r).cloned().collect();
            for x in w_adm_data.vias.union(&t1_to_check) {
                let x_adm_data = self.adm_data.get(&x).unwrap();
                for y in x_adm_data.t1.intersection(&self.l){
                    if !u.is_an_endpoint_in_pack(&y){
                        //TODO use editgraph to check N(u)
                        if u.n_in_r.contains(&x) {
                            //If there is a shorter path to y add it
                            u.add_t2_to_packing(&y,x);
                        }else {
                            u.add_t3_to_packing(&y, x, &w);
                        }
                        return;
                    }
                }
            }
        }
    }

    /*
        Find an augmenting path to see if packing of u can be extended
     */
    fn stage_2_update(&mut self, u: &mut AdmData){}


    /*
        Update data structures now that vertx v is moving to R
    */
    fn update(&mut self, v:&Vertex){
        let mut v_adm_data = self.adm_data.remove(&v).unwrap();
        let t = self.identify(&v_adm_data);
        self.compute_vias(&mut v_adm_data);
        self.adm_data.insert(*v, v_adm_data);

        for u in t{
            if self.candidates.contains(&u){
                continue;
            }
            let mut u_adm_data = self.adm_data.remove(&u).unwrap();
            self.simple_update(&mut u_adm_data, *v);

            //check if a disjoint path can be added to packing of u
            if u_adm_data.size_of_packing() <= self.p{
                self.stage_1_update(&mut u_adm_data);
            }
            //If size of packing is still p after simple update
            //then do augmenting path to see if packing of u can be extended
            if u_adm_data.size_of_packing() <= self.p{
                self.stage_2_update(&mut u_adm_data);
            }

            //If size of packing is still p then add to candidates
            if u_adm_data.size_of_packing() <= self.p{
                self.candidates.insert(u);
            }
            self.adm_data.insert(u, u_adm_data);
        }
    }

    pub fn is_all_vertices_in_r_or_candidates(&self) -> bool {
        return self.r.len() + self.candidates.len() == self.num_of_vertices;
    }

    pub fn get_next_v_in_ordering(&mut self) -> Option<Vertex>{
        let v = self.candidates.iter().next();

        match v {
            Some(&v) => {
                self.candidates.remove(&v);
                self.l.remove(&v);
                self.r.insert(v);
                self.update(&v);
                Some(v)
            }
            None => None
        }

    }
}