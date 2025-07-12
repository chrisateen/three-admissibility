use graphbench::editgraph::EditGraph;
use graphbench::graph::{Graph, Vertex, VertexMap, VertexSet};
use crate::adm_data::AdmData;
use crate::vias::Vias;

pub(crate) struct AdmGraph<'a> {
    l: VertexSet,
    r: VertexSet,
    pub candidates: VertexSet,
    adm_data: VertexMap<AdmData>,
    num_of_vertices: usize,
    p: usize,
    graph: &'a EditGraph,
    vias: Vias
}

impl<'a> AdmGraph<'a> {
    pub(crate) fn new(graph: &'a EditGraph, p: usize) -> Self {
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
            graph,
            vias: Vias::new(p),
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
        v.delete_packing(); //clear 3-packing of v as we now want to store a 2-packing for v

        for u in  v.n_in_r.clone() {
            let u_adm_data = self.adm_data.get(&u).unwrap();

            for w in u_adm_data.t1.difference(&self.l) {
                if !v.is_v_in_pack(&w){
                    v.add_t2_to_packing(&w, &u);
                }
                self.vias.add_a_via(v.id, *w, u);
            }

            for w in v.t1.difference(&self.l) {
                // note that as v is moving to R
                // from u we can get to w via v
                self.vias.add_a_via(u, *w, v.id);
            }
        }
    }

    /*
        Fetches all the vertices in L that is in T1, T2 and T3 of v
    */
    fn collect_targets(&self, v: &AdmData) -> VertexSet{
        let mut t: VertexSet = v.t1.intersection(&self.l).cloned().collect();
        let t1_t2: VertexSet  = v.t1.union(&v.t2).cloned().collect();

        for u in t1_t2.intersection(&self.r) {
            if *u == v.id { continue;}
            let u_adm_data = self.adm_data.get(&u).unwrap();

            for w in u_adm_data.t1.intersection(&self.l) {
                t.insert(*w);
            }

            if !self.graph.adjacent(&u, &v.id){
                continue;
            };

            //if u in N(v) then we need to use 2-packing of u to get vertices t3 of v
            for w in u_adm_data.t1.intersection(&self.r) {
                if *w == v.id { continue;}
                let w_adm_data = self.adm_data.get(&w).unwrap();

                for x in w_adm_data.t1.intersection(&self.l) {
                    t.insert(*x);
                }
            }
        }
        t
    }

    /*
        Simple update of a packing
    */
    fn simple_update(&mut self, u: &mut AdmData, v: Vertex){
        if !u.is_v_in_pack(&v){
            return;
        }

        let path: Vec<Vertex> = u.remove_v_from_packing(&v);
        let w = if path.len() > 0 { path[0] } else { v };
        let w_adm_data = self.adm_data.get(&w).unwrap();

        //check if we can add a path of length 2
        for x in w_adm_data.t1.intersection(&self.l){
            if !u.is_v_in_pack(x){
                u.add_t2_to_packing(x, &w);
                return;
            }
        }

        //check if we can add a path of length 3
        for x in w_adm_data.t1.intersection(&self.r){
            if u.t2.contains(x) | u.t3.contains(x) | (u.id == v) {
                continue;
            }
            let x_adm_data = self.adm_data.get(&x).unwrap();
            for y in x_adm_data.t1.intersection(&self.l){
                if !u.is_v_in_pack(&y){
                    //Check if there is a shorter path to y (i.e u,x,y)
                    // if so add the shorter path instead
                    if self.graph.adjacent(&u.id, x) {
                        u.add_t2_to_packing(&y,x);
                    }else {
                        u.add_t3_to_packing(&y, x, &w);
                    }
                    return;
                }
            }
        }
    }

    /*
        Try to see if a disjoint path can be added to the packing of u
     */
    fn stage_1_update(&mut self, u: &mut AdmData){
        let t1 : VertexSet = u.t1.intersection(&self.l).cloned().collect();
        let t3_and_t2 : VertexSet = self.collect_targets(u).difference(&t1).cloned().collect();

        for w in u.n_in_r.clone(){
            if u.is_v_in_pack(&w){continue;}

            for y in &t3_and_t2{
                if u.is_v_in_pack(&y){continue;}
                let w_adm_data = self.adm_data.get(&w).unwrap();

                //first check if w,y is a path
                if w_adm_data.t1.contains(y){
                    u.add_t2_to_packing(y, &w);
                    return;
                }

                //if w,y is not a path check vias between w and y
                let x_vias = self.vias.get_vias(w, *y);
                if x_vias == None { continue; }

                for x in x_vias.unwrap(){
                    if u.is_v_in_pack(x){ continue;}

                    //first check if there is a shorter path x,y
                    if self.graph.adjacent(&u.id, x) {
                        u.add_t2_to_packing(y,x);
                    }else {
                        u.add_t3_to_packing(y,x, &w);
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
        let t = self.collect_targets(&v_adm_data);
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

    /*
        Checks if all the vertices are in r or can be added to r
     */
    pub fn is_all_vertices_in_r_or_candidates(&self) -> bool {
        self.r.len() + self.candidates.len() == self.num_of_vertices
    }

    /*
        Gets next vertex in the ordering
     */
    pub fn get_next_v_in_ordering(&mut self) -> Option<Vertex>{
        let v = self.candidates.iter().next();;

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

#[cfg(test)]
mod test_adm_graph {
    use super::*;
    use graphbench::editgraph::EditGraph;
    use graphbench::graph::{EdgeSet, MutableGraph};

    fn create_test_graph(edges: EdgeSet) -> EditGraph {
        let mut graph = EditGraph::new();
        for (u, v) in edges.iter() {
            graph.add_edge(u, v);
        }

        graph
    }

    #[test]
    fn initialise_candidates_should_add_vertices_with_degree_p_or_less_to_candidates() {
        let edges: EdgeSet = [(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7)]
            .iter()
            .cloned()
            .collect();
        let graph = create_test_graph(edges);
        let mut adm_graph = AdmGraph::new(&graph, 2);

        adm_graph.initialise_candidates();

        assert_eq!(
            adm_graph.candidates,
            [3, 4, 5, 6, 7].iter().cloned().collect()
        );
    }

    #[test]
    fn is_all_vertices_in_r_or_candidates_returns_true_if_all_vertices_are_in_r_or_candidates() {
        let graph = EditGraph::new();
        let mut adm_graph = AdmGraph::new(&graph, 1);
        adm_graph.r = vec![1,2,3].into_iter().collect();
        adm_graph.candidates = vec![4,5].into_iter().collect();
        adm_graph.num_of_vertices = 5;

        assert_eq!(adm_graph.is_all_vertices_in_r_or_candidates(), true);
    }

    #[test]
    fn is_all_vertices_in_r_or_candidates_returns_false_if_some_vertices_are_not_in_r_or_candidates() {
        let graph = EditGraph::new();
        let mut adm_graph = AdmGraph::new(&graph, 1);
        adm_graph.r = vec![1,2,3].into_iter().collect();
        adm_graph.candidates = vec![4,5].into_iter().collect();
        adm_graph.num_of_vertices = 10;

        assert_eq!(adm_graph.is_all_vertices_in_r_or_candidates(), false);
    }

    #[test]
    fn get_next_v_returns_v_in_candidates() {
        let edges: EdgeSet = [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4)]
            .iter()
            .cloned()
            .collect();
        let graph = EditGraph::new();
        let graph = create_test_graph(edges);
        let mut adm_graph = AdmGraph::new(&graph, 2);
        adm_graph.initialise_candidates();

        let next = adm_graph.get_next_v_in_ordering().unwrap();

        assert_eq!(adm_graph.r.contains(&next), true);
        assert_eq!(adm_graph.candidates.contains(&next), false);
        assert_eq!(adm_graph.l.contains(&next), false);
    }

    #[test]
    fn get_next_v_returns_none_if_no_candidates() {
        let edges: EdgeSet = [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4)]
            .iter()
            .cloned()
            .collect(); //Each vertex has degree > 1
        let graph = EditGraph::new();
        let graph = create_test_graph(edges);
        let mut adm_graph = AdmGraph::new(&graph, 1);
        adm_graph.initialise_candidates();

        let next = adm_graph.get_next_v_in_ordering();

        assert_eq!(next, None);
    }
}