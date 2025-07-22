use crate::adm_data::{self, AdmData, Path};
use crate::flow_network::FlowNetwork;
use crate::vias::Vias;
use graphbench::editgraph::EditGraph;
use graphbench::graph::{Graph, MutableGraph, Vertex, VertexMap, VertexSet};
use graphbench::iterators::EdgeIterable;

pub(crate) struct AdmGraph<'a> {
    pub l: VertexSet,
    pub r: VertexSet,
    pub candidates: VertexSet,
    pub adm_data: VertexMap<AdmData>,
    num_of_vertices: usize,
    p: usize,
    pub graph: &'a EditGraph,
    pub vias: Vias,
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
    fn compute_vias(&mut self, v_data: &mut AdmData, v_targets:&VertexSet) {
        v_data.delete_packing(); //clear 3-packing of v as we now want to store a 2-packing for v
        //TODO: Can we remove all vertices from v's T1 that are not in L?
        v_data.n_in_r.extend(v_data.t1.iter().filter(|x|self.r.contains(x)).cloned());
        v_data.t1.retain(|x| self.l.contains(x));

        let v = v_data.id;

        // Update via data structure
        let v_right_neighbours = v_data.n_in_r.clone();
        for &x in v_right_neighbours.iter() { // Right neighbours of v
            for &y in v_data.t1.iter() { // Left neighbours of v
                self.vias.add_a_via(x, y, v);
            }
            let x_data = self.adm_data.get(&x).unwrap();
            for &y in x_data.t1.iter() {
                self.vias.add_a_via(v, y, x);
            }
        }

        // Update other 2-packings
        println!("Updating 2-packing for neighbours {:?}", v_right_neighbours);
        for &u in v_right_neighbours.iter() {
            let u_data = self.adm_data.get_mut(&u).unwrap();
            debug_assert!(u_data.t1.contains(&v), "T1 for u={:?} does not contain {:?}", u, v); 

            u_data.remove_v_from_packing(&v);

            for y in v_data.t1.iter() {
                if u_data.is_v_in_pack(y) {
                    continue
                }
                debug_assert!(self.l.contains(y));
                debug_assert!(v_data.t1.contains(y));
                debug_assert!(!u_data.t1.contains(y));
                debug_assert!(self.graph.adjacent(&u, &v));
                debug_assert!(!self.graph.adjacent(&u, y)); // Cordless
                u_data.add_t2_to_packing(&v, y);
                break                
            }
        }

        // Compute 2-packing for v
        println!("Finding maximal 2-packing for targets {:?}", v_targets);
        for &y in v_targets {
            if self.graph.adjacent(&y, &v) {
                debug_assert!(v_data.t1.contains(&y)); 
                continue
            }

            debug_assert!(self.r.contains(&v));
            debug_assert!(self.l.contains(&y));
            if let Some(vias) = self.vias.get_vias(v, y) {
                for x in vias.iter() {
                    if v_data.is_v_in_pack(&x) {
                        continue
                    }
                    debug_assert!(self.graph.adjacent(&v, &x));
                    debug_assert!(!self.graph.adjacent(&v, &y)); // Cordless
                    v_data.add_t2_to_packing(&x, &y);
                    break;
                }
            } 
        }
        println!("2-packing for v={}: T1 = {:?}, T2 = {:?}", v, v_data.t1, v_data.packing);
        println!("  actual left neighbours: {:?}", self.graph.neighbours(&v).filter(|x| self.l.contains(x)).collect::<Vec<&Vertex>>())
    }

    /*
        Fetches all the vertices in L that is in T1, T2 and T3 of v
    */
    fn collect_targets(&self, v_data: &AdmData) -> VertexSet {
        let v = v_data.id;
        let mut targets: VertexSet = v_data.t1.intersection(&self.l).cloned().collect();
        let t1_t2: VertexSet = v_data.t1.union(&v_data.t2).cloned().collect();

        // println!("Collect targets for {v}:");
        // println!("  T1 = {:?}", targets);
        // println!("  Packing = {:?}", v_data.packing);

        for w in t1_t2.intersection(&self.r) {
            if *w == v {
                continue; // We already added t1 to t
            }
            
            let w_data = self.adm_data.get(w).unwrap();

            // Add left neigbhours of w to targets
            targets.extend(w_data.t1.intersection(&self.l));
            // println!("  Adding left neighbours of {}: {:?}", w, w_data.t1.intersection(&self.l).cloned().collect::<VertexSet>() );

            // The remaining code is only for vertices at distance 1 from v
            if !self.graph.adjacent(w, &v) {
                continue;
            };

            // To find the targets at distance three, we use the 2-packing stored for
            // w.
            for x in w_data.t1.intersection(&self.r) {
                if *x == v {
                    continue;
                }
                let x_data = self.adm_data.get(x).unwrap();
                // println!("  Adding left neighbours of {}: {:?}", x, x_data.t1.intersection(&self.l).cloned().collect::<VertexSet>() );                
                targets.extend(x_data.t1.intersection(&self.l));
            }
        }
        targets.remove(&v);
        targets
    }

    /*
        Simple update of a packing
    */
    fn simple_update(&mut self, u_data: &mut AdmData, v: Vertex) {
        if !u_data.is_v_in_pack(&v) {
            return;
        }
        let u = u_data.id;

        println!("\nSimple update for u={u}");
        println!("  T1 = {:?} ", u_data.t1.intersection(&self.l).cloned().collect::<VertexSet>());
        println!("  Packing = {:?}", u_data.packing);

        let removed_path = u_data.remove_v_from_packing(&v);
        println!("Removed path {:?}", removed_path);
        u_data.debug_check_consistency(self); // DEBUG

        let (w, removed_len) = match removed_path {
            Some(path) => match path {
                Path::TwoPath(x, _) => (x, 2),
                Path::ThreePath(u, _, _) => (u, 3),
            },
            None => (v, 1),
        };
        let w_data = self.adm_data.get(&w).unwrap();
        debug_assert!(self.r.contains(&w));

        // Check if we can add a path of length 2
        println!("  Attempting to add path of length 2 through {:?}", w_data.t1.intersection(&self.l).cloned().collect::<VertexSet>());
        for x in w_data.t1.intersection(&self.l) {
            if *x == u {
                continue;
            }
            if !u_data.is_v_in_pack(x) {
                debug_assert!(self.l.contains(x));
                debug_assert!(!self.graph.adjacent(&u, x)); // Cordless
                u_data.add_t2_to_packing(&w, x);
                return;
            }
        }

        // Check if we can add a path of length 3
        println!("  Attempting to add path of length 3 through {:?}", w_data.t1.intersection(&self.r).cloned().collect::<VertexSet>());
        for x in w_data.t1.intersection(&self.r) {
            if u_data.is_v_in_pack(x) || (u == *x) {
                continue;
            }

            let x_adm_data = self.adm_data.get(x).unwrap();
            for y in x_adm_data.t1.intersection(&self.l) {
                if *y == u {
                    continue;
                }
                if !u_data.is_v_in_pack(y) {
                    //Check if there is a shorter path to y (i.e u,x,y)
                    // if so add the shorter path instead
                    if self.graph.adjacent(&u, x) {
                        debug_assert!(self.r.contains(x));
                        debug_assert!(self.l.contains(y));
                        debug_assert!(self.graph.adjacent(&u, &x));
                        debug_assert!(!self.graph.adjacent(&u, y)); // Cordless
                        u_data.add_t2_to_packing(x, y); 
                    } else {
                        debug_assert!(self.r.contains(&w));
                        debug_assert!(self.r.contains(x));
                        debug_assert!(self.l.contains(y));
                        debug_assert!(self.graph.adjacent(&u, &w));
                        debug_assert!(!self.graph.adjacent(&u, y)); // Cordless
                        debug_assert!(!self.graph.adjacent(&u, x)); // Cordless
                        u_data.add_t3_to_packing(&w, x, y);
                    }
                    return;
                }
            }
        }

        // Check if we can add a path of length 3 through v by
        // finding a via x in Vias[v][u] and a vertex y in N_L(v)
        if removed_len != 3 {
            return
        }

        let v_data = self.adm_data.get(&v).unwrap();
        debug_assert!(self.r.contains(&v));
        debug_assert!(self.l.contains(&u));
        let vias = if let Some(vias) = self.vias.get_vias(v, u) {
            vias
        } else {
            return // No vias, no path
        };

        for &y in v_data.t1.intersection(&self.l) {
            if u_data.is_v_in_pack(&y) || (u == y) {
                continue // Already part of Pack(u)
            }
            for &x in vias {
                if u_data.is_v_in_pack(&x) || (u == x) {
                    continue; // Already part of pack(u)
                }

                debug_assert!(self.graph.adjacent(&u, &x)); 
                debug_assert!(!self.graph.adjacent(&u, &y)); // Cordless
                debug_assert!(!self.graph.adjacent(&u, &v)); // Cordless
                u_data.add_t3_to_packing(&x, &v, &y);
            }
        }
    }

    /*
       Try to see if a disjoint path can be added to the packing of u
    */
    fn stage_1_update(&mut self, u_data: &mut AdmData, u_targets: &VertexSet) {
        let u = u_data.id;
        let t1: VertexSet = u_data.t1.intersection(&self.l).cloned().collect();
        let t_23: VertexSet = u_targets.difference(&t1).cloned().collect();

        println!("\nStage 1 update for u={u}");
        for w in u_data.n_in_r.clone() {
            debug_assert!(self.r.contains(&w));
            if u_data.is_v_in_pack(&w) {
                continue;
            }

            for y in &t_23 {
                debug_assert!(self.l.contains(y));
                if *y == u || u_data.is_v_in_pack(y){
                    continue;
                }

                let w_data = self.adm_data.get(&w).unwrap();

                // First check if uwy is a path
                if w_data.t1.contains(y) {
                    debug_assert!(self.l.contains(y));
                    debug_assert!(self.graph.adjacent(&u, &w));
                    debug_assert!(self.graph.adjacent(&w, &y));
                    debug_assert!(!self.graph.adjacent(&u, y)); // Cordless
                    u_data.add_t2_to_packing(&w, y);
                    return
                }

                // If uwy is not a path check vias between w and y
                // to see whether a path uwxy exists.
                debug_assert!(self.r.contains(&w));
                debug_assert!(self.l.contains(y));
                let x_vias = self.vias.get_vias(w, *y);
                if x_vias.is_none() {
                    continue;
                }

                for x in x_vias.unwrap() {
                    if u_data.is_v_in_pack(x) {
                        continue;
                    }

                    // First check if there is a shorter path uxy to reach y
                    if self.graph.adjacent(&u, x) {
                        debug_assert!(self.l.contains(y));
                        debug_assert!(self.graph.adjacent(&u, &x));
                        debug_assert!(self.graph.adjacent(&x, &y));
                        debug_assert!(!self.graph.adjacent(&u, y)); // Cordless
                        u_data.add_t2_to_packing(x, y);
                    } else {
                        debug_assert!(self.l.contains(y));
                        debug_assert!(self.graph.adjacent(&u, &w));
                        debug_assert!(self.graph.adjacent(&w, x));
                        debug_assert!(self.graph.adjacent(x, y));
                        debug_assert!(!self.graph.adjacent(&u, y)); // Cordless
                        debug_assert!(!self.graph.adjacent(&u, x)); // Cordless
                        u_data.add_t3_to_packing(&w, x, y);
                    }

                    // In either case, we extended the size of the packing by one
                    return
                }
            }
        }
    }

    /*
       Find an augmenting path to see if packing of u can be extended
    */
    fn stage_2_update(&mut self, u_data: &mut AdmData, targets: &VertexSet) {
        println!("\nStage 2 update for u={}", u_data.id);
        let mut flow_network = FlowNetwork::new(u_data.id);
        flow_network.construct_flow_network(self, u_data, targets);
        flow_network.augmenting_path(u_data, &self);
    }

    /*
        Update data structures now that vertx v is moving to R
    */
    fn update(&mut self, v: &Vertex) {
        let mut v_adm_data = self.adm_data.remove(v).unwrap();
        debug_assert!(v_adm_data.size_of_packing() <= self.p);
        let t = self.collect_targets(&v_adm_data);
        
        self.compute_vias(&mut v_adm_data, &t);
        self.adm_data.insert(*v, v_adm_data);

        debug_assert!(!self.l.contains(v));
        debug_assert!(self.r.contains(v));

        println!("Updating 3-packing for targets {:?}", t);
        for u in t {
            debug_assert!(self.l.contains(&u));

            if self.candidates.contains(&u) {
                // TEMP: Update packings in candidates as well.
                let mut u_adm_data = self.adm_data.remove(&u).unwrap();
                self.simple_update(&mut u_adm_data, *v);   
                self.adm_data.insert(u, u_adm_data);             
                continue;
            }

            let mut u_adm_data = self.adm_data.remove(&u).unwrap();
            self.simple_update(&mut u_adm_data, *v);
            u_adm_data.debug_check_consistency(self);


            //check if a disjoint path can be added to packing of u
            if u_adm_data.size_of_packing() <= self.p {
                let targets_u = self.collect_targets(&u_adm_data);
                self.stage_1_update(&mut u_adm_data, &targets_u);

                //If size of packing is still p after simple update
                //then do augmenting path to see if packing of u can be extended
                if u_adm_data.size_of_packing() <= self.p {
                    self.stage_2_update(&mut u_adm_data, &targets_u);
                }
            }

            //If size of packing is still p then add to candidates
            if u_adm_data.size_of_packing() <= self.p {
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

    pub fn debug_check_consistency(&self) {
        if !cfg!(debug_assertions) {
            return
        }

        let mut abort = false;
        let mut log = "".to_string();
        println!("\nDEBUG: Checking consistency");

        // Check adm vertex data for consistency
        for v_data in self.adm_data.values() {
            let v = v_data.id;
            if self.l.contains(&v) {
                let targets = self.collect_targets(v_data);

                let left_1:VertexSet = self.graph.neighbours(&v).filter(|x| self.l.contains(x)).cloned().collect();
                let right_1:VertexSet = self.graph.neighbours(&v).filter(|x| self.r.contains(x)).cloned().collect();
                let right_2:VertexSet = self.graph.neighbourhood(right_1.iter()).intersection(&self.r).cloned().collect();
                let right_all:VertexSet = right_1.union(&right_2).cloned().collect();
                
                let mut targets_check:VertexSet = self.graph.neighbourhood(right_all.iter()).into_iter().filter(|x| *x != v && self.l.contains(x)).collect();
                targets_check.extend(left_1);

                // Sorted list for nicer output
                let mut targets:Vec<Vertex> = targets.into_iter().collect();
                targets.sort();
                let mut targets_check:Vec<Vertex> = targets_check.into_iter().collect();
                targets_check.sort();

                if targets != targets_check {
                    log += &format!("\nTargets3 for vertex {v} (L) do not match:");
                    log += &format!("\n  collect targets = {:?}", targets);
                    log += &format!("\n     real targets = {:?}", targets_check);
                    log += &format!("\n 3-Packing = {:?} ", v_data.packing);
                    abort = true;
                }
            } else {
                // TODO: check 2-packing
                let left_1:VertexSet = self.graph.neighbours(&v).filter(|x| self.l.contains(x)).cloned().collect();
                let right_1:VertexSet = self.graph.neighbours(&v).filter(|x| self.r.contains(x)).cloned().collect();

                // We have to do this by hand, collect_targets always collects up to distance 3
                let mut targets:VertexSet = VertexSet::default();
                targets.extend(v_data.t1.intersection(&self.l));
                for (_, path) in v_data.packing.iter() {
                    match path {
                        Path::TwoPath(x, _) => {
                            let x_data = self.adm_data.get(x).unwrap();
                            targets.extend(x_data.t1.intersection(&self.l));
                        },
                        Path::ThreePath(_, _, _) => unreachable!(),
                    }
                }

                let mut targets_check:VertexSet = self.graph.neighbourhood(right_1.iter()).into_iter().filter(|x| self.l.contains(x)).collect();
                targets_check.extend(left_1);

                // Sorted list for nicer output
                let mut targets:Vec<Vertex> = targets.into_iter().collect();
                targets.sort();
                let mut targets_check:Vec<Vertex> = targets_check.into_iter().collect();
                targets_check.sort();

                if targets != targets_check {
                    log += &format!("\nTargets2 for vertex {v} (R) do not match:");
                    log += &format!("\n  2-pack targets = {:?}", targets);
                    log += &format!("\n    real targets = {:?}", targets_check);
                    log += &format!("\n 2-Packing = {:?} ", v_data.packing);
                    abort = true;
                }
            };
        }

        if abort {
            println!("Issues found:");
            println!("{log}");
            panic!();
        }

        // Also check basics
        for v_data in self.adm_data.values() {
            v_data.debug_check_consistency(&self);
        }        
    }

    /*
       Gets next vertex in the ordering
    */
    pub fn get_next_v_in_ordering(&mut self) -> Option<Vertex> {
        let v = self.candidates.iter().next();

        println!();
        println!("=== Moving {:?} from L to R ===" , v);
        println!();

        match v {
            Some(&v) => {
                debug_assert!(self.l.contains(&v));
                debug_assert!(!self.r.contains(&v));
                self.candidates.remove(&v);
                self.l.remove(&v);
                self.r.insert(v);

                println!("  L = {:?}", self.l);
                println!("  R = {:?}", self.r);

                self.update(&v);
                Some(v)
            }
            None => None,
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
        adm_graph.r = vec![1, 2, 3].into_iter().collect();
        adm_graph.candidates = vec![4, 5].into_iter().collect();
        adm_graph.num_of_vertices = 5;

        assert!(adm_graph.is_all_vertices_in_r_or_candidates());
    }

    #[test]
    fn is_all_vertices_in_r_or_candidates_returns_false_if_some_vertices_are_not_in_r_or_candidates()
     {
        let graph = EditGraph::new();
        let mut adm_graph = AdmGraph::new(&graph, 1);
        adm_graph.r = vec![1, 2, 3].into_iter().collect();
        adm_graph.candidates = vec![4, 5].into_iter().collect();
        adm_graph.num_of_vertices = 10;

        assert!(!adm_graph.is_all_vertices_in_r_or_candidates());
    }

    #[test]
    fn get_next_v_returns_v_in_candidates() {
        let edges: EdgeSet = [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4)]
            .iter()
            .cloned()
            .collect();
        let graph = create_test_graph(edges);
        let mut adm_graph = AdmGraph::new(&graph, 2);
        adm_graph.initialise_candidates();

        let next = adm_graph.get_next_v_in_ordering().unwrap();

        assert!(adm_graph.r.contains(&next));
        assert!(!adm_graph.candidates.contains(&next));
        assert!(!adm_graph.l.contains(&next));
    }

    #[test]
    fn get_next_v_returns_none_if_no_candidates() {
        let edges: EdgeSet = [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4)]
            .iter()
            .cloned()
            .collect(); //Each vertex has degree > 1
        let graph = create_test_graph(edges);
        let mut adm_graph = AdmGraph::new(&graph, 1);
        adm_graph.initialise_candidates();

        let next = adm_graph.get_next_v_in_ordering();

        assert_eq!(next, None);
    }
}
