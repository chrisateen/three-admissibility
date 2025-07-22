use crate::adm_data::{AdmData, Path};
use crate::adm_graph::AdmGraph;
use crate::vias::Vias;
use graphbench::editgraph::EditGraph;
use graphbench::graph::{EdgeSet, Graph, MutableGraph, Vertex, VertexMap, VertexSet};
use graphbench::algorithms::*;

use std::collections::{HashMap, HashSet, VecDeque};
use std::default;
use std::ops::Neg;

pub struct FlowNetwork {
    pub id: Vertex,
    pub arcs: VertexMap<VertexSet>,
    pub s1: VertexSet,    // t1_r in packing
    pub s2: VertexSet,    // t2_r in packing
    pub t_in: VertexSet,  // targets in packing
    pub t_out: VertexSet, //target not in packing
}

impl FlowNetwork {
    pub fn new(id: Vertex) -> Self {
        FlowNetwork {
            id,
            arcs: VertexMap::default(),
            s1: VertexSet::default(),
            s2: VertexSet::default(),
            t_in: VertexSet::default(),
            t_out: VertexSet::default(),
        }
    }

    /*
        Add edges from pack of self.id
        edges added in a direction that points away from root
    */
    fn add_pack_edges(&mut self, packing: &VertexMap<Path>) {
        let mut root_neighbours = self.arcs.remove(&self.id).unwrap_or_default();

        for w in packing.keys() {
            self.t_in.insert(*w);
            let path = packing.get(w).unwrap();

            match path {
                Path::TwoPath(s1, t2) => {
                    let s1_neighbours = self.arcs.entry(*s1).or_default();
                    root_neighbours.insert(*s1);
                    self.s1.insert(*s1);
                    s1_neighbours.insert(*t2);
                    debug_assert_ne!(*t2, self.id);
                }
                Path::ThreePath(s1, s2, t3) => {
                    let mut s1_neighbours = self.arcs.remove(s1).unwrap_or_default();
                    let mut s2_neighbours = self.arcs.remove(s2).unwrap_or_default();

                    root_neighbours.insert(*s1);
                    self.s1.insert(*s1);
                    s1_neighbours.insert(*s2);
                    self.s2.insert(*s2);
                    s2_neighbours.insert(*t3);

                    debug_assert!(!s1_neighbours.contains(&self.id));
                    debug_assert!(!s2_neighbours.contains(&self.id));

                    self.arcs.insert(*s1, s1_neighbours);
                    self.arcs.insert(*s2, s2_neighbours);
                }
            }
        }

        debug_assert!(!root_neighbours.contains(&self.id));
        self.arcs.insert(self.id, root_neighbours);
    }

    /*
        Add all edges between vertices in packing
        edges added in a direction that points away from root
    */
    fn add_edges_between(&mut self, graph: &EditGraph) {
        let s2_union_t: VertexSet = self.s2.union(&self.t_in).cloned().collect();
        for v in &self.s1 {
            for w in &s2_union_t {
                if graph.adjacent(v, w) {
                    //edges are added as s1 -> s2 or s1 -> t
                    let v_neighbours = self.arcs.entry(*v).or_default();
                    v_neighbours.insert(*w);
                    debug_assert_ne!(*w, self.id);
                }
            }
        }

        for v in &self.s2 {
            for w in &self.t_in {
                if graph.adjacent(v, w) {
                    //edges are added as s2 -> t
                    let v_neighbours = self.arcs.entry(*v).or_default();
                    v_neighbours.insert(*w);
                    debug_assert_ne!(*w, self.id);
                }
            }
        }
    }

    /*
        Add all edges between N_R(v) that are not in the packing
         and vertices in pack that are in s2 and t
    */
    fn add_direct_edge_from_n_r(&mut self, n_in_r: &VertexSet, graph: &EditGraph) {
        let s2_union_t: VertexSet = self.s2.union(&self.t_in).cloned().collect();
        let mut root_neighbours = self.arcs.remove(&self.id).unwrap_or_default();

        for y in &s2_union_t {
            for x in n_in_r.difference(&self.s1) {
                if graph.adjacent(x, y) {
                    // if self.id == 30 && *y == 70 {
                    //     panic!();
                    // }
                    // if self.id == 30 && *x == 70 {
                    //     panic!(); <--- this one 
                    // }                    

                    // Add edge from root -> N_R(v) / x
                    root_neighbours.insert(*x);

                    // Edges are added as N_R(v)/ x -> s2 or N_R(v) -> t
                    let x_neighbours = self.arcs.entry(*x).or_default();
                    x_neighbours.insert(*y);
                    debug_assert_ne!(*y, self.id);
                    break;
                }
            }
        }

        self.arcs.insert(self.id, root_neighbours);
    }

    /*
        Add all edges N_R(v) and t through a via
    */
    fn add_edge_with_via_from_n_r(&mut self, n_in_r: &VertexSet, vias: &Vias) {
        let mut root_neighbours = self.arcs.remove(&self.id).unwrap_or_default();
        let s: VertexSet = self.s1.union(&self.s2).cloned().collect();

        'outer: for y in &self.t_in {
            for x in n_in_r.difference(&self.s1) {
                let vias_x_y = vias.get_vias(*x, *y);
                if vias_x_y.is_some() {
                    let eligible_via = vias_x_y.unwrap().difference(&s).next();
                    if let Some(w) = eligible_via {
                        //add edge from root -> N_R(v) / x
                        root_neighbours.insert(*x);
                        debug_assert_ne!(*x, self.id);

                        //add edge from N_R(v) -> via / w
                        let x_neighbours = self.arcs.entry(*x).or_default();
                        x_neighbours.insert(*w);
                        debug_assert_ne!(*w, self.id);

                        //add edge from via -> t / y
                        let w_neighbours = self.arcs.entry(*w).or_default();
                        w_neighbours.insert(*y);
                        debug_assert_ne!(*y, self.id);

                        continue 'outer;
                    }
                }
            }
        }
        self.arcs.insert(self.id, root_neighbours);
    }

    /*
       Add a via between each s1 and target
    */
    fn add_vias_from_s1(&mut self, vias: &Vias) {
        let s: VertexSet = self.s1.union(&self.s2).cloned().collect();

        'outer: for v in &self.s1 {
            for w in &self.t_in {
                let vias_v_w = vias.get_vias(*v, *w);
                if vias_v_w.is_some() {
                    let eligible_via = vias_v_w.unwrap().difference(&s).next();
                    if let Some(x) = eligible_via {
                        //edges are added as s1 -> via
                        let v_neighbours = self.arcs.entry(*v).or_default();
                        v_neighbours.insert(*x);
                        debug_assert_ne!(*x, self.id);

                        //then via -> t
                        let x_neighbours = self.arcs.entry(*x).or_default();
                        x_neighbours.insert(*w);
                        debug_assert_ne!(*w, self.id);

                        continue 'outer;
                    }
                }
            }
        }
    }

    /*
        For each s1 and s2 add a target that is not already in the packing
    */
    fn add_extra_targets(
        &mut self,
        adm_data: &VertexMap<AdmData>,
        u_data: &AdmData,        
        targets: &VertexSet,
        vias: &Vias,
        l: &VertexSet,
    ) {
        let s1_s2: VertexSet = self.s1.union(&self.s2).cloned().collect();
        let u = u_data.id;

        'outer: for w in &s1_s2 {
            let w_data = adm_data.get(w).unwrap();
            for y in w_data.t1.intersection(l) { // Left neighbours w of v 
                if *y == u ||  u_data.t1.contains(y)  {
                    continue // w is equal to u or is a neighbour of root u
                }

                if !self.t_in.contains(y) {
                    //adds edge s2 -> target outside packing or s1 -> target outside packing
                    self.arcs.get_mut(w).unwrap().insert(*y);
                    debug_assert_ne!(*y, u);
                    self.t_out.insert(*y);
                    continue 'outer;
                }
            }

            if self.s2.contains(w) {
                continue 'outer;
            }

            for w in targets.difference(&self.t_in) {
                let v_w_vias = vias.get_vias(*w, *w);
                if v_w_vias.is_none() {
                    continue;
                }
                let eligible_via = v_w_vias.unwrap().difference(&s1_s2).next();
                match eligible_via {
                    None => {}
                    Some(x) => {
                        //adds edge s1 (v) -> via (x)
                        self.arcs.get_mut(w).unwrap().insert(*x);
                        self.arcs.entry(*x).or_default().insert(*w);
                        debug_assert_ne!(*w, u_data.id);
                        self.t_out.insert(*w);
                        continue 'outer;
                    }
                }
            }
        }
    }

    /*
        Spilt edges in network to prepare for augmenting path
        create a duplicate vertex which a negation of the original vertex
    */
    fn split_edges_in_network(&self) -> HashMap<i32, HashSet<i32>> {
        let mut edges = HashMap::default();

        for (v, neighbours) in self.arcs.iter() {
            let negative_set = neighbours.iter().map(|&x| -(x as i32)).collect();
            let negative_v = -(*v as i32) as u32;
            edges.insert(*v as i32, negative_set);

            if *v != self.id {
                edges.insert(
                    negative_v as i32,
                    [*v as i32].iter().copied().collect::<HashSet<i32>>(),
                );
            }
        }

        for v in self.t_in.union(&self.t_out) {
            let negative_v: Vertex = -(*v as i32) as u32;
            edges.insert(
                negative_v as i32,
                [*v as i32].iter().copied().collect::<HashSet<i32>>(),
            );
        }

        edges
    }

    /*
       Reverse direction of edges in original packing
    */
    fn set_edges_direction(
        &self,
        edges: HashMap<i32, HashSet<i32>>,
        u: &AdmData,
    ) -> HashMap<i32, HashSet<i32>> {
        let mut flow = edges.clone();
        let root = self.id as i32;
        let mut root_neighbours = flow.remove(&root).unwrap();

        for w in u.packing.keys() {
            let path = u.packing.get(w).unwrap();

            match path {
                Path::TwoPath(s1, t2) => {
                    let (s1, t2) = (*s1 as i32, *t2 as i32);
                    let (neg_s1, neg_t2) = (-s1, -t2);

                    //reverse root -> -s1 -> s1 -> -t2 -> t2
                    // to t2 -> -t2 -> s1 -> -s1 -> root
                    root_neighbours.remove(&neg_s1);

                    let negative_s1_neighbours = flow.get_mut(&neg_s1).unwrap();
                    negative_s1_neighbours.insert(root);
                    negative_s1_neighbours.remove(&s1);

                    let s1_neighbours = flow.get_mut(&s1).unwrap();
                    s1_neighbours.insert(neg_s1);
                    s1_neighbours.remove(&neg_t2);

                    let negative_t2_neighbours = flow.get_mut(&neg_t2).unwrap();
                    negative_t2_neighbours.insert(s1);
                    negative_t2_neighbours.remove(&t2);

                    let t2_neighbours = flow.entry(t2).or_default();
                    t2_neighbours.insert(neg_t2);
                }
                Path::ThreePath(s1, s2, t3) => {
                    let (s1, s2, t3) = (*s1 as i32, *s2 as i32, *t3 as i32);
                    let (neg_s1, neg_s2, neg_t3) = (-s1, -s2, -t3);

                    //reverse root -> -s1 -> s1 -> -s2 -> s2 -> -t3 -> t3
                    // to t3 -> -t3 -> s2 -> -s2 -> s1 -> -s1 -> root
                    root_neighbours.remove(&neg_s1);

                    let negative_s1_neighbours = flow.get_mut(&neg_s1).unwrap();
                    negative_s1_neighbours.insert(root);
                    negative_s1_neighbours.remove(&s1);

                    let s1_neighbours = flow.get_mut(&s1).unwrap();
                    s1_neighbours.insert(neg_s1);
                    s1_neighbours.remove(&neg_s2);

                    let negative_s2_neighbours = flow.get_mut(&neg_s2).unwrap();
                    negative_s2_neighbours.insert(s1);
                    negative_s2_neighbours.remove(&s2);

                    let s2_neighbours = flow.get_mut(&s2).unwrap();
                    s2_neighbours.insert(neg_s2);
                    s2_neighbours.remove(&neg_t3);

                    let negative_t3_neighbours = flow.get_mut(&neg_t3).unwrap();
                    negative_t3_neighbours.insert(s2);
                    negative_t3_neighbours.remove(&t3);

                    let t3_neighbours = flow.entry(t3).or_default();
                    t3_neighbours.insert(neg_t3);
                }
            }
        }
        flow.insert(root, root_neighbours);
        flow
    }

    /*
        Used to find the shortest path from v (a neighbour of the root vertex) to vertex in t_out
    */
    fn bfs(&self, flow: &HashMap<i32, HashSet<i32>>) -> Option<Vec<i32>> {
        //Keeps a track of the path from endpoint to root
        //Hashmap stores child -> parent to make it easy to get the right path from sink
        let mut path: HashMap<i32, i32> = HashMap::default();
        let mut visited: HashSet<i32> = HashSet::default();
        let mut queue: VecDeque<i32> = VecDeque::new();
        let root = self.id as i32;

        visited.insert(root);
        queue.push_back(root);

        while let Some(u) = queue.pop_front() {
            let u_neighbours = if let Some(u_neighbours) = flow.get(&u) {
                u_neighbours
            } else {
                continue
            };

            for w in u_neighbours {
                if !visited.contains(w) {
                    queue.push_back(*w);
                    visited.insert(*w);
                    path.insert(*w, u);
                }
                if self.t_out.contains(&(*w as Vertex)) {
                    let mut result_path = vec![*w];
                    let mut current = w;
                    while *current != root {
                        let next = path.get(current).unwrap();
                        result_path.push(*next);
                        current = next;
                    }
                    result_path.reverse();
                    return Some(result_path);
                }
            }
        }
        None
    }


    /*
        Update packing based on results of augmenting path
    */
    fn update_packing(&self, u: &mut AdmData, aug_path:Vec<Vertex>, adm_graph:&AdmGraph) {
        // Collect edges from augmenting path. To ensure consitency, we store edge in sorted vertex-order
        let aug_path_edges:EdgeSet = aug_path.windows(2).map(|e| if e[0] < e[1] {(e[0],e[1])} else {(e[1],e[0])}).collect();

        println!("Aug path = {:?}", aug_path);
        println!("Aug path edges = {:?}", aug_path_edges);

        // Collect packing edges
        let root = u.id;        
        println!("Root = {root}");

        println!("{:?}", u.packing);

        let mut packing_edges: EdgeSet = EdgeSet::default();
        for path in u.packing.values() {
            match path {
                Path::TwoPath(x, y) => { // root -> x -> y
                    debug_assert!(adm_graph.l.contains(y));
                    packing_edges.insert( if root < *x { (root,*x) } else { (*x,root)} );
                    packing_edges.insert( if x < y { (*x,*y) } else { (*y,*x)} );

                }
                Path::ThreePath(w, x, y) => { // root -> w -> x -> y
                    debug_assert!(adm_graph.l.contains(y));
                    packing_edges.insert( if root < *w { (root,*w) } else { (*w,root)} );
                    packing_edges.insert( if w < x { (*w,*x) } else { (*x,*w)} );
                    packing_edges.insert(  if x < y { (*x,*y) } else { (*y,*x)} );
                }
            };
        }     

        println!("Packing edges = {:?}", packing_edges);
        debug_assert!(aug_path_edges.iter().all(|(x,y)| x < y));
        debug_assert!(packing_edges.iter().all(|(x,y)| x < y));

        // Construct aux. graph        
        let mut aux = EditGraph::new();
        for (x,y) in aug_path_edges.symmetric_difference(&packing_edges) {
            aux.add_edge(x, y);
        }

        debug_assert!(aux.contains(&root));
        println!("Aux graph = {:?}", aux);        

        // Collect start of paths (neighbours of root)
        u.delete_packing(); // Delete old packing
        u.t1.retain(|x| adm_graph.l.contains(x));


        let root_neighbours:VertexSet = aux.neighbours(&root).cloned().collect();
        aux.remove_vertex(&root);
        for x in root_neighbours {
            debug_assert_eq!(aux.degree(&x), 1); // We cannot have a path of length 1 here!
            let y = *aux.neighbours(&x).next().unwrap();
            aux.remove_vertex(&x);

            if aux.degree(&y) == 0 {
                // Found a path of length 2: x y
                // debug_assert!()
                println!("Adding path {}-{}-{}", u.id, x, y);
                debug_assert!(adm_graph.graph.adjacent(&u.id, &x));
                debug_assert!(!adm_graph.graph.adjacent(&u.id, &y)); // Cordless
                u.add_t2_to_packing(&x, &y);
            } else if aux.degree(&y) == 1 {
                // Path of length 3: x y z
                let z = aux.neighbours(&y).next().unwrap();
                println!("Adding path {}-{}-{}-{}", u.id, x, y, z);
                debug_assert!(adm_graph.graph.adjacent(&u.id, &x));
                debug_assert!(!adm_graph.graph.adjacent(&u.id, &y)); // Cordless
                debug_assert!(!adm_graph.graph.adjacent(&u.id, &z)); // Cordless
                u.add_t3_to_packing(&x, &y, z);
            } else {
                unreachable!("This should not happen");
            }
        }
    }

    pub fn construct_flow_network(
        &mut self,
        adm_graph: &AdmGraph,
        u: &AdmData,
        targets: &VertexSet,
    ) {
        // DEBUG
        u.debug_check_consistency(adm_graph);

        self.add_pack_edges(&u.packing);
        self.add_edges_between(adm_graph.graph);
        self.add_direct_edge_from_n_r(&u.n_in_r, adm_graph.graph);

        // polbooks crash: In the previous method, the vertex 70
        // is added to the flown network for u = 30


        self.add_edge_with_via_from_n_r(&u.n_in_r, &adm_graph.vias);
        self.add_vias_from_s1(&adm_graph.vias);
        self.add_extra_targets(&adm_graph.adm_data, u, targets, &adm_graph.vias, &adm_graph.l);
    }

    pub fn augmenting_path(&self, u: &mut AdmData, adm_graph:&AdmGraph) {
        println!("Root = {}", u.id);
        println!("Arcs = {:?}", self.arcs);
        println!("T in = {:?}", self.t_in);
        println!("T out = {:?}", self.t_out);

        let split_edges = self.split_edges_in_network();
        let flow = self.set_edges_direction(split_edges, u);

        match self.bfs(&flow) {
            None => {}
            Some(path) => {
                println!("{:?}", path);
                let mut path_without_duplicates = Vec::with_capacity(path.len());
                let mut last_vertex = None;
                for v in path {
                    let v = if v < 0 { -v } else { v } as u32;
                    if last_vertex == Some(v) {
                        continue;
                    }
                    last_vertex = Some(v);
                    path_without_duplicates.push(v);
                }

                assert!(adm_graph.l.contains(path_without_duplicates.last().unwrap()));
                self.update_packing(u, path_without_duplicates, adm_graph);
                return;
            }
        }
    }
}

#[cfg(test)]
mod test_flow_network {
    use super::*;
    use graphbench::graph::{EdgeSet, MutableGraph};
    use crate::adm_data::Path;

    fn create_test_graph(edges: EdgeSet) -> EditGraph {
        let mut graph = EditGraph::new();
        for (u, v) in edges.iter() {
            graph.add_edge(u, v);
        }

        graph
    }

    #[test]
    fn add_pack_edges_add_edges_in_correct_order() {
        let mut packing: VertexMap<Path> = VertexMap::default();
        packing.insert(4, Path::ThreePath(2, 3, 4));
        packing.insert(6, Path::TwoPath(5, 6));
        let mut network = FlowNetwork::new(1);

        network.add_pack_edges(&packing);

        assert_eq!(
            *network.arcs.get(&1).unwrap(),
            [2, 5].into_iter().collect()
        );
        assert_eq!(*network.arcs.get(&2).unwrap(), [3].into_iter().collect());
        assert_eq!(*network.arcs.get(&3).unwrap(), [4].into_iter().collect());
        assert_eq!(*network.arcs.get(&5).unwrap(), [6].into_iter().collect());
        assert!(!network.arcs.contains_key(&4));
        assert!(!network.arcs.contains_key(&6));
        assert_eq!(network.s1, [2, 5].into_iter().collect());
        assert_eq!(network.s2, [3].into_iter().collect());
        assert_eq!(network.t_in, [4, 6].into_iter().collect());
    }

    #[test]
    fn add_edges_between_adds_edges() {
        let edges: EdgeSet = [(1, 2), (1, 5), (2, 3), (2, 6), (3, 4), (5, 3), (5, 6)]
            .iter()
            .cloned()
            .collect();
        let graph = create_test_graph(edges);
        let mut network = FlowNetwork {
            id: 1,
            s1: [2, 5].into_iter().collect(),
            s2: [3].into_iter().collect(),
            t_in: [4, 6].into_iter().collect(),
            t_out: VertexSet::default(),
            arcs: VertexMap::default(),
        };
        network
            .arcs
            .insert(1, [2, 5].into_iter().collect());
        network.arcs.insert(2, [3].into_iter().collect());
        network.arcs.insert(3, [4].into_iter().collect());
        network.arcs.insert(5, [6].into_iter().collect());

        network.add_edges_between(&graph);

        assert_eq!(
            *network.arcs.get(&2).unwrap(),
            [3, 6].into_iter().collect()
        );
        assert_eq!(
            *network.arcs.get(&5).unwrap(),
            [3, 6].into_iter().collect()
        );
    }

    #[test]
    fn add_direct_edge_from_n_r_adds_edge_between_n_r_and_s2() {
        let edges: EdgeSet = [(1, 2), (1, 5), (2, 3), (3, 4), (5, 3), (6, 3)]
            .iter()
            .cloned()
            .collect();
        let graph = create_test_graph(edges);
        let mut network = FlowNetwork {
            id: 1,
            s1: [2].into_iter().collect(),
            s2: [3].into_iter().collect(),
            t_in: [4].into_iter().collect(),
            t_out: VertexSet::default(),
            arcs: VertexMap::default(),
        };
        network.arcs.insert(1, [2].into_iter().collect());
        network.arcs.insert(2, [3].into_iter().collect());
        network.arcs.insert(3, [4].into_iter().collect());
        let n_in_r: VertexSet = [2, 5, 6].into_iter().collect();

        network.add_direct_edge_from_n_r(&n_in_r, &graph);

        assert_eq!(*network.arcs.get(&5).unwrap(), [3].into_iter().collect());
        assert_eq!(
            *network.arcs.get(&1).unwrap(),
            [2, 5].into_iter().collect()
        );
        //ensures only 1 edge is added from n_r to s2
        assert!(!network.arcs.contains_key(&6));
    }

    #[test]
    fn add_direct_edge_from_n_r_adds_edge_between_n_r_and_target() {
        let edges: EdgeSet = [(1, 2), (1, 5), (2, 3), (3, 4), (4, 5), (4, 6)]
            .iter()
            .cloned()
            .collect();
        let graph = create_test_graph(edges);
        let mut network = FlowNetwork {
            id: 1,
            s1: [2].into_iter().collect(),
            s2: [3].into_iter().collect(),
            t_in: [4].into_iter().collect(),
            t_out: VertexSet::default(),
            arcs: VertexMap::default(),
        };
        network.arcs.insert(1, [2].into_iter().collect());
        network.arcs.insert(2, [3].into_iter().collect());
        network.arcs.insert(3, [4].into_iter().collect());
        let n_in_r: VertexSet = [2, 5].into_iter().collect();

        network.add_direct_edge_from_n_r(&n_in_r, &graph);

        assert_eq!(*network.arcs.get(&5).unwrap(), [4].into_iter().collect());
        assert_eq!(
            *network.arcs.get(&1).unwrap(),
            [2, 5].into_iter().collect()
        );
        //ensures only 1 edge is added from n_r to t
        assert!(!network.arcs.contains_key(&6));
    }

    #[test]
    fn add_edge_with_via_from_n_r_adds_edge_with_via_between_n_r_and_target() {
        let mut network = FlowNetwork {
            id: 1,
            s1: [2].into_iter().collect(),
            s2: [3].into_iter().collect(),
            t_in: [4].into_iter().collect(),
            t_out: VertexSet::default(),
            arcs: VertexMap::default(),
        };
        network.arcs.insert(1, [2].into_iter().collect());
        network.arcs.insert(2, [3].into_iter().collect());
        network.arcs.insert(3, [4].into_iter().collect());
        let n_in_r: VertexSet = [2, 5, 7].into_iter().collect();
        let mut vias = Vias::new(10);
        let mut four_vias_1: VertexMap<VertexSet> = VertexMap::default();
        four_vias_1.insert(4, [6].into_iter().collect());
        let mut four_vias_2: VertexMap<VertexSet> = VertexMap::default();
        four_vias_2.insert(4, [8].into_iter().collect());
        vias.vias.insert(5, four_vias_1);
        vias.vias.insert(7, four_vias_2);

        network.add_edge_with_via_from_n_r(&n_in_r, &vias);

        assert_eq!(*network.arcs.get(&5).unwrap(), [6].into_iter().collect());
        assert_eq!(*network.arcs.get(&6).unwrap(), [4].into_iter().collect());
        assert_eq!(
            *network.arcs.get(&1).unwrap(),
            [2, 5].into_iter().collect()
        );
        //ensures only one of such edge is added for a target
        assert!(!network.arcs.contains_key(&7));
    }

    #[test]
    fn add_vias_from_s1_adds_edge_with_via_between_s1_and_target() {
        let mut network = FlowNetwork {
            id: 1,
            s1: [2].into_iter().collect(),
            s2: [3].into_iter().collect(),
            t_in: [4].into_iter().collect(),
            t_out: VertexSet::default(),
            arcs: VertexMap::default(),
        };
        network.arcs.insert(1, [2].into_iter().collect());
        network.arcs.insert(2, [3].into_iter().collect());
        network.arcs.insert(3, [4].into_iter().collect());
        let mut vias = Vias::new(10);
        let mut four_vias: VertexMap<VertexSet> = VertexMap::default();
        four_vias.insert(4, [5, 6].into_iter().collect());
        vias.vias.insert(2, four_vias);

        network.add_vias_from_s1(&vias);

        assert_eq!(
            *network.arcs.get(&2).unwrap(),
            [5, 3].into_iter().collect()
        );
        assert_eq!(*network.arcs.get(&5).unwrap(), [4].into_iter().collect());
        //ensures only one of such edge is added for a target
        assert!(!network.arcs.contains_key(&6));
    }

    #[test]
    fn add_extra_targets_add_direct_targets_from_s() {
        let mut network = FlowNetwork {
            id: 1,
            s1: [2].into_iter().collect(),
            s2: [3].into_iter().collect(),
            t_in: [4].into_iter().collect(),
            t_out: VertexSet::default(),
            arcs: VertexMap::default(),
        };
        network.arcs.insert(1, [2].into_iter().collect());
        network.arcs.insert(2, [3].into_iter().collect());
        network.arcs.insert(3, [4].into_iter().collect());
        let targets: VertexSet = [4, 5, 6, 7, 8].into_iter().collect();
        let mut adm: VertexMap<AdmData> = VertexMap::default();
        adm.insert(
            3,
            AdmData::new(3, [4, 5, 6].into_iter().collect()),
        );
        adm.insert(
            2,
            AdmData::new(2, [7, 8].into_iter().collect()),
        );

        network.add_extra_targets(&adm, &targets, &Vias::new(10), &targets.clone());

        assert_eq!(
            *network.arcs.get(&3).unwrap(),
            [4, 5].into_iter().collect()
        );
        assert_eq!(
            *network.arcs.get(&2).unwrap(),
            [3, 8].into_iter().collect()
        );
    }

    #[test]
    fn add_extra_targets_add_indirect_targets_from_s1() {
        let mut network = FlowNetwork {
            id: 1,
            s1: [2].into_iter().collect(),
            s2: [3].into_iter().collect(),
            t_in: [4].into_iter().collect(),
            t_out: VertexSet::default(),
            arcs: VertexMap::default(),
        };
        network.arcs.insert(1, [2].into_iter().collect());
        network.arcs.insert(2, [3].into_iter().collect());
        network.arcs.insert(3, [4].into_iter().collect());
        let targets: VertexSet = [4, 6, 8].into_iter().collect();
        let mut adm: VertexMap<AdmData> = VertexMap::default();
        adm.insert(
            3,
            AdmData::new(3, [4].into_iter().collect()),
        );
        adm.insert(2, AdmData::new(2, VertexSet::default()));
        let mut vias = Vias::new(10);
        let mut two_vias: VertexMap<VertexSet> = VertexMap::default();
        two_vias.insert(8, [7].into_iter().collect());
        two_vias.insert(6, [5].into_iter().collect());
        vias.vias.insert(2, two_vias);

        network.add_extra_targets(&adm, &targets, &vias, &targets.clone());

        assert_eq!(*network.arcs.get(&3).unwrap(), [4].into_iter().collect());
        assert_eq!(
            *network.arcs.get(&2).unwrap(),
            [3, 7].into_iter().collect()
        );
        assert_eq!(*network.arcs.get(&7).unwrap(), [8].into_iter().collect());
        //ensures only one of such edge is added for a target
        assert!(!network.arcs.contains_key(&5));
    }

    #[test]
    fn add_extra_targets_does_not_add_indirect_targets_from_s1_if_a_direct_extra_target_exists() {
        let mut network = FlowNetwork {
            id: 1,
            s1: [2].into_iter().collect(),
            s2: [3].into_iter().collect(),
            t_in: [4].into_iter().collect(),
            t_out: VertexSet::default(),
            arcs: VertexMap::default(),
        };
        network.arcs.insert(1, [2].into_iter().collect());
        network.arcs.insert(2, [3].into_iter().collect());
        network.arcs.insert(3, [4].into_iter().collect());
        let targets: VertexSet = [4, 6, 7].into_iter().collect();
        let mut adm: VertexMap<AdmData> = VertexMap::default();
        adm.insert(
            3,
            AdmData::new(3, [4].into_iter().collect()),
        );
        adm.insert(
            2,
            AdmData::new(2, [7].into_iter().collect()),
        );
        let mut vias = Vias::new(10);
        let mut two_vias: VertexMap<VertexSet> = VertexMap::default();
        two_vias.insert(6, [5].into_iter().collect());
        vias.vias.insert(2, two_vias);

        network.add_extra_targets(&adm, &targets, &vias, &targets.clone());

        assert_eq!(*network.arcs.get(&3).unwrap(), [4].into_iter().collect());
        assert_eq!(
            *network.arcs.get(&2).unwrap(),
            [3, 7].into_iter().collect()
        );
        //ensures only one of such edge is added for a target
        assert!(!network.arcs.contains_key(&5));
    }

    #[test]
    fn split_edges_in_network_creates_new_graph_with_split_edges() {
        let mut network = FlowNetwork {
            id: 1,
            s1: [2, 6].into_iter().collect(),
            s2: [3].into_iter().collect(),
            t_in: [4].into_iter().collect(),
            t_out: [5].into_iter().collect(),
            arcs: VertexMap::default(),
        };
        network
            .arcs
            .insert(1, [2, 6].into_iter().collect());
        network
            .arcs
            .insert(2, [3, 5].into_iter().collect());
        network.arcs.insert(3, [4].into_iter().collect());
        network.arcs.insert(6, [4].into_iter().collect());

        let result = network.split_edges_in_network();

        assert!(!result.contains_key(&(-1)));
        assert_eq!(*result.get(&1).unwrap(), [-2, -6].into_iter().collect());
        assert_eq!(*result.get(&(-2)).unwrap(), [2].into_iter().collect());
        assert_eq!(*result.get(&(-6)).unwrap(), [6].into_iter().collect());
        assert_eq!(*result.get(&2).unwrap(), [-3, -5].into_iter().collect());
        assert_eq!(*result.get(&(-5)).unwrap(), [5].into_iter().collect());
        assert_eq!(*result.get(&3).unwrap(), [-4].into_iter().collect());
        assert_eq!(*result.get(&(-3)).unwrap(), [3].into_iter().collect());
        assert_eq!(*result.get(&(-4)).unwrap(), [4].into_iter().collect());
        assert_eq!(*result.get(&6).unwrap(), [-4].into_iter().collect());
    }

    #[test]
    fn set_edges_direction_reverses_edges_in_packing() {
        let u = 1;
        let network = FlowNetwork::new(u);
        let mut u_adm_data = AdmData::new(u, VertexSet::default());
        u_adm_data
            .packing
            .insert(4, Path::ThreePath(2, 3, 4));
        let mut flow: HashMap<i32, HashSet<i32>> = HashMap::default();
        flow.insert(1, [-2, -6].into_iter().collect());
        flow.insert(-2, [2].into_iter().collect());
        flow.insert(-6, [6].into_iter().collect());
        flow.insert(2, [-3, -5].into_iter().collect());
        flow.insert(-5, [5].into_iter().collect());
        flow.insert(-3, [3].into_iter().collect());
        flow.insert(3, [-4].into_iter().collect());
        flow.insert(-4, [4].into_iter().collect());
        flow.insert(6, [-4].into_iter().collect());

        let result = network.set_edges_direction(flow, &u_adm_data);

        assert_eq!(*result.get(&1).unwrap(), [-6].into_iter().collect());
        assert_eq!(*result.get(&(-2)).unwrap(), [1].into_iter().collect());
        assert_eq!(*result.get(&(-6)).unwrap(), [6].into_iter().collect());
        assert_eq!(*result.get(&2).unwrap(), [-2, -5].into_iter().collect());
        assert_eq!(*result.get(&(-5)).unwrap(), [5].into_iter().collect());
        assert_eq!(*result.get(&(-3)).unwrap(), [2].into_iter().collect());
        assert_eq!(*result.get(&3).unwrap(), [-3].into_iter().collect());
        assert_eq!(*result.get(&(-4)).unwrap(), [3].into_iter().collect());
        assert_eq!(*result.get(&6).unwrap(), [-4].into_iter().collect());
        assert_eq!(*result.get(&4).unwrap(), [-4].into_iter().collect());
    }

    #[test]
    fn bfs_returns_path_if_there_is_a_path_from_root_to_t_out() {
        let mut network = FlowNetwork::new(1);
        network.t_out.insert(10);
        let mut flow: HashMap<i32, HashSet<i32>> = HashMap::default();
        flow.insert(1, [2].into_iter().collect());
        flow.insert(2, [3].into_iter().collect());
        flow.insert(3, [4].into_iter().collect());
        flow.insert(4, [5].into_iter().collect());
        flow.insert(5, [6, 9].into_iter().collect());
        flow.insert(6, [1].into_iter().collect());
        flow.insert(9, [8].into_iter().collect());
        flow.insert(8, [7, 10].into_iter().collect());
        flow.insert(7, [1].into_iter().collect());

        let result = network.bfs(2, &flow).unwrap();

        assert_eq!(result, vec![1, 2, 3, 4, 5, 9, 8, 10]);
    }

    #[test]
    fn bfs_returns_none_if_there_is_not_a_path_from_root_to_t_out() {
        let mut network = FlowNetwork::new(1);
        let mut flow: HashMap<i32, HashSet<i32>> = HashMap::default();
        flow.insert(1, [5].into_iter().collect());
        flow.insert(5, [6].into_iter().collect());
        flow.insert(6, [4].into_iter().collect());
        flow.insert(4, [3].into_iter().collect());
        flow.insert(3, [2].into_iter().collect());
        flow.insert(2, [1].into_iter().collect());

        let result = network.bfs(5, &flow);

        assert!(result.is_none());
    }

    #[test]
    fn get_new_packing_edges_returns_correct_edges() {
        let mut network = FlowNetwork::new(1);
        network.t_in = [4].into_iter().collect();
        network.t_out = [6].into_iter().collect();
        network.arcs.insert(1, [2, 5].into_iter().collect());
        network.arcs.insert(2, [3, 6].into_iter().collect());
        network.arcs.insert(3, [4].into_iter().collect());
        network.arcs.insert(5, [3].into_iter().collect());
        let aug_path_vertices: VertexSet = [1, 5, 3, 2, 6].into_iter().collect();
        let mut flow: HashMap<i32, HashSet<i32>> = HashMap::default();
        flow.insert(1, [-5].into_iter().collect());
        flow.insert(-2, [1].into_iter().collect());
        flow.insert(2, [-2, -6].into_iter().collect());
        flow.insert(-3, [2].into_iter().collect());
        flow.insert(3, [-3].into_iter().collect());
        flow.insert(-4, [3].into_iter().collect());
        flow.insert(4, [-4].into_iter().collect());
        flow.insert(-5, [5].into_iter().collect());
        flow.insert(5, [-3].into_iter().collect());
        flow.insert(-6, [6].into_iter().collect());

        let (new_edges, new_t) = network.get_new_packing_edges(aug_path_vertices, flow);

        assert_eq!(new_t, [4, 6].into_iter().collect());
        assert_eq!(*new_edges.get(&2).unwrap(), 1);
        assert_eq!(*new_edges.get(&5).unwrap(), 1);
        assert_eq!(*new_edges.get(&6).unwrap(), 2);
        assert_eq!(*new_edges.get(&3).unwrap(), 5);
        assert_eq!(*new_edges.get(&4).unwrap(), 3);
    }
    
    #[test]
    fn update_packing_updates_packing_correctly(){
        let network = FlowNetwork::new(1);
        let mut u_adm_data = AdmData::new(1, VertexSet::default());
        let mut new_edges : VertexMap<Vertex> = VertexMap::default();
        let new_t : VertexSet = [4,6].into_iter().collect();
        new_edges.insert(4,3);
        new_edges.insert(3,2);
        new_edges.insert(2,1);
        new_edges.insert(5,1);
        new_edges.insert(6,5);

        network.update_packing(&mut u_adm_data, new_edges, &new_t);

        assert_eq!(*u_adm_data.packing.get(&6).unwrap(), Path::TwoPath(5, 6));
        assert_eq!(*u_adm_data.packing.get(&4).unwrap(), Path::ThreePath(2,3,4));
    }
}
