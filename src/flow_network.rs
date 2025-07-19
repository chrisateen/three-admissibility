use graphbench::editgraph::EditGraph;
use graphbench::graph::{Graph, Vertex, VertexMap, VertexSet};
use crate::adm_data::{AdmData, Path};
use crate::adm_graph::AdmGraph;
use crate::vias::Vias;

pub struct FlowNetwork {
    pub id: Vertex,
    pub edges: VertexMap<VertexSet>,
    pub s1: VertexSet,
    pub s2: VertexSet,
    pub t_in: VertexSet, // targets in packing
    pub t_out: VertexSet, //target not in packing
}

impl FlowNetwork {
    pub fn new(id: Vertex) -> Self {
        FlowNetwork {
            id,
            edges: VertexMap::default(),
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
        let mut root_neighbours = self.edges.remove(&self.id).unwrap_or_default();

        for w in packing.keys() {
            self.t_in.insert(*w);
            let path = packing.get(w).unwrap();

            match path {
                Path::TwoPath(s1, t2) => {
                    let s1_neighbours = self.edges.entry(*s1).or_default();
                    root_neighbours.insert(*s1);
                    self.s1.insert(*s1);
                    s1_neighbours.insert(*t2);
                },
                Path::ThreePath(s1, s2 , t3) => {
                    let mut s1_neighbours = self.edges.remove(s1).unwrap_or_default();
                    let mut s2_neighbours = self.edges.remove(s2).unwrap_or_default();

                    root_neighbours.insert(*s1);
                    self.s1.insert(*s1);
                    s1_neighbours.insert(*s2);
                    self.s2.insert(*s2);
                    s2_neighbours.insert(*t3);

                    self.edges.insert(*s1, s1_neighbours);
                    self.edges.insert(*s2, s2_neighbours);
                },
            }
        }

        self.edges.insert(self.id, root_neighbours);
    }

    /*
        Add all edges between vertices in packing
        edges added in a direction that points away from root
    */
    fn add_edges_between(&mut self, graph: &EditGraph) {
        let s2_union_t: VertexSet = self.s2.union(&self.t_in).cloned().collect();
        for v in &self.s1 {
            for w in &s2_union_t {
                if graph.adjacent(v,w){
                    //edges are added as s1 -> s2 or s1 -> t
                    let v_neighbours = self.edges.entry(*v).or_default();
                    v_neighbours.insert(*w);
                }
            }
        }

        for v in &self.s2 {
            for w in &self.t_in {
                if graph.adjacent(v,w){
                    //edges are added as s2 -> t
                    let v_neighbours = self.edges.entry(*v).or_default();
                    v_neighbours.insert(*w);
                }
            }
        }
    }

    /*
        Add all edges between N_R(v) that are not in the packing
         and vertices in pack that are in s2 and t
         Also add vertices between N_R(v) and t through a via
    */
    fn add_edges_containing_n_of_r(&mut self, n_in_r: &VertexSet,  graph: &EditGraph, vias: &Vias) {
        let s2_union_t: VertexSet = self.s2.union(&self.t_in).cloned().collect();
        let s : VertexSet = self.s1.union(&self.s2).cloned().collect();
        let mut root_neighbours = self.edges.remove(&self.id).unwrap_or_default();

        for v in n_in_r.difference(&self.s1) {
            let mut v_neighbours = self.edges.remove(v).unwrap_or_default();

            //Check if N_R(v) is connected directly to a vertex in s2_union_t
            for w in &s2_union_t {
                if graph.adjacent(v, w){
                    //edges are added as N_R(v) -> s2 or N_R(v) -> t
                    v_neighbours.insert(*w);
                    //add edge from root -> N_R(v) as well
                    root_neighbours.insert(*v);
                }

                //Check if we can add a via between N_R(v) and a target
                if self.t_in.contains(w){
                    let vias_v_w = vias.get_vias(*v,*w);
                    if vias_v_w.is_some() {
                        let eligible_via = vias_v_w.unwrap().difference(&s).next();
                        if let Some(x) = eligible_via {
                            let x_neighbours = self.edges.entry(*x).or_default();
                            //edges are added as N_R(v) -> via
                            //then via -> t
                            v_neighbours.insert(*x);
                            x_neighbours.insert(*w);
                            //add edge from root -> N_R(v) as well
                            root_neighbours.insert(*v);
                        }
                    }
                }
            }

            self.edges.insert(*v, v_neighbours);
        }

        self.edges.insert(self.id, root_neighbours);
    }

    /*
       Add a via between each s1 and target
    */
    fn add_vias_from_s1(&mut self, vias: &Vias) {
        let s : VertexSet = self.s1.union(&self.s2).cloned().collect();

        for v in &self.s1 {
            let mut v_neighbours = self.edges.remove(v).unwrap_or_default();

            for w in &self.t_in {
                let vias_v_w = vias.get_vias(*v, *w);

                if vias_v_w.is_some() {
                    let eligible_via = vias_v_w.unwrap().difference(&s).next();
                    if let Some(x) = eligible_via {
                        let x_neighbours = self.edges.entry(*x).or_default();
                        //edges are added as s1 -> via
                        v_neighbours.insert(*x);
                        //then via -> t
                        x_neighbours.insert(*w);
                    }
                }
            }
            self.edges.insert(*v, v_neighbours);
        }
    }

    /*
        For each s1 and s2 add a target that is not already in the packing
    */
    fn add_extra_targets(&mut self, adm_data: &VertexMap<AdmData>, targets: &VertexSet, vias: &Vias, l: &VertexSet) {
        let s1_s2 : VertexSet = self.s1.union(&self.s2).cloned().collect();

        'outer: for v in &s1_s2 {
            let v_adm_data = adm_data.get(v).unwrap();
            for w in v_adm_data.t1.intersection(l) {
                if !self.t_in.contains(w){
                    //adds edge s2 -> target outside packing or s1 -> target outside packing
                    self.edges.get_mut(v).unwrap().insert(*w);
                    self.t_out.insert(*w);
                    continue 'outer;
                }
            }

            if self.s2.contains(v) {
                continue 'outer;
            }

            for w in targets.difference(&self.t_in){
                let v_w_vias = vias.get_vias(*v, *w);
                if v_w_vias.is_none() {
                    continue;
                }
                let eligible_via = v_w_vias.unwrap().difference(&s1_s2).next();
                match eligible_via {
                    None => {}
                    Some(x) => {
                        //adds edge s1 (v) -> via (x)
                        self.edges.get_mut(v).unwrap().insert(*x);
                        self.edges.entry(*x).or_default().insert(*w);
                        self.t_out.insert(*w);
                        continue 'outer;
                    }
                }
            }
        }

    }

    pub fn construct_flow_network(&mut self, adm_graph: &AdmGraph, u: &AdmData, targets: &VertexSet){
        self.add_pack_edges(&u.packing);
        self.add_edges_between(adm_graph.graph);
        self.add_edges_containing_n_of_r(&u.n_in_r, adm_graph.graph, &adm_graph.vias);
        self.add_vias_from_s1(&adm_graph.vias);
        self.add_extra_targets(&adm_graph.adm_data, targets, &adm_graph.vias, &adm_graph.l);
    }

}

#[cfg(test)]
mod test_flow_network {
    use graphbench::graph::{EdgeSet, MutableGraph};
    use super::*;

    fn vertex(v: usize) -> Vertex {
        v as Vertex
    }

    fn create_test_graph(edges: EdgeSet) -> EditGraph {
        let mut graph = EditGraph::new();
        for (u, v) in edges.iter() {
            graph.add_edge(u, v);
        }

        graph
    }

    #[test]
    fn add_pack_edges_add_edges_in_correct_order() {
        let mut packing : VertexMap<Path> = VertexMap::default();
        packing.insert(vertex(4), Path::ThreePath(vertex(2), vertex(3), vertex(4)));
        packing.insert(6, Path::TwoPath(vertex(5), vertex(6)));
        let mut network = FlowNetwork::new(1);

        network.add_pack_edges(&packing);

        assert_eq!(*network.edges.get(&1).unwrap(), [2,5].into_iter().collect());
        assert_eq!(*network.edges.get(&2).unwrap(), [3].into_iter().collect());
        assert_eq!(*network.edges.get(&3).unwrap(), [4].into_iter().collect());
        assert_eq!(*network.edges.get(&5).unwrap(), [6].into_iter().collect());
        assert!(!network.edges.contains_key(&4));
        assert!(!network.edges.contains_key(&6));
        assert_eq!(network.s1, [2,5].into_iter().collect());
        assert_eq!(network.s2, [3].into_iter().collect());
        assert_eq!(network.t_in, [4,6].into_iter().collect());
    }

    #[test]
    fn add_edges_between_adds_edges() {
        let edges: EdgeSet = [(1, 2), (1, 5), (2, 3), (2,6),(3,4), (5,3), (5,6)]
            .iter()
            .cloned()
            .collect();
        let graph = create_test_graph(edges);
        let mut network = FlowNetwork {
            id: vertex(1),
            s1: [2,5].into_iter().collect(),
            s2: [3].into_iter().collect(),
            t_in: [4,6].into_iter().collect(),
            t_out: VertexSet::default(),
            edges: VertexMap::default(),
        };
        network.edges.insert(vertex(1), [2,5].into_iter().collect());
        network.edges.insert(vertex(2), [3].into_iter().collect());
        network.edges.insert(vertex(3), [4].into_iter().collect());
        network.edges.insert(vertex(5), [6].into_iter().collect());
        network.s1 = [2,5].into_iter().collect();
        network.s2 = [3].into_iter().collect();
        network.t_in = [4,6].into_iter().collect();

        network.add_edges_between(&graph);

        assert_eq!(*network.edges.get(&2).unwrap(), [3,6].into_iter().collect());
        assert_eq!(*network.edges.get(&5).unwrap(), [3,6].into_iter().collect());
    }
}