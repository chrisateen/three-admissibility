use std::collections::HashSet;
use graphbench::editgraph::EditGraph;
use graphbench::graph::{Graph, Vertex, VertexMap, VertexSet};
use crate::vias::Vias;

pub struct AugmentingPath {
    pub id: Vertex,
    pub edges: VertexMap<VertexSet>,
    pub s1: VertexSet,
    pub s2: VertexSet,
    pub t: VertexSet,
}

impl AugmentingPath {
    pub fn new(id: Vertex) -> Self {
        AugmentingPath{
            id,
            edges: VertexMap::default(),
            s1: VertexSet::default(),
            s2: VertexSet::default(),
            t: VertexSet::default(),
        }
    }

    /*
        Add edges from pack of self.id
        edges added in a direction that points away from root
    */
    pub fn add_pack_edges(&mut self, packing: &VertexMap<Vec<Vertex>>) {
        let mut root_neighbours = self.edges.remove(&self.id).unwrap_or_default();

        for w in packing.keys() {
            self.t.insert(*w);
            let edge = packing.get(w).unwrap();
            if edge.is_empty() || edge.len() > 2 {
                panic!("Invalid edge length of {} in pack of {}", edge.len(), self.id);
            }
            let t1 = edge.first().unwrap();
            let t1_neighbours = self.edges.entry(*t1).or_default();
            let t2 = edge.get(1);

            root_neighbours.insert(*t1);
            self.s1.insert(*t1);

            match t2 {
                Some(t2) => {
                    t1_neighbours.insert(*t2);
                    let t2_neighbours = self.edges.entry(*t2).or_default();
                    t2_neighbours.insert(*w);
                    self.s2.insert(*t2);
                }
                None => {
                    t1_neighbours.insert(*w);
                }
            };
        }

        self.edges.insert(self.id, root_neighbours);
    }

    /*
        Add all edges between vertices in packing
        edges added in a direction that points away from root
    */
    pub fn add_edges_between(&mut self, graph: &EditGraph) {
        let s2_union_t: VertexSet = self.s2.union(&self.t).cloned().collect();
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
            for w in &self.t {
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
    pub fn add_edges_containing_n_of_r(&mut self, n_in_r: VertexSet,  graph: &EditGraph, vias: &Vias) {
        let s2_union_t: VertexSet = self.s2.union(&self.t).cloned().collect();
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
                if self.t.contains(w){
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
        TODO: If t in t2 should we be adding this?
    */
    pub fn add_vias_from_s1(&mut self, vias: &Vias) {
        let s : VertexSet = self.s1.union(&self.s2).cloned().collect();

        for v in &self.s1 {
            let mut v_neighbours = self.edges.remove(v).unwrap_or_default();

            for w in &self.t {
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

}

#[cfg(test)]
mod test_augmenting_path {
    use super::*;

    fn vertex(v: usize) -> Vertex {
        v as Vertex
    }

    #[test]
    fn add_pack_edges_add_edges_in_correct_order() {
        let mut packing : VertexMap<Vec<Vertex>> = VertexMap::default();
        packing.insert(vertex(4), Vec::from([vertex(2), vertex(3)]));
        packing.insert(6, Vec::from([vertex(5)]));
        let mut augmenting_path = AugmentingPath::new(1);

        augmenting_path.add_pack_edges(&packing);

        assert_eq!(*augmenting_path.edges.get(&1).unwrap(), [2,5].into_iter().collect());
        assert_eq!(*augmenting_path.edges.get(&2).unwrap(), [3].into_iter().collect());
        assert_eq!(*augmenting_path.edges.get(&3).unwrap(), [4].into_iter().collect());
        assert_eq!(*augmenting_path.edges.get(&5).unwrap(), [6].into_iter().collect());
        assert!(!augmenting_path.edges.contains_key(&4));
        assert!(!augmenting_path.edges.contains_key(&6));
        assert_eq!(augmenting_path.s1, [2,5].into_iter().collect());
    }
}