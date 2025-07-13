use graphbench::editgraph::EditGraph;
use graphbench::graph::{Graph, Vertex, VertexMap, VertexSet};
use crate::vias::Vias;

pub struct AugmentingPath {
    pub id: Vertex,
    pub edges: VertexMap<VertexSet>,
    pub s1: VertexSet,
    pub s2: VertexSet,
}

impl AugmentingPath {
    pub fn new(id: Vertex) -> Self {
        AugmentingPath{
            id,
            edges: VertexMap::default(),
            s1: VertexSet::default(),
            s2: VertexSet::default(),
        }
    }

    /*
        Add edges from pack of self.id
        pack edges added in a direction that points away from self.id/root
    */
    pub fn add_pack_edges(&mut self, packing: &VertexMap<Vec<Vertex>>) {
        let mut root_neighbours = self.edges.remove(&self.id).unwrap_or_default();

        for w in packing.keys() {
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
        Add all edges between vertices in packing and the targets
    */
    pub fn add_edges_between(&mut self, t2_t3: &VertexSet, graph: &EditGraph) {
        for v in &self.s1 {
            for w in self.s2.union(t2_t3) {
                if graph.adjacent(v,w){
                    //edges are added as s2 -> s1 or target -> s1
                    let w_neighbours = self.edges.entry(*w).or_default();
                    w_neighbours.insert(*v);
                }
            }
        }

        for v in &self.s2 {
            for w in t2_t3 {
                if graph.adjacent(v,w){
                    let w_neighbours = self.edges.entry(*w).or_default();
                    w_neighbours.insert(*v);
                }
            }
        }
    }

    /*
        Add outside edges using vias
    */
    pub fn add_outside_edges(&mut self, t2_t3: &VertexSet, graph: &EditGraph, vias: &Vias) {
        let mut z1 = VertexSet::default();
        let mut z2 = VertexSet::default();

        for x in self.s2.union(t2_t3){
            let x_u_vias = vias.get_vias(*x, self.id);
            if x_u_vias.is_none(){
                continue;
            }

            //TODO should we loop until we find the first eligible via
            let mut eligible_vias = x_u_vias
                .unwrap()
                .iter()
                .filter( |x| !self.s1.contains(x) && !self.s2.contains(x) && !z1.contains(x));
            z1.insert(*eligible_vias.next().unwrap());
        }

        for y in t2_t3 {
            //TODO do we loop over all N(u) cap R
        }

        //TODO do we need to add vertices between each R vertex
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