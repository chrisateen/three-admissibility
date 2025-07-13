use graphbench::graph::{Vertex, VertexMap, VertexSet};

pub struct AugmentingPath {
    pub id: Vertex,
    pub edges: VertexMap<VertexSet>,
}

impl AugmentingPath {
    pub fn new(id: Vertex) -> Self {
        AugmentingPath{
            id,
            edges: VertexMap::default(),
        }
    }

    /*
        Add edges from pack of self.id
        pack edges added in a direction that points away from self.id/root
    */
    pub fn add_pack_edges(&mut self, packing: VertexMap<Vec<Vertex>>) {
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

            match t2 {
                Some(t2) => {
                    t1_neighbours.insert(*t2);
                    let t2_neighbours = self.edges.entry(*t2).or_default();
                    t2_neighbours.insert(*w);
                }
                None => {
                    t1_neighbours.insert(*w);
                }
            };
        }

        self.edges.insert(self.id, root_neighbours);
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

        augmenting_path.add_pack_edges(packing);

        assert_eq!(*augmenting_path.edges.get(&1).unwrap(), [2,5].into_iter().collect());
        assert_eq!(*augmenting_path.edges.get(&2).unwrap(), [3].into_iter().collect());
        assert_eq!(*augmenting_path.edges.get(&3).unwrap(), [4].into_iter().collect());
        assert_eq!(*augmenting_path.edges.get(&5).unwrap(), [6].into_iter().collect());
        assert!(!augmenting_path.edges.contains_key(&4));
        assert!(!augmenting_path.edges.contains_key(&6));
    }
}