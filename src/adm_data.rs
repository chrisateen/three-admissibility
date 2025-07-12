use graphbench::graph::{Vertex, VertexMap, VertexSet};

pub struct AdmData {
    pub id: Vertex,
    pub n_in_r: VertexSet,
    pub t1: VertexSet,
    pub t2: VertexSet,
    pub t3: VertexSet,
    pub packing: VertexMap<Vec<Vertex>>,
}

impl AdmData {
    pub fn new(v: Vertex, v_neighbours: VertexSet) -> Self {
        AdmData {
            id: v,
            n_in_r: VertexSet::default(),
            t1: v_neighbours,
            t2: VertexSet::default(),
            t3: VertexSet::default(),
            packing: VertexMap::default(),
        }
    }

    pub fn can_add_t2_path_to_pack(&mut self, t2: &Vertex, t1: &Vertex) -> bool {
        !self.packing.contains_key(t2) & !self.t1.contains(t1)
    }

    // pub fn can_add_t3_path_to_pack(&mut self, t3: &Vertex, t2: &Vertex, t1: &Vertex) -> bool {
    //     !self.packing.contains_key(t3) & !self.t1.contains(t1) & !self.t2.contains(t2)
    // }
    pub fn add_t2_to_packing(&mut self, t2: &Vertex, t1: &Vertex) {
        self.packing.insert(*t2, vec![*t1]);
        self.t1.insert(*t1);
        self.t2.insert(*t2);
    }

    pub fn add_t3_to_packing(&mut self, t3: &Vertex, t1: &Vertex,  t2: &Vertex) {
        self.packing.insert(*t3, vec![*t1, *t2]);
        self.t1.insert(*t1);
        self.t2.insert(*t2);
        self.t3.insert(*t3);
    }

    pub fn delete_packing(&mut self) {
        self.packing.clear();
        self.t2.clear();
        self.t3.clear();
    }
    
    pub fn is_v_in_pack(&self, v: &Vertex) -> bool {
        self.t1.contains(v) || self.t2.contains(v) || self.t3.contains(v)
    }

    pub fn remove_v_from_packing(&mut self, v: &Vertex) ->  Vec<Vertex> {
        if self.t1.contains(v) & !self.n_in_r.contains(v){
            self.t1.remove(v);
            self.n_in_r.insert(*v);
            return vec![];
        }

        let p = self.packing.remove(v).unwrap();
        self.t1.remove(&p[0]);
        if p.len() == 1 {
            self.t2.remove(v);
        }else {
            self.t2.remove(&p[1]);
            self.t3.remove(v);
        }
        p
    }

    pub fn size_of_packing(&self) -> usize{
        let num_t1_l = self.t1.difference(&self.n_in_r).count();
        num_t1_l + self.packing.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use graphbench::graph::{VertexSet};

    fn vertex(v: usize) -> Vertex {
        v as Vertex
    }

    #[test]
    fn can_add_t2_path_to_pack_returns_true_if_vertices_are_not_in_pack() {
        let mut adm = AdmData::new(vertex(1),VertexSet::default());;
        adm.packing.insert(vertex(3), vec![2]);
        adm.t1.insert(vertex(2));
        adm.t2.insert(vertex(3));

        assert_eq!(adm.can_add_t2_path_to_pack(&vertex(5), &vertex(4)), true);
        assert_eq!(adm.can_add_t2_path_to_pack(&vertex(7), &vertex(6)), true);
    }

    #[test]
    fn can_add_t2_path_to_pack_returns_false_if_vertices_are_in_pack() {
        let mut adm = AdmData::new(vertex(1),VertexSet::default());;
        adm.packing.insert(vertex(3), vec![2]);
        adm.t1.insert(vertex(2));
        adm.t2.insert(vertex(3));

        assert_eq!(adm.can_add_t2_path_to_pack(&vertex(3), &vertex(4)), false); // 3 is already in t2
        assert_eq!(adm.can_add_t2_path_to_pack(&vertex(5), &vertex(2)), false); // 2 is already in t1
    }

    #[test]
    fn add_t2_to_packing_updates_correctly() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());

        adm.add_t2_to_packing(&vertex(2), &vertex(3));

        assert_eq!(adm.packing.len(),1);
        assert_eq!(adm.packing.contains_key(&vertex(2)), true);
        assert_eq!(adm.packing.get(&vertex(2)), Some(&vec![3]));
        assert_eq!(adm.t1.contains(&vertex(3)), true);
        assert_eq!(adm.t2.contains(&vertex(2)), true);
    }

    #[test]
    fn add_t3_to_packing_updates_correctly() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());

        adm.add_t3_to_packing(&vertex(4), &vertex(5),&vertex(6));

        assert_eq!(adm.packing.len(),1);
        assert_eq!(adm.packing.contains_key(&vertex(4)), true);
        assert_eq!(adm.packing.get(&vertex(4)), Some(&vec![5,6]));
        assert_eq!(adm.t1.contains(&vertex(5)), true);
        assert_eq!(adm.t2.contains(&vertex(6)), true);
    }

    #[test]
    fn delete_packing_clears_packing_t2_and_t3() {
        let mut adm = AdmData::new(vertex(1),VertexSet::default());
        adm.packing.insert(vertex(4), vec![2,3]);
        adm.t1.insert(vertex(2));
        adm.t2.insert(vertex(3));
        adm.t3.insert(vertex(4));

        adm.delete_packing();

        assert_eq!(adm.packing.is_empty(), true);
        assert_eq!(adm.t2.is_empty(), true);
        assert_eq!(adm.t3.is_empty(), true);
    }

    #[test]
    fn test_is_an_endpoint_in_pack_returns_true_if_v_in_t1_and_in_l() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.t1.insert(vertex(2));

        assert_eq!(adm.is_an_endpoint_in_pack(&vertex(2)), true);
    }

    #[test]
    fn test_is_an_endpoint_in_pack_returns_true_if_v_is_a_key_in_pack() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.packing.insert(vertex(4), vec![2,3]);

        assert_eq!(adm.is_an_endpoint_in_pack(&vertex(4)), true);
    }

    #[test]
    fn test_is_an_endpoint_in_pack_returns_false_if_v_in_t1_and_in_r() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.t1.insert(vertex(2));
        adm.n_in_r.insert(vertex(2));

        assert_eq!(adm.is_an_endpoint_in_pack(&vertex(2)), false);
    }

    #[test]
    fn remove_v_from_packing_removes_t1_vertex() {
        let t1_l : VertexSet = vec![2,3].into_iter().collect();
        let mut adm = AdmData::new(vertex(1), t1_l);

        let removed = adm.remove_v_from_packing(&vertex(2));
        assert_eq!(removed, vec![]);
        assert_eq!(adm.t1.contains(&vertex(2)), false);
    }

    #[test]
    fn remove_v_from_packing_removes_t2_vertex() {
        let mut adm = AdmData::new(vertex(1),VertexSet::default());;
        adm.packing.insert(vertex(3), vec![2]);
        adm.t1.insert(vertex(2));
        adm.t2.insert(vertex(3));

        let removed = adm.remove_v_from_packing(&vertex(3));
        assert_eq!(removed, vec![2]);
        assert_eq!(adm.t1.contains(&vertex(2)), false);
        assert_eq!(adm.t2.contains(&vertex(3)), false);
    }

    #[test]
    fn remove_v_from_packing_removes_t3_vertex() {
        let mut adm = AdmData::new(vertex(1),VertexSet::default());;
        adm.packing.insert(vertex(4), vec![5,6]);
        adm.t1.insert(vertex(5));
        adm.t2.insert(vertex(6));

        let removed = adm.remove_v_from_packing(&vertex(4));
        assert_eq!(removed, vec![5,6]);
        assert_eq!(adm.t1.contains(&vertex(5)), false);
        assert_eq!(adm.t2.contains(&vertex(6)), false);
        assert_eq!(adm.t3.contains(&vertex(4)), false);
    }

    #[test]
    fn test_size_of_packing() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.t1 = vec![2,3,4,5].iter().cloned().collect();
        adm.n_in_r = vec![2,3].iter().cloned().collect();
        adm.packing.insert(vertex(6),vec![2]);
        adm.packing.insert(vertex(8), vec![3,7]);

        //vertex 4,5 in t1 and not in r, vertex 6,7 in t3
        assert_eq!(adm.size_of_packing(), 4);
    }
}