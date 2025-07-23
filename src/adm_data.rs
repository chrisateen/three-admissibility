use graphbench::graph::{Graph, Vertex, VertexMap, VertexSet};

use crate::adm_graph::AdmGraph;
use crate::debug_println;

#[derive(Debug, PartialEq)]
pub enum Path {
    TwoPath(Vertex, Vertex),           //t1_r. t2_l
    ThreePath(Vertex, Vertex, Vertex), //t1_r, t2_r, t3_l
}

pub struct AdmData {
    pub id: Vertex,
    pub n_in_r: VertexSet,
    pub t1: VertexSet,
    pub t2: VertexSet,
    pub t3: VertexSet,
    pub packing: VertexMap<Path>,
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

    pub fn debug_check_consistency(&self, adm_graph: &AdmGraph) {
        let neighbours: VertexSet = adm_graph.graph.neighbours(&self.id).cloned().collect();
        let left_neighbours: VertexSet = neighbours
            .iter()
            .filter(|x| adm_graph.l.contains(x))
            .cloned()
            .collect();

        // Left neighbours of this vertex should all be contained in t1
        debug_assert_eq!(left_neighbours.difference(&self.t1).count(), 0);

        for x in self.t1.iter() {
            debug_assert!(neighbours.contains(x));
        }
        for x in self.n_in_r.iter() {
            debug_assert!(neighbours.contains(x));
        }

        if adm_graph.l.contains(&self.id) && !adm_graph.candidates.contains(&self.id) {
            // This is a 3-packing
            // Ensure that vertices of packaging are on correct sides
            let mut all_good = true;
            for path in self.packing.values() {
                match path {
                    Path::TwoPath(s1, t2) => {
                        all_good &= adm_graph.r.contains(s1);
                        all_good &= adm_graph.l.contains(t2);
                    }
                    Path::ThreePath(s1, s2, t3) => {
                        all_good &= adm_graph.r.contains(s1);
                        all_good &= adm_graph.r.contains(s2);
                        all_good &= adm_graph.l.contains(t3);
                    }
                }
            }

            if !all_good {
                debug_println!("Inconsistency in packing rooted at {}", self.id);
                debug_println!("Left = {:?}", adm_graph.l);
                debug_println!("Candidates = {:?}", adm_graph.candidates);
                debug_println!("Right = {:?}", adm_graph.r);
                for path in self.packing.values() {
                    match path {
                        Path::TwoPath(s1, t2) => {
                            debug_println!("Path {}-{}-{}", self.id, s1, t2);
                            debug_println!(
                                "   vertex {} is in R: {}",
                                s1,
                                adm_graph.r.contains(s1)
                            );
                            debug_println!(
                                "   vertex {} is in L: {}",
                                t2,
                                adm_graph.l.contains(t2)
                            );
                        }
                        Path::ThreePath(s1, s2, t3) => {
                            debug_println!("Path {}-{}-{}-{}", self.id, s1, s2, t3);
                            debug_println!(
                                "   vertex {} is in R: {}",
                                s1,
                                adm_graph.r.contains(s1)
                            );
                            debug_println!(
                                "   vertex {} is in R: {}",
                                s2,
                                adm_graph.r.contains(s2)
                            );
                            debug_println!(
                                "   vertex {} is in L: {}",
                                t3,
                                adm_graph.l.contains(t3)
                            );
                        }
                    }
                }

                panic!();
            }
        } else if adm_graph.r.contains(&self.id) {
            debug_assert_eq!(self.t3.len(), 0); // There should be no T3 vertices
        }
    }

    pub fn add_t2_to_packing(&mut self, t1: &Vertex, t2: &Vertex) {
        self.packing.insert(*t2, Path::TwoPath(*t1, *t2));

        debug_assert_ne!(*t1, self.id);
        debug_assert_ne!(*t2, self.id);
        debug_assert!(!self.t1.contains(t1));
        debug_assert!(!self.t2.contains(t1));
        debug_assert!(!self.t1.contains(t2));
        debug_assert!(!self.t2.contains(t2));

        self.t1.insert(*t1);
        self.t2.insert(*t2);
    }

    pub fn add_t3_to_packing(&mut self, t1: &Vertex, t2: &Vertex, t3: &Vertex) {
        self.packing.insert(*t3, Path::ThreePath(*t1, *t2, *t3));

        debug_assert_ne!(*t1, self.id);
        debug_assert_ne!(*t2, self.id);
        debug_assert_ne!(*t3, self.id);
        debug_assert!(!self.t1.contains(t1));
        debug_assert!(!self.t2.contains(t1));
        debug_assert!(!self.t3.contains(t1));
        debug_assert!(!self.t1.contains(t2));
        debug_assert!(!self.t2.contains(t2));
        debug_assert!(!self.t3.contains(t2));
        debug_assert!(!self.t1.contains(t3));
        debug_assert!(!self.t2.contains(t3));
        debug_assert!(!self.t3.contains(t3));

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

    pub fn remove_v_from_packing(&mut self, v: &Vertex) -> Option<Path> {
        if self.t1.contains(v) & !self.n_in_r.contains(v) {
            self.t1.remove(v);
            self.n_in_r.insert(*v);
            return None;
        }

        let path = self.packing.remove(v).unwrap();

        match path {
            Path::TwoPath(x, y) => {
                self.t1.remove(&x);
                self.t2.remove(&y);
            }
            Path::ThreePath(u, w, y) => {
                self.t1.remove(&u);
                self.t2.remove(&w);
                self.t3.remove(&y);
            }
        }
        Some(path)
    }

    pub fn size_of_packing(&self) -> usize {
        let num_t1_l = self.t1.difference(&self.n_in_r).count();
        num_t1_l + self.packing.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use graphbench::graph::VertexSet;

    fn vertex(v: usize) -> Vertex {
        v as Vertex
    }

    #[test]
    fn add_t2_to_packing_updates_correctly() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());

        adm.add_t2_to_packing(&vertex(3), &vertex(2));

        assert_eq!(adm.packing.len(), 1);
        assert!(adm.packing.contains_key(&vertex(2)));
        assert_eq!(
            adm.packing.get(&vertex(2)),
            Some(&Path::TwoPath(vertex(3), vertex(2)))
        );
        assert!(adm.t1.contains(&vertex(3)));
        assert!(adm.t2.contains(&vertex(2)));
    }

    #[test]
    fn add_t3_to_packing_updates_correctly() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());

        adm.add_t3_to_packing(&vertex(5), &vertex(6), &vertex(4));

        assert_eq!(adm.packing.len(), 1);
        assert!(adm.packing.contains_key(&vertex(4)));
        assert_eq!(
            adm.packing.get(&vertex(4)),
            Some(&Path::ThreePath(vertex(5), vertex(6), vertex(4)))
        );
        assert!(adm.t1.contains(&vertex(5)));
        assert!(adm.t2.contains(&vertex(6)));
    }

    #[test]
    fn delete_packing_clears_packing_t2_and_t3() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.packing
            .insert(vertex(4), Path::ThreePath(vertex(2), vertex(3), vertex(4)));
        adm.t1.insert(vertex(2));
        adm.t2.insert(vertex(3));
        adm.t3.insert(vertex(4));

        adm.delete_packing();

        assert!(adm.packing.is_empty());
        assert!(adm.t2.is_empty());
        assert!(adm.t3.is_empty());
    }

    #[test]
    fn test_is_v_in_pack_returns_true_if_v_in_t1() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.t1.insert(vertex(2));

        assert!(adm.is_v_in_pack(&vertex(2)));
    }

    #[test]
    fn test_is_v_in_pack_returns_true_if_v_in_t2() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.t2.insert(vertex(2));

        assert!(adm.is_v_in_pack(&vertex(2)));
    }

    #[test]
    fn test_is_v_in_pack_returns_true_if_v_in_t3() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.t3.insert(vertex(2));

        assert!(adm.is_v_in_pack(&vertex(2)));
    }

    #[test]
    fn test_is_v_in_pack_returns_false_if_v_is_not_in_t() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.t1.insert(vertex(2));
        adm.t2.insert(vertex(3));
        adm.t3.insert(vertex(4));

        assert!(!adm.is_v_in_pack(&vertex(5)));
    }

    #[test]
    fn remove_v_from_packing_removes_t1_vertex() {
        let t1_l: VertexSet = vec![2, 3].into_iter().collect();
        let mut adm = AdmData::new(vertex(1), t1_l);

        let removed = adm.remove_v_from_packing(&vertex(2));
        assert_eq!(removed, None);
        assert!(!adm.t1.contains(&vertex(2)));
    }

    #[test]
    fn remove_v_from_packing_removes_t2_vertex() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.packing
            .insert(vertex(3), Path::TwoPath(vertex(2), vertex(3)));
        adm.t1.insert(vertex(2));
        adm.t2.insert(vertex(3));

        let removed = adm.remove_v_from_packing(&vertex(3));
        assert_eq!(removed, Some(Path::TwoPath(vertex(2), vertex(3))));
        assert!(!adm.t1.contains(&vertex(2)));
        assert!(!adm.t2.contains(&vertex(3)));
    }

    #[test]
    fn remove_v_from_packing_removes_t3_vertex() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.packing
            .insert(vertex(4), Path::ThreePath(vertex(5), vertex(6), vertex(4)));
        adm.t1.insert(vertex(5));
        adm.t2.insert(vertex(6));

        let removed = adm.remove_v_from_packing(&vertex(4));
        assert_eq!(
            removed,
            Some(Path::ThreePath(vertex(5), vertex(6), vertex(4)))
        );
        assert!(!adm.t1.contains(&vertex(5)));
        assert!(!adm.t2.contains(&vertex(6)));
        assert!(!adm.t3.contains(&vertex(4)));
    }

    #[test]
    fn test_size_of_packing() {
        let mut adm = AdmData::new(vertex(1), VertexSet::default());
        adm.t1 = [2, 3, 4, 5].iter().cloned().collect();
        adm.n_in_r = [2, 3].iter().cloned().collect();
        adm.packing
            .insert(vertex(6), Path::TwoPath(vertex(2), vertex(6)));
        adm.packing
            .insert(vertex(8), Path::ThreePath(vertex(8), vertex(3), vertex(7)));

        //vertex 4,5 in t1 and not in r, vertex 6,7 in t3
        assert_eq!(adm.size_of_packing(), 4);
    }
}
