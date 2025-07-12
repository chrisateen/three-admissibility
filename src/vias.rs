use graphbench::graph::{Vertex, VertexMap, VertexSet};

pub struct Vias {
    //Key vias[v_in_r][v_in_t2_l_of_v]
    pub vias: VertexMap<VertexMap<VertexSet>>,
    p: usize
}

impl Vias {
    pub fn new(p : usize) -> Self {
        Vias {
            vias: VertexMap::default(),
            p
        }
    }

    pub fn add_a_via(&mut self, v:  Vertex, t2_v: Vertex, via: Vertex) -> bool {
        let v_entries = self.vias.entry(v).or_insert(VertexMap::default());
        let t2_v_entries = v_entries.entry(t2_v).or_insert(VertexSet::default());
        
        let max_vias = (2*self.p + 1);
        
        if t2_v_entries.len() > max_vias{
            panic!("Number of vias for {} is too large", t2_v);
        }
        if t2_v_entries.len() == max_vias{
            return false;
        }
        t2_v_entries.insert(via);
        true
    }

    //TODO to delete vias for a vertex v moving to R, we need to identify all the t2 of v
    // and for each of them delete v from vias
    pub fn remove_vias(&mut self, v: Vertex) {

    }

    pub fn get_vias(&self, v:  Vertex, t2_v: Vertex) -> Option<&VertexSet> {
        let v_entries = self.vias.get(&v);
        match v_entries {
            Some(v_entries) => {
                let vias = v_entries.get(&t2_v);
                match vias {
                    Some(vias) => Some(vias),
                    None => None
                }
            }
            None => None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use graphbench::graph::{Vertex, VertexMap, VertexSet};

    #[test]
    fn add_a_via_adds_a_via() {
        let mut vias = Vias::new(1);
        let v: Vertex = 1;
        let t2_v: Vertex = 3;
        let via: Vertex = 2;

        let via_added = vias.add_a_via(v, t2_v, via);
        let result = vias.get_vias(v, t2_v).unwrap();

        assert_eq!(via_added, true);
        assert_eq!(result.contains(&via), true);
    }

    #[test]
    fn add_a_via_does_not_add_a_via_if_t2_has_enough_vias() {
        let mut vias = Vias::new(1);
        let v: Vertex = 1;
        let t2_v: Vertex = 3;
        let mut t2_vias : VertexMap<VertexSet> = VertexMap::default();
        t2_vias.insert(t2_v, [4,5,6].into_iter().collect());
        vias.vias.insert(v, t2_vias);
        let via: Vertex = 2;

        let via_added = vias.add_a_via(v, t2_v, via);
        let result = vias.vias.get(&v).unwrap().get(&t2_v).unwrap();

        assert_eq!(via_added, false);
        assert_eq!(result.contains(&via), false);
    }

    #[test]
    #[should_panic(expected = "Number of vias for 3 is too large")]
    fn add_a_via_throws_if_t2_has_too_much_vias() {
        let mut vias = Vias::new(1);
        let v: Vertex = 1;
        let t2_v: Vertex = 3;
        let mut t2_vias : VertexMap<VertexSet> = VertexMap::default();
        t2_vias.insert(t2_v, [4,5,6,7].into_iter().collect());
        vias.vias.insert(v, t2_vias);
        let via: Vertex = 2;

        vias.add_a_via(v, t2_v, via);
    }

    #[test]
    fn get_vias_returns_vias_for_v(){
        let mut vias = Vias::new(1);
        let v: Vertex = 1;
        let t2_v: Vertex = 3;
        let expected: VertexSet = [4,5,6,7].into_iter().collect();
        let mut t2_vias : VertexMap<VertexSet> = VertexMap::default();
        t2_vias.insert(t2_v, expected.clone());
        vias.vias.insert(v, t2_vias);

        let result = vias.get_vias(v, t2_v).unwrap();
        assert_eq!(result.clone(), expected);
    }

    #[test]
    fn get_vias_returns_none_if_no_vias_for_t2(){
        let vias = Vias::new(1);
        let v: Vertex = 1;
        let t2_v: Vertex = 3;

        let result = vias.get_vias(v, t2_v);
        assert_eq!(result, None);
    }

    #[test]
    fn get_vias_returns_none_if_no_vias_for_v(){
        let mut vias = Vias::new(1);
        let v: Vertex = 1;
        let t2_v: Vertex = 3;
        let mut t2_vias : VertexMap<VertexSet> = VertexMap::default();
        t2_vias.insert(t2_v, [4,5,6,7].into_iter().collect());
        vias.vias.insert(2 as Vertex, t2_vias);

        let result = vias.get_vias(v, t2_v);
        assert_eq!(result, None);
    }
}