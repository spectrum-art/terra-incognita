//! Connected-component labeling on boolean pixel masks.
//!
//! Supports 4-connectivity and 8-connectivity. Returns a label grid
//! (0 = background, 1..N = component indices) and a list of component sizes.

/// Connectivity mode for connected-component labeling.
#[derive(Clone, Copy)]
pub enum Connectivity {
    Four,
    Eight,
}

/// Result of connected-component labeling.
pub struct Components {
    /// Label grid, row-major. 0 = background, 1..n = component labels.
    pub labels: Vec<u32>,
    /// Size (pixel count) of each component. Index 0 = component label 1, etc.
    pub sizes: Vec<usize>,
}

/// Label connected components in a boolean mask using BFS.
///
/// `mask[r * width + c]` is true for foreground pixels.
pub fn label_components(mask: &[bool], width: usize, height: usize, connectivity: Connectivity) -> Components {
    let n = width * height;
    let mut labels = vec![0u32; n];
    let mut sizes = Vec::new();
    let mut current_label = 0u32;

    for start in 0..n {
        if !mask[start] || labels[start] != 0 {
            continue;
        }
        current_label += 1;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        labels[start] = current_label;
        let mut size = 0usize;

        while let Some(idx) = queue.pop_front() {
            size += 1;
            let r = (idx / width) as i32;
            let c = (idx % width) as i32;

            let neighbors: &[(i32, i32)] = match connectivity {
                Connectivity::Four => &[(-1, 0), (1, 0), (0, -1), (0, 1)],
                Connectivity::Eight => &[
                    (-1, -1), (-1, 0), (-1, 1),
                    (0, -1),           (0, 1),
                    (1, -1),  (1, 0),  (1, 1),
                ],
            };

            for &(dr, dc) in neighbors {
                let nr = r + dr;
                let nc = c + dc;
                if nr < 0 || nr >= height as i32 || nc < 0 || nc >= width as i32 {
                    continue;
                }
                let ni = nr as usize * width + nc as usize;
                if mask[ni] && labels[ni] == 0 {
                    labels[ni] = current_label;
                    queue.push_back(ni);
                }
            }
        }
        sizes.push(size);
    }

    Components { labels, sizes }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_mask_no_components() {
        let mask = vec![false; 9];
        let c = label_components(&mask, 3, 3, Connectivity::Four);
        assert!(c.sizes.is_empty());
        assert!(c.labels.iter().all(|&l| l == 0));
    }

    #[test]
    fn single_pixel_component() {
        let mut mask = vec![false; 9];
        mask[4] = true; // centre of 3×3
        let c = label_components(&mask, 3, 3, Connectivity::Four);
        assert_eq!(c.sizes.len(), 1);
        assert_eq!(c.sizes[0], 1);
        assert_eq!(c.labels[4], 1);
    }

    #[test]
    fn two_separate_components_four() {
        // 3×3: corners set, 4-connectivity so they are separate.
        // Pattern:
        //   1 0 1
        //   0 0 0
        //   1 0 1
        let mut mask = vec![false; 9];
        mask[0] = true; mask[2] = true;
        mask[6] = true; mask[8] = true;
        let c = label_components(&mask, 3, 3, Connectivity::Four);
        assert_eq!(c.sizes.len(), 4, "four separate corner pixels");
        for &s in &c.sizes {
            assert_eq!(s, 1);
        }
    }

    #[test]
    fn two_components_become_one_eight() {
        // Two pixels at (0,0) and (1,1) — diagonally adjacent.
        // 4-connectivity: separate; 8-connectivity: connected.
        let mut mask = vec![false; 4]; // 2×2 grid
        mask[0] = true; // (0,0)
        mask[3] = true; // (1,1)
        let c4 = label_components(&mask, 2, 2, Connectivity::Four);
        assert_eq!(c4.sizes.len(), 2, "4-conn: two separate components");
        let c8 = label_components(&mask, 2, 2, Connectivity::Eight);
        assert_eq!(c8.sizes.len(), 1, "8-conn: one connected component");
        assert_eq!(c8.sizes[0], 2);
    }

    #[test]
    fn horizontal_strip_four_connected() {
        // 1×5 mask: all true → one component of size 5.
        let mask = vec![true; 5];
        let c = label_components(&mask, 5, 1, Connectivity::Four);
        assert_eq!(c.sizes.len(), 1);
        assert_eq!(c.sizes[0], 5);
    }

    #[test]
    fn two_separated_rows() {
        // 3×3: rows 0 and 2 all true, row 1 all false.
        let mut mask = vec![false; 9];
        for c in 0..3 { mask[c] = true; mask[6 + c] = true; }
        let comp = label_components(&mask, 3, 3, Connectivity::Four);
        assert_eq!(comp.sizes.len(), 2, "two rows = two components");
        assert_eq!(comp.sizes[0], 3);
        assert_eq!(comp.sizes[1], 3);
    }

    #[test]
    fn sizes_sum_equals_foreground_count() {
        let mut mask = vec![false; 16];
        // L-shape: rows 0-2 col 0, row 2 col 0-3.
        for r in 0..3 { mask[r * 4] = true; }
        for c in 0..4 { mask[2 * 4 + c] = true; }
        let comp = label_components(&mask, 4, 4, Connectivity::Four);
        let total: usize = comp.sizes.iter().sum();
        let fg: usize = mask.iter().filter(|&&b| b).count();
        assert_eq!(total, fg);
    }
}
