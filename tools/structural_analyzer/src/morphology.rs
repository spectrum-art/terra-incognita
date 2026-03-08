//! Morphological operations on boolean pixel masks.
//!
//! Provides dilation, erosion, closing (dilation then erosion), and
//! Zhang-Suen iterative thinning (skeletonization). These operate on
//! flat row-major boolean arrays of size `width × height`.

// ── Structuring element ────────────────────────────────────────────────────

/// Precompute the (dr, dc) offsets within a circular structuring element
/// of the given radius. Includes the origin (0, 0).
fn structuring_element(radius: i32) -> Vec<(i32, i32)> {
    let mut offsets = Vec::new();
    for dr in -radius..=radius {
        for dc in -radius..=radius {
            if dr * dr + dc * dc <= radius * radius {
                offsets.push((dr, dc));
            }
        }
    }
    offsets
}

// ── Dilation ───────────────────────────────────────────────────────────────

/// Morphological dilation: for each ON pixel in `mask`, set all pixels within
/// the circular structuring element of the given `radius` as ON in the output.
pub fn dilate(mask: &[bool], width: usize, height: usize, radius: i32) -> Vec<bool> {
    let n = width * height;
    let mut output = vec![false; n];
    let se = structuring_element(radius);

    for (idx, &on) in mask.iter().enumerate().take(n) {
        if !on {
            continue;
        }
        let r = (idx / width) as i32;
        let c = (idx % width) as i32;
        for &(dr, dc) in &se {
            let nr = r + dr;
            let nc = c + dc;
            if nr >= 0 && nr < height as i32 && nc >= 0 && nc < width as i32 {
                output[nr as usize * width + nc as usize] = true;
            }
        }
    }

    output
}

// ── Erosion ────────────────────────────────────────────────────────────────

/// Morphological erosion: start with a copy of `mask`. For each OFF pixel in
/// `mask`, set all pixels within the circular structuring element of the given
/// `radius` as OFF in the output.
pub fn erode(mask: &[bool], width: usize, height: usize, radius: i32) -> Vec<bool> {
    let n = width * height;
    let mut output = mask.to_vec();
    let se = structuring_element(radius);

    for (idx, &on) in mask.iter().enumerate().take(n) {
        if on {
            continue; // only propagate from OFF pixels
        }
        let r = (idx / width) as i32;
        let c = (idx % width) as i32;
        for &(dr, dc) in &se {
            let nr = r + dr;
            let nc = c + dc;
            if nr >= 0 && nr < height as i32 && nc >= 0 && nc < width as i32 {
                output[nr as usize * width + nc as usize] = false;
            }
        }
    }

    output
}

// ── Closing ────────────────────────────────────────────────────────────────

/// Morphological closing: dilation followed by erosion with the same radius.
/// Fills gaps smaller than `radius` between ON pixels.
pub fn close(mask: &[bool], width: usize, height: usize, radius: i32) -> Vec<bool> {
    let dilated = dilate(mask, width, height, radius);
    erode(&dilated, width, height, radius)
}

// ── Skeletonization (Zhang-Suen) ────────────────────────────────────────────

/// Iterative thinning using the Zhang-Suen algorithm.
///
/// Neighbors in order: P2=N, P3=NE, P4=E, P5=SE, P6=S, P7=SW, P8=W, P9=NW
/// (row decreases going north; col increases going east).
/// B(P) = number of ON neighbors.
/// A(P) = number of 0→1 transitions in the ordered sequence of 8 neighbors.
pub fn skeletonize(mask: &[bool], width: usize, height: usize) -> Vec<bool> {
    let n = width * height;
    let mut current = mask.to_vec();
    let mut changed = true;

    while changed {
        changed = false;

        // Sub-iteration 1: remove if B in [2,6], A==1, P2*P4*P6==0, P4*P6*P8==0
        let to_remove1 = sub_iteration(&current, width, height, false);
        for idx in 0..n {
            if to_remove1[idx] && current[idx] {
                current[idx] = false;
                changed = true;
            }
        }

        // Sub-iteration 2: remove if B in [2,6], A==1, P2*P4*P8==0, P2*P6*P8==0
        let to_remove2 = sub_iteration(&current, width, height, true);
        for idx in 0..n {
            if to_remove2[idx] && current[idx] {
                current[idx] = false;
                changed = true;
            }
        }
    }

    current
}

/// Compute a removal mask for one Zhang-Suen sub-iteration.
/// `second_iter` selects between the two sets of conditions.
fn sub_iteration(mask: &[bool], width: usize, height: usize, second_iter: bool) -> Vec<bool> {
    let n = width * height;
    let mut to_remove = vec![false; n];

    for r in 1..height as i32 - 1 {
        for c in 1..width as i32 - 1 {
            let idx = r as usize * width + c as usize;
            if !mask[idx] {
                continue;
            }

            // Neighbors in Zhang-Suen order: P2..P9
            // P2=N, P3=NE, P4=E, P5=SE, P6=S, P7=SW, P8=W, P9=NW
            let neighbors = [
                mask[(r - 1) as usize * width + c as usize],       // P2 = N
                mask[(r - 1) as usize * width + (c + 1) as usize], // P3 = NE
                mask[r as usize * width + (c + 1) as usize],       // P4 = E
                mask[(r + 1) as usize * width + (c + 1) as usize], // P5 = SE
                mask[(r + 1) as usize * width + c as usize],       // P6 = S
                mask[(r + 1) as usize * width + (c - 1) as usize], // P7 = SW
                mask[r as usize * width + (c - 1) as usize],       // P8 = W
                mask[(r - 1) as usize * width + (c - 1) as usize], // P9 = NW
            ];

            // B(P): count of ON neighbors
            let b: u32 = neighbors.iter().filter(|&&v| v).count() as u32;
            if !(2..=6).contains(&b) {
                continue;
            }

            // A(P): count of 0→1 transitions in cyclic order P2..P9, P2
            let mut a = 0u32;
            for i in 0..8 {
                let curr = neighbors[i];
                let next = neighbors[(i + 1) % 8];
                if !curr && next {
                    a += 1;
                }
            }
            if a != 1 {
                continue;
            }

            let p2 = neighbors[0];
            let p4 = neighbors[2];
            let p6 = neighbors[4];
            let p8 = neighbors[6];

            let remove = if !second_iter {
                // Sub-iteration 1: not (p2 AND p4 AND p6) AND not (p4 AND p6 AND p8)
                !p4 || !p6 || !p2 && !p8
            } else {
                // Sub-iteration 2: not (p2 AND p4 AND p8) AND not (p2 AND p6 AND p8)
                !p2 || !p8 || !p4 && !p6
            };

            if remove {
                to_remove[idx] = true;
            }
        }
    }

    to_remove
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Dilation tests ─────────────────────────────────────────────────────────

    #[test]
    fn dilate_single_pixel_expands() {
        // Single ON pixel in a 9×9 grid at (4,4). Dilation with R=2 should
        // produce many ON pixels.
        let w = 9usize;
        let h = 9usize;
        let mut mask = vec![false; w * h];
        mask[4 * w + 4] = true;
        let out = dilate(&mask, w, h, 2);
        let count = out.iter().filter(|&&v| v).count();
        assert!(count > 1, "expected dilation to expand single pixel, got {}", count);
        // With R=2, SE has offsets dr²+dc²<=4 → at least a 3×3 cross and more.
        assert!(count >= 9, "expected at least 9 ON pixels after R=2 dilation, got {}", count);
    }

    #[test]
    fn dilate_preserves_original_pixels() {
        let w = 9usize;
        let h = 9usize;
        let mut mask = vec![false; w * h];
        mask[4 * w + 4] = true;
        let out = dilate(&mask, w, h, 3);
        // Original ON pixel must still be ON.
        assert!(out[4 * w + 4], "original pixel must remain ON after dilation");
    }

    // ── Erosion tests ──────────────────────────────────────────────────────────

    #[test]
    fn erode_isolated_pixel_removed() {
        // Single ON pixel in 9×9 at (4,4). Erosion with R=2 should turn it OFF.
        let w = 9usize;
        let h = 9usize;
        let mut mask = vec![false; w * h];
        mask[4 * w + 4] = true;
        let out = erode(&mask, w, h, 2);
        let count = out.iter().filter(|&&v| v).count();
        assert_eq!(count, 0, "isolated pixel should be eroded away");
    }

    #[test]
    fn erode_large_solid_block_shrinks() {
        // Solid 11×11 block in a 15×15 grid. After erosion with R=2, it shrinks
        // inward by 2 pixels on each side.
        let w = 15usize;
        let h = 15usize;
        let mut mask = vec![false; w * h];
        for r in 2..13 {
            for c in 2..13 {
                mask[r * w + c] = true;
            }
        }
        let before = mask.iter().filter(|&&v| v).count();
        let out = erode(&mask, w, h, 2);
        let after = out.iter().filter(|&&v| v).count();
        assert!(after < before, "erosion should shrink the block");
    }

    // ── Closing tests ──────────────────────────────────────────────────────────

    #[test]
    fn close_fills_small_gap() {
        // Two ON pixels separated by one OFF pixel: [T, F, T] in a 1×5 grid.
        // Closing with R=1 should fill the gap (the OFF pixel connects to both ON pixels
        // within radius 1).
        let w = 5usize;
        let h = 1usize;
        let mask = vec![true, false, true, false, false];
        let out = close(&mask, w, h, 1);
        // After dilation: T,T,T,T,F — after erosion the gap at index 1 should be filled.
        assert!(out[1], "closing with R=1 should fill single-pixel gap at index 1");
    }

    #[test]
    fn close_preserves_isolated_large_feature() {
        // A large solid block: closing should not dramatically change it.
        let w = 20usize;
        let h = 20usize;
        let mut mask = vec![false; w * h];
        for r in 5..15 {
            for c in 5..15 {
                mask[r * w + c] = true;
            }
        }
        let out = close(&mask, w, h, 2);
        // The center of the block should still be ON.
        assert!(out[10 * w + 10], "center of large block must stay ON after closing");
    }

    // ── Skeletonization tests ──────────────────────────────────────────────────

    #[test]
    fn skeletonize_solid_rectangle_becomes_line() {
        // A 9-row-wide solid block embedded in a larger grid of background pixels.
        // The block spans rows 2..11, cols 2..22 in a 15×25 grid.
        // Zhang-Suen thinning requires background pixels adjacent to the block.
        let w = 25usize;
        let h = 15usize;
        let mut mask = vec![false; w * h];
        let block_r_start = 2usize;
        let block_r_end = 12usize; // exclusive
        let block_c_start = 2usize;
        let block_c_end = 23usize; // exclusive → 21 wide, 10 tall
        let block_pixels = (block_r_end - block_r_start) * (block_c_end - block_c_start);
        for r in block_r_start..block_r_end {
            for c in block_c_start..block_c_end {
                mask[r * w + c] = true;
            }
        }
        let out = skeletonize(&mask, w, h);
        let after = out.iter().filter(|&&v| v).count();
        // The skeleton should be significantly thinner than the original block.
        assert!(
            after < block_pixels / 2,
            "skeleton of solid rectangle should be much smaller than original; block_pixels={}, after={}",
            block_pixels,
            after
        );
    }

    #[test]
    fn skeletonize_single_line_unchanged() {
        // A 1-pixel-wide horizontal line through a 7×20 grid.
        // After skeletonization, most pixels in the line should survive.
        let w = 20usize;
        let h = 7usize;
        let mut mask = vec![false; w * h];
        let row = 3;
        for c in 0..w {
            mask[row * w + c] = true;
        }
        let out = skeletonize(&mask, w, h);
        let surviving: usize = (0..w).filter(|&c| out[row * w + c]).count();
        // At least half the interior pixels should survive.
        assert!(
            surviving > w / 2,
            "single-width line should mostly survive skeletonization; surviving={}/{}",
            surviving,
            w
        );
    }

    #[test]
    fn skeletonize_empty_mask_stays_empty() {
        let mask = vec![false; 9];
        let out = skeletonize(&mask, 3, 3);
        assert!(out.iter().all(|&v| !v), "empty mask must stay empty after skeletonize");
    }
}
