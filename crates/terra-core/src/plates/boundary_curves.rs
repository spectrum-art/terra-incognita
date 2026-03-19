//! Boundary polyline extraction from raster plate-boundary fields.

use std::collections::{HashMap, HashSet, VecDeque};
use std::f64::consts::PI;

use crate::plates::plate_dynamics::BoundaryCharacter;
use crate::sphere::{great_circle_distance_rad, Vec3};

const EARTH_RADIUS_KM: f64 = 6371.0;
const MIN_COMPONENT_PIXELS: usize = 5;
const CHAIKIN_ITERS: usize = 2;

/// A single vertex along a boundary polyline, with interpolated metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryVertex {
    /// Position in grid coordinates (fractional — sub-pixel precision after smoothing)
    pub x: f64,
    pub y: f64,
    /// Position in spherical coordinates (radians)
    pub lat: f64,
    pub lon: f64,
    /// Convergent rate at this vertex (cm/yr, positive = convergent)
    pub convergent_rate: f32,
    /// Transform rate at this vertex (cm/yr)
    pub transform_rate: f32,
    /// Unit tangent vector along the polyline at this vertex (east, north components)
    pub tangent: (f32, f32),
    /// Unit normal vector perpendicular to the polyline (east, north components)
    /// Points toward the overriding plate side for convergent boundaries.
    pub normal: (f32, f32),
}

/// A single connected boundary segment extracted as a smooth ordered polyline.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryPolyline {
    /// The two plate IDs this boundary separates.
    pub plate_a: u32,
    pub plate_b: u32,
    /// Dominant boundary type along the segment.
    pub dominant_character: BoundaryType,
    /// Ordered vertices from one end to the other or around a closed loop.
    pub vertices: Vec<BoundaryVertex>,
    /// Whether the polyline forms a closed loop.
    pub is_closed: bool,
    /// Cumulative arc-length in kilometres at each vertex.
    pub arc_lengths: Vec<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryType {
    Convergent,
    Divergent,
    Transform,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct RasterBoundaryComponent {
    plate_a: u8,
    plate_b: u8,
    pixels: Vec<usize>,
}

#[derive(Clone, Copy, Debug)]
struct ControlVertex {
    x: f64,
    y: f64,
    convergent_rate: f32,
    transform_rate: f32,
    reference_normal: (f32, f32),
}

/// Extract smooth ordered boundary polylines from raster plate-boundary data.
pub fn extract_boundary_polylines(
    boundary_field: &[BoundaryCharacter],
    is_boundary: &[bool],
    plate_ids: &[u8],
    width: usize,
    height: usize,
) -> Vec<BoundaryPolyline> {
    let components =
        extract_components_from_mask(boundary_field, is_boundary, plate_ids, width, height);

    let mut polylines = Vec::new();
    for component in components
        .into_iter()
        .filter(|component| component.pixels.len() >= MIN_COMPONENT_PIXELS)
    {
        let skeleton = thin_component(&component, width, height);
        if skeleton.is_empty() {
            continue;
        }
        let (ordered_pixels, is_closed) = order_skeleton_pixels(&skeleton, width, height);
        if ordered_pixels.is_empty() {
            continue;
        }
        let dominant_character = dominant_character(boundary_field, &component.pixels);
        let control_vertices =
            ordered_control_vertices(&ordered_pixels, boundary_field, plate_ids, width, height);
        let smoothed = chaikin_smooth(&control_vertices, is_closed, CHAIKIN_ITERS);
        let (vertices, arc_lengths) = build_polyline_vertices(&smoothed, is_closed, width, height);

        polylines.push(BoundaryPolyline {
            plate_a: u32::from(component.plate_a),
            plate_b: u32::from(component.plate_b),
            dominant_character,
            vertices,
            is_closed,
            arc_lengths,
        });
    }

    if cfg!(debug_assertions) && std::env::var_os("TERRA_DEBUG_BOUNDARY_CURVES").is_some() {
        eprintln!("boundary polylines extracted: {}", polylines.len());
        for polyline in &polylines {
            let mean_convergent = if polyline.vertices.is_empty() {
                0.0
            } else {
                polyline
                    .vertices
                    .iter()
                    .map(|vertex| vertex.convergent_rate)
                    .sum::<f32>()
                    / polyline.vertices.len() as f32
            };
            let total_arc_km = polyline.arc_lengths.last().copied().unwrap_or(0.0);
            eprintln!(
                "  pair {}-{} {:?}: {} verts, {:.1} km, mean conv {:.2}",
                polyline.plate_a,
                polyline.plate_b,
                polyline.dominant_character,
                polyline.vertices.len(),
                total_arc_km,
                mean_convergent
            );
        }
    }

    polylines
}

fn extract_components_from_mask(
    boundary_field: &[BoundaryCharacter],
    mask: &[bool],
    plate_ids: &[u8],
    width: usize,
    height: usize,
) -> Vec<RasterBoundaryComponent> {
    let mut visited = vec![false; mask.len()];
    let mut components = Vec::new();

    for start in 0..mask.len() {
        if visited[start] || !mask[start] {
            continue;
        }
        let pair = ordered_pair(plate_ids[start], boundary_field[start].neighbor_plate);
        let mut queue = VecDeque::from([start]);
        let mut pixels = Vec::new();
        visited[start] = true;

        while let Some(idx) = queue.pop_front() {
            pixels.push(idx);
            for neighbor in neighbors8(idx, width, height) {
                if visited[neighbor] || !mask[neighbor] {
                    continue;
                }
                if ordered_pair(plate_ids[neighbor], boundary_field[neighbor].neighbor_plate)
                    != pair
                {
                    continue;
                }
                visited[neighbor] = true;
                queue.push_back(neighbor);
            }
        }

        pixels.sort_unstable();
        components.push(RasterBoundaryComponent {
            plate_a: pair.0,
            plate_b: pair.1,
            pixels,
        });
    }

    components.sort_by_key(|component| (component.plate_a, component.plate_b, component.pixels[0]));
    components
}

fn thin_component(component: &RasterBoundaryComponent, width: usize, height: usize) -> Vec<usize> {
    let local = component_local_grid(&component.pixels, width, height);
    let skeleton = zhang_suen_thin(local.mask, local.mask_width, local.mask_height);
    let mut result = Vec::new();

    for local_y in 0..local.mask_height {
        for local_x in 0..local.mask_width {
            let mask_idx = local_y * local.mask_width + local_x;
            if !skeleton[mask_idx] {
                continue;
            }
            let x_unwrapped = local.min_x + local_x as isize - 1;
            let y = local.min_y + local_y as isize - 1;
            if !(0..height as isize).contains(&y) {
                continue;
            }
            let x = x_unwrapped.rem_euclid(width as isize) as usize;
            result.push(y as usize * width + x);
        }
    }

    result.sort_unstable();
    result.dedup();
    simplify_skeleton_to_path(&result, width, height)
}

fn order_skeleton_pixels(skeleton: &[usize], width: usize, height: usize) -> (Vec<usize>, bool) {
    let index_by_pixel = skeleton
        .iter()
        .copied()
        .enumerate()
        .map(|(i, idx)| (idx, i))
        .collect::<HashMap<_, _>>();
    let mut adjacency = vec![Vec::new(); skeleton.len()];

    for (i, &idx) in skeleton.iter().enumerate() {
        for neighbor in neighbors8(idx, width, height) {
            if let Some(&neighbor_idx) = index_by_pixel.get(&neighbor) {
                adjacency[i].push(neighbor_idx);
            }
        }
        adjacency[i].sort_unstable();
    }

    let start = adjacency
        .iter()
        .enumerate()
        .find_map(|(i, neighbors)| (neighbors.len() <= 1).then_some(i))
        .unwrap_or(0);

    let mut ordered = Vec::with_capacity(skeleton.len());
    let mut visited = vec![false; skeleton.len()];
    let mut previous = None;
    let mut current = start;

    loop {
        ordered.push(skeleton[current]);
        visited[current] = true;

        let mut candidates = adjacency[current]
            .iter()
            .copied()
            .filter(|&neighbor| !visited[neighbor])
            .collect::<Vec<_>>();
        if candidates.is_empty() {
            break;
        }
        candidates.sort_by(|&a, &b| {
            let score_a = branch_score(previous, current, a, skeleton, width);
            let score_b = branch_score(previous, current, b, skeleton, width);
            score_b
                .partial_cmp(&score_a)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| skeleton[a].cmp(&skeleton[b]))
        });
        previous = Some(current);
        current = candidates[0];
    }

    let is_closed = ordered.len() == skeleton.len()
        && adjacency[current]
            .iter()
            .any(|&neighbor| skeleton[neighbor] == ordered[0]);

    (ordered, is_closed)
}

fn branch_score(
    previous: Option<usize>,
    current: usize,
    candidate: usize,
    skeleton: &[usize],
    width: usize,
) -> f64 {
    let current_pos = idx_to_xy(skeleton[current], width);
    let candidate_pos = idx_to_xy(skeleton[candidate], width);
    let forward = wrapped_grid_delta(current_pos, candidate_pos, width);
    if let Some(previous) = previous {
        let previous_pos = idx_to_xy(skeleton[previous], width);
        let incoming = wrapped_grid_delta(previous_pos, current_pos, width);
        let incoming_norm = (incoming.0 * incoming.0 + incoming.1 * incoming.1).sqrt();
        let forward_norm = (forward.0 * forward.0 + forward.1 * forward.1).sqrt();
        if incoming_norm <= f64::EPSILON || forward_norm <= f64::EPSILON {
            return -1.0;
        }
        (incoming.0 * forward.0 + incoming.1 * forward.1) / (incoming_norm * forward_norm)
    } else {
        0.0
    }
}

fn ordered_control_vertices(
    ordered_pixels: &[usize],
    boundary_field: &[BoundaryCharacter],
    plate_ids: &[u8],
    width: usize,
    height: usize,
) -> Vec<ControlVertex> {
    let mut vertices = Vec::with_capacity(ordered_pixels.len());
    let mut x = 0.0_f64;

    for (i, &idx) in ordered_pixels.iter().enumerate() {
        let (col, row) = idx_to_xy(idx, width);
        if i == 0 {
            x = col as f64;
        } else {
            let (prev_col, _) = idx_to_xy(ordered_pixels[i - 1], width);
            let delta = wrapped_grid_delta((prev_col, row), (col, row), width);
            x += delta.0;
        }

        let character = boundary_field[idx];
        let own_plate = plate_ids[idx];
        let mut reference_normal = (character.normal_east, character.normal_north);
        if character.overriding_plate == own_plate {
            reference_normal = (-reference_normal.0, -reference_normal.1);
        }
        reference_normal = normalize_2d(reference_normal);

        vertices.push(ControlVertex {
            x,
            y: row as f64,
            convergent_rate: character.convergent_rate,
            transform_rate: character.transform_rate,
            reference_normal,
        });
    }

    let _ = height;
    vertices
}

fn chaikin_smooth(
    points: &[ControlVertex],
    is_closed: bool,
    iterations: usize,
) -> Vec<ControlVertex> {
    if points.len() < 2 {
        return points.to_vec();
    }

    let mut current = points.to_vec();
    for _ in 0..iterations {
        let n = current.len();
        let mut next = Vec::with_capacity(n * 2);

        if !is_closed {
            next.push(current[0]);
        }

        let segments = if is_closed { n } else { n - 1 };
        for i in 0..segments {
            let a = current[i];
            let b = current[(i + 1) % n];
            next.push(lerp_control(a, b, 0.25));
            next.push(lerp_control(a, b, 0.75));
        }

        if !is_closed {
            next.push(current[n - 1]);
        }
        current = next;
    }

    current
}

fn build_polyline_vertices(
    controls: &[ControlVertex],
    is_closed: bool,
    width: usize,
    height: usize,
) -> (Vec<BoundaryVertex>, Vec<f64>) {
    let points = controls
        .iter()
        .map(|control| control_to_point(*control, width, height))
        .collect::<Vec<_>>();

    let mut vertices = Vec::with_capacity(controls.len());
    for (i, &control) in controls.iter().enumerate() {
        let point = points[i];
        let tangent = tangent_at(i, &points, is_closed);
        let mut normal = (-tangent.1, tangent.0);
        if dot_2d(normal, control.reference_normal) < 0.0 {
            normal = (-normal.0, -normal.1);
        }
        let wrapped_x = control.x.rem_euclid(width as f64);
        let lat = PI * 0.5 - (control.y + 0.5) * PI / height as f64;
        let lon = -PI + (wrapped_x + 0.5) * (2.0 * PI / width as f64);
        let _ = point;
        vertices.push(BoundaryVertex {
            x: wrapped_x,
            y: control.y,
            lat,
            lon,
            convergent_rate: control.convergent_rate,
            transform_rate: control.transform_rate,
            tangent,
            normal,
        });
    }

    let mut arc_lengths = Vec::with_capacity(vertices.len());
    let mut cumulative = 0.0_f64;
    for i in 0..vertices.len() {
        if i > 0 {
            cumulative += great_circle_distance_rad(points[i - 1], points[i]) * EARTH_RADIUS_KM;
        }
        arc_lengths.push(cumulative);
    }

    (vertices, arc_lengths)
}

fn dominant_character(boundary_field: &[BoundaryCharacter], pixels: &[usize]) -> BoundaryType {
    let mean_convergent = pixels
        .iter()
        .map(|&idx| boundary_field[idx].convergent_rate)
        .sum::<f32>()
        / pixels.len().max(1) as f32;

    if mean_convergent > 1.0 {
        BoundaryType::Convergent
    } else if mean_convergent < -1.0 {
        BoundaryType::Divergent
    } else {
        BoundaryType::Transform
    }
}

fn control_to_point(control: ControlVertex, width: usize, height: usize) -> Vec3 {
    let wrapped_x = control.x.rem_euclid(width as f64);
    let lat_deg = 90.0 - (control.y + 0.5) * 180.0 / height as f64;
    let lon_deg = -180.0 + (wrapped_x + 0.5) * 360.0 / width as f64;
    Vec3::from_latlon(lat_deg, lon_deg)
}

fn tangent_at(index: usize, points: &[Vec3], is_closed: bool) -> (f32, f32) {
    let point = points[index];
    let (east, north) = local_east_north(point);

    let direction = if points.len() == 1 {
        (1.0, 0.0)
    } else if index == 0 && !is_closed {
        project_to_tangent(point, points[index + 1], east, north)
    } else if index + 1 == points.len() && !is_closed {
        let prev = project_to_tangent(point, points[index - 1], east, north);
        (-prev.0, -prev.1)
    } else {
        let prev = if index == 0 {
            points[points.len() - 1]
        } else {
            points[index - 1]
        };
        let next = if index + 1 == points.len() {
            points[0]
        } else {
            points[index + 1]
        };
        let prev_vec = project_to_tangent(point, prev, east, north);
        let next_vec = project_to_tangent(point, next, east, north);
        (next_vec.0 - prev_vec.0, next_vec.1 - prev_vec.1)
    };

    normalize_2d((direction.0 as f32, direction.1 as f32))
}

fn project_to_tangent(point: Vec3, sample: Vec3, east: Vec3, north: Vec3) -> (f64, f64) {
    let tangent = Vec3::new(
        sample.x - point.x * point.dot(sample),
        sample.y - point.y * point.dot(sample),
        sample.z - point.z * point.dot(sample),
    );
    (tangent.dot(east), tangent.dot(north))
}

fn local_east_north(point: Vec3) -> (Vec3, Vec3) {
    let east_raw = if point.x.abs() + point.y.abs() > 1e-12 {
        Vec3::new(-point.y, point.x, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let east = east_raw.normalize();
    let north = point.cross(east).normalize();
    (east, north)
}

fn lerp_control(a: ControlVertex, b: ControlVertex, t: f64) -> ControlVertex {
    let mut reference_normal = (
        ((1.0 - t) as f32) * a.reference_normal.0 + (t as f32) * b.reference_normal.0,
        ((1.0 - t) as f32) * a.reference_normal.1 + (t as f32) * b.reference_normal.1,
    );
    reference_normal = normalize_2d(reference_normal);
    ControlVertex {
        x: a.x + (b.x - a.x) * t,
        y: a.y + (b.y - a.y) * t,
        convergent_rate: a.convergent_rate + (b.convergent_rate - a.convergent_rate) * t as f32,
        transform_rate: a.transform_rate + (b.transform_rate - a.transform_rate) * t as f32,
        reference_normal,
    }
}

fn component_local_grid(pixels: &[usize], width: usize, height: usize) -> LocalComponentGrid {
    let mut local_xy = HashMap::<usize, (isize, isize)>::new();
    let mut queue = VecDeque::from([pixels[0]]);
    let (start_x, start_y) = idx_to_xy(pixels[0], width);
    local_xy.insert(pixels[0], (start_x as isize, start_y as isize));

    let pixel_set = pixels.iter().copied().collect::<HashSet<_>>();
    while let Some(idx) = queue.pop_front() {
        let &(x, y) = local_xy.get(&idx).expect("seeded local coordinate");
        for neighbor in neighbors8(idx, width, height) {
            if !pixel_set.contains(&neighbor) || local_xy.contains_key(&neighbor) {
                continue;
            }
            let (neighbor_col, neighbor_row) = idx_to_xy(neighbor, width);
            let (current_col, current_row) = idx_to_xy(idx, width);
            let delta = wrapped_grid_delta(
                (current_col, current_row),
                (neighbor_col, neighbor_row),
                width,
            );
            local_xy.insert(neighbor, (x + delta.0 as isize, y + delta.1 as isize));
            queue.push_back(neighbor);
        }
    }

    let min_x = local_xy.values().map(|&(x, _)| x).min().unwrap_or(0);
    let max_x = local_xy.values().map(|&(x, _)| x).max().unwrap_or(0);
    let min_y = local_xy.values().map(|&(_, y)| y).min().unwrap_or(0);
    let max_y = local_xy.values().map(|&(_, y)| y).max().unwrap_or(0);
    let mask_width = (max_x - min_x + 3) as usize;
    let mask_height = (max_y - min_y + 3) as usize;
    let mut mask = vec![false; mask_width * mask_height];

    for &(x, y) in local_xy.values() {
        let local_x = (x - min_x + 1) as usize;
        let local_y = (y - min_y + 1) as usize;
        mask[local_y * mask_width + local_x] = true;
    }

    LocalComponentGrid {
        min_x,
        min_y,
        mask_width,
        mask_height,
        mask,
    }
}

struct LocalComponentGrid {
    min_x: isize,
    min_y: isize,
    mask_width: usize,
    mask_height: usize,
    mask: Vec<bool>,
}

fn zhang_suen_thin(mut mask: Vec<bool>, width: usize, height: usize) -> Vec<bool> {
    let neighbors = |idx: usize, mask: &[bool]| -> [bool; 8] {
        let row = idx / width;
        let col = idx % width;
        let p2 = row > 0 && mask[(row - 1) * width + col];
        let p3 = row > 0 && col + 1 < width && mask[(row - 1) * width + (col + 1)];
        let p4 = col + 1 < width && mask[row * width + (col + 1)];
        let p5 = row + 1 < height && col + 1 < width && mask[(row + 1) * width + (col + 1)];
        let p6 = row + 1 < height && mask[(row + 1) * width + col];
        let p7 = row + 1 < height && col > 0 && mask[(row + 1) * width + (col - 1)];
        let p8 = col > 0 && mask[row * width + (col - 1)];
        let p9 = row > 0 && col > 0 && mask[(row - 1) * width + (col - 1)];
        [p2, p3, p4, p5, p6, p7, p8, p9]
    };

    loop {
        let mut changed = false;
        for step in 0..2 {
            let mut to_remove = Vec::new();
            for idx in 0..mask.len() {
                if !mask[idx] {
                    continue;
                }
                let row = idx / width;
                let col = idx % width;
                if row == 0 || row + 1 == height || col == 0 || col + 1 == width {
                    continue;
                }
                let ns = neighbors(idx, &mask);
                let count = ns.iter().filter(|&&v| v).count();
                if !(2..=6).contains(&count) {
                    continue;
                }
                let transitions = ns
                    .iter()
                    .copied()
                    .zip(ns.iter().copied().cycle().skip(1))
                    .take(8)
                    .filter(|&(a, b)| !a && b)
                    .count();
                if transitions != 1 {
                    continue;
                }
                let (p2, _p3, p4, _p5, p6, _p7, p8, _p9) =
                    (ns[0], ns[1], ns[2], ns[3], ns[4], ns[5], ns[6], ns[7]);
                let passes = if step == 0 {
                    !p4 || !p6 || (!p2 && !p8)
                } else {
                    !p2 || !p8 || (!p4 && !p6)
                };
                if passes {
                    to_remove.push(idx);
                }
            }
            if !to_remove.is_empty() {
                changed = true;
                for idx in to_remove {
                    mask[idx] = false;
                }
            }
        }
        if !changed {
            break;
        }
    }

    mask
}

fn simplify_skeleton_to_path(skeleton: &[usize], width: usize, height: usize) -> Vec<usize> {
    if skeleton.len() <= 2 {
        return skeleton.to_vec();
    }

    let index_by_pixel = skeleton
        .iter()
        .copied()
        .enumerate()
        .map(|(i, idx)| (idx, i))
        .collect::<HashMap<_, _>>();
    let mut adjacency = vec![Vec::new(); skeleton.len()];
    for (i, &idx) in skeleton.iter().enumerate() {
        for neighbor in neighbors8(idx, width, height) {
            if let Some(&neighbor_idx) = index_by_pixel.get(&neighbor) {
                adjacency[i].push(neighbor_idx);
            }
        }
        adjacency[i].sort_unstable();
        adjacency[i].dedup();
    }

    let mut visited = vec![false; skeleton.len()];
    let mut best_path = Vec::new();

    for start in 0..skeleton.len() {
        if visited[start] {
            continue;
        }
        let mut queue = VecDeque::from([start]);
        let mut component = Vec::new();
        visited[start] = true;
        while let Some(node) = queue.pop_front() {
            component.push(node);
            for &neighbor in &adjacency[node] {
                if visited[neighbor] {
                    continue;
                }
                visited[neighbor] = true;
                queue.push_back(neighbor);
            }
        }

        let mut endpoints = component
            .iter()
            .copied()
            .filter(|&node| adjacency[node].len() <= 1)
            .collect::<Vec<_>>();
        if endpoints.is_empty() {
            endpoints = component.clone();
        }

        let mut component_best = Vec::new();
        for &endpoint in &endpoints {
            let (path, _) = farthest_path_from(endpoint, &component, &adjacency);
            if path.len() > component_best.len() {
                component_best = path;
            }
        }

        if component_best.len() > best_path.len() {
            best_path = component_best;
        }
    }

    let mut pixels = best_path
        .into_iter()
        .map(|node| skeleton[node])
        .collect::<Vec<_>>();
    pixels.sort_unstable();
    pixels
}

fn farthest_path_from(
    start: usize,
    component: &[usize],
    adjacency: &[Vec<usize>],
) -> (Vec<usize>, usize) {
    let mut queue = VecDeque::from([start]);
    let mut distance = HashMap::from([(start, 0usize)]);
    let mut parent = HashMap::new();

    while let Some(node) = queue.pop_front() {
        let node_distance = distance[&node];
        for &neighbor in &adjacency[node] {
            if !component.contains(&neighbor) || distance.contains_key(&neighbor) {
                continue;
            }
            distance.insert(neighbor, node_distance + 1);
            parent.insert(neighbor, node);
            queue.push_back(neighbor);
        }
    }

    let mut farthest = start;
    for &node in component {
        if distance.get(&node).copied().unwrap_or(0) > distance.get(&farthest).copied().unwrap_or(0)
        {
            farthest = node;
        }
    }

    let mut path = Vec::new();
    let mut current = farthest;
    path.push(current);
    while let Some(&prev) = parent.get(&current) {
        current = prev;
        path.push(current);
    }
    path.reverse();
    (path, farthest)
}

fn neighbors8(idx: usize, width: usize, height: usize) -> Vec<usize> {
    let row = idx / width;
    let col = idx % width;
    let mut neighbors = Vec::with_capacity(8);
    for dr in -1isize..=1 {
        for dc in -1isize..=1 {
            if dr == 0 && dc == 0 {
                continue;
            }
            let rr = row as isize + dr;
            if rr < 0 || rr >= height as isize {
                continue;
            }
            let cc = ((col as isize + dc).rem_euclid(width as isize)) as usize;
            let neighbor = rr as usize * width + cc;
            if !neighbors.contains(&neighbor) {
                neighbors.push(neighbor);
            }
        }
    }
    neighbors
}

fn idx_to_xy(idx: usize, width: usize) -> (usize, usize) {
    (idx % width, idx / width)
}

fn wrapped_grid_delta(from: (usize, usize), to: (usize, usize), width: usize) -> (f64, f64) {
    let mut dx = to.0 as isize - from.0 as isize;
    let half_width = (width / 2) as isize;
    if dx > half_width {
        dx -= width as isize;
    } else if dx < -half_width {
        dx += width as isize;
    }
    let dy = to.1 as isize - from.1 as isize;
    (dx as f64, dy as f64)
}

fn normalize_2d(vector: (f32, f32)) -> (f32, f32) {
    let length = (vector.0 * vector.0 + vector.1 * vector.1).sqrt();
    if length <= f32::EPSILON {
        (1.0, 0.0)
    } else {
        (vector.0 / length, vector.1 / length)
    }
}

fn dot_2d(a: (f32, f32), b: (f32, f32)) -> f32 {
    a.0 * b.0 + a.1 * b.1
}

fn ordered_pair(a: u8, b: u8) -> (u8, u8) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::OnceLock;

    use crate::plates::plate_dynamics::compute_plate_dynamics;
    use crate::plates::plate_generation::generate_plate_geometry;

    const TEST_WIDTH: usize = 1024;
    const TEST_HEIGHT: usize = 512;
    const TEST_PLATES: usize = 13;
    const TEST_WARP_DEG: f64 = 7.0;

    struct Seed42Fixture {
        boundary_field: Vec<BoundaryCharacter>,
        is_boundary: Vec<bool>,
        plate_ids: Vec<u8>,
        convergent_components_raw: Vec<RasterBoundaryComponent>,
        convergent_components_filtered: Vec<RasterBoundaryComponent>,
        largest_convergent_raw: RasterBoundaryComponent,
        largest_convergent_thinned: Vec<usize>,
        largest_convergent_chain: Vec<usize>,
        largest_polyline: BoundaryPolyline,
    }

    fn seed42_fixture() -> &'static Seed42Fixture {
        static FIXTURE: OnceLock<Seed42Fixture> = OnceLock::new();
        FIXTURE.get_or_init(|| {
            let geometry =
                generate_plate_geometry(TEST_PLATES, 42, TEST_WARP_DEG, TEST_WIDTH, TEST_HEIGHT);
            let dynamics = compute_plate_dynamics(&geometry, 0.5, 42);
            let convergent_mask = dynamics
                .is_boundary
                .iter()
                .enumerate()
                .map(|(idx, &is_boundary)| {
                    is_boundary && dynamics.boundary_field[idx].convergent_rate > 1.0
                })
                .collect::<Vec<_>>();
            let convergent_components_raw = extract_components_from_mask(
                &dynamics.boundary_field,
                &convergent_mask,
                &geometry.plate_ids,
                TEST_WIDTH,
                TEST_HEIGHT,
            );
            let convergent_components_filtered = convergent_components_raw
                .iter()
                .filter(|component| component.pixels.len() >= MIN_COMPONENT_PIXELS)
                .cloned()
                .collect::<Vec<_>>();
            let largest_convergent_raw = convergent_components_raw
                .iter()
                .filter(|component| ordered_pair(component.plate_a, component.plate_b) == (3, 10))
                .max_by_key(|component| component.pixels.len())
                .expect("pair 3-10 convergent component")
                .clone();
            let largest_convergent_thinned =
                thin_component(&largest_convergent_raw, TEST_WIDTH, TEST_HEIGHT);
            let (largest_convergent_chain, is_closed) =
                order_skeleton_pixels(&largest_convergent_thinned, TEST_WIDTH, TEST_HEIGHT);
            let control_vertices = ordered_control_vertices(
                &largest_convergent_chain,
                &dynamics.boundary_field,
                &geometry.plate_ids,
                TEST_WIDTH,
                TEST_HEIGHT,
            );
            let smoothed = chaikin_smooth(&control_vertices, is_closed, CHAIKIN_ITERS);
            let (vertices, arc_lengths) =
                build_polyline_vertices(&smoothed, is_closed, TEST_WIDTH, TEST_HEIGHT);
            let largest_polyline = BoundaryPolyline {
                plate_a: 3,
                plate_b: 10,
                dominant_character: dominant_character(
                    &dynamics.boundary_field,
                    &largest_convergent_raw.pixels,
                ),
                vertices,
                is_closed,
                arc_lengths,
            };

            Seed42Fixture {
                boundary_field: dynamics.boundary_field,
                is_boundary: dynamics.is_boundary,
                plate_ids: geometry.plate_ids,
                convergent_components_raw,
                convergent_components_filtered,
                largest_convergent_raw,
                largest_convergent_thinned,
                largest_convergent_chain,
                largest_polyline,
            }
        })
    }

    #[test]
    fn component_extraction_matches_seed42_probe() {
        let fixture = seed42_fixture();
        assert_eq!(fixture.convergent_components_raw.len(), 17);
        assert_eq!(fixture.convergent_components_filtered.len(), 15);
    }

    #[test]
    fn thinning_produces_clean_skeleton_for_largest_convergent_component() {
        let fixture = seed42_fixture();
        assert_eq!(
            ordered_pair(
                fixture.largest_convergent_raw.plate_a,
                fixture.largest_convergent_raw.plate_b
            ),
            (3, 10)
        );
        assert_eq!(fixture.largest_convergent_raw.pixels.len(), 476);
        for &idx in &fixture.largest_convergent_thinned {
            let degree = neighbors8(idx, TEST_WIDTH, TEST_HEIGHT)
                .into_iter()
                .filter(|neighbor| fixture.largest_convergent_thinned.contains(neighbor))
                .count();
            assert!(degree <= 2, "skeleton degree {degree} at idx {idx}");
        }
        assert!(
            (110..=160).contains(&fixture.largest_convergent_thinned.len()),
            "expected thinned length in 110..=160, got {}",
            fixture.largest_convergent_thinned.len()
        );
    }

    #[test]
    fn chain_ordering_preserves_skeleton_connectivity() {
        let fixture = seed42_fixture();
        assert_eq!(
            fixture.largest_convergent_chain.len(),
            fixture.largest_convergent_thinned.len()
        );
        for window in fixture.largest_convergent_chain.windows(2) {
            assert!(
                neighbors8(window[0], TEST_WIDTH, TEST_HEIGHT).contains(&window[1]),
                "chain step {:?} is not 8-connected",
                window
            );
        }
    }

    #[test]
    fn chaikin_smoothing_increases_vertices_and_keeps_tangents_continuous() {
        let fixture = seed42_fixture();
        assert!(
            fixture.largest_polyline.vertices.len() > fixture.largest_convergent_chain.len(),
            "Chaikin smoothing should add vertices"
        );
        for window in fixture.largest_polyline.vertices.windows(2) {
            let a = window[0].tangent;
            let b = window[1].tangent;
            let dot = dot_2d(a, b).clamp(-1.0, 1.0);
            let angle_deg = dot.acos().to_degrees();
            assert!(
                angle_deg <= 45.0,
                "adjacent tangents deviate by {angle_deg:.2} degrees"
            );
        }
    }

    #[test]
    fn convergent_polyline_metadata_tracks_original_component_mean() {
        let fixture = seed42_fixture();
        assert_eq!(
            fixture.largest_polyline.dominant_character,
            BoundaryType::Convergent
        );
        assert!(fixture
            .largest_polyline
            .vertices
            .iter()
            .all(|vertex| vertex.convergent_rate > 0.0));

        let original_mean = fixture
            .largest_convergent_raw
            .pixels
            .iter()
            .map(|&idx| fixture.boundary_field[idx].convergent_rate)
            .sum::<f32>()
            / fixture.largest_convergent_raw.pixels.len() as f32;
        let polyline_mean = fixture
            .largest_polyline
            .vertices
            .iter()
            .map(|vertex| vertex.convergent_rate)
            .sum::<f32>()
            / fixture.largest_polyline.vertices.len() as f32;
        let relative_error = (polyline_mean - original_mean).abs() / original_mean.abs().max(1e-6);
        assert!(
            relative_error <= 0.10,
            "mean convergent rate drifted too far: original={original_mean:.3} polyline={polyline_mean:.3}"
        );
    }

    #[test]
    fn full_extraction_returns_polylines() {
        let fixture = seed42_fixture();
        let polylines = extract_boundary_polylines(
            &fixture.boundary_field,
            &fixture.is_boundary,
            &fixture.plate_ids,
            TEST_WIDTH,
            TEST_HEIGHT,
        );
        assert!(!polylines.is_empty());
    }
}
