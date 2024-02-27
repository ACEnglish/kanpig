// Written by chatgpt because the kmeans crate is broken because simd is??
use rand::seq::SliceRandom;
use std::cmp::Ordering;

// Define a point in multidimensional space
pub type Point = Vec<f32>;

// Define a cluster centroid
pub type Centroid = Point;

// Define a cluster containing points
#[derive(Debug)]
pub struct Cluster {
    pub centroid: Centroid,
    pub points: Vec<Point>,
    pub points_idx: Vec<usize>,
}

impl Cluster {
    fn new(centroid: Centroid) -> Self {
        Self {
            centroid,
            points: Vec::new(),
            points_idx: Vec::new(),
        }
    }

    fn update_centroid(&mut self) {
        if !self.points.is_empty() {
            let dim = self.points[0].len();
            let mut new_centroid = vec![0.0; dim];

            for point in &self.points {
                for (i, &coord) in point.iter().enumerate() {
                    new_centroid[i] += coord;
                }
            }

            for coord in &mut new_centroid {
                *coord /= self.points.len() as f32;
            }

            self.centroid = new_centroid;
        }
    }
}

fn distance(p1: &Point, p2: &Point) -> f32 {
    p1.iter()
        .zip(p2)
        .map(|(&a, &b)| (a - b).powi(2))
        .sum::<f32>()
        .sqrt()
}

pub fn kmeans(data: &[Point], k: usize) -> Vec<Cluster> {
    // Initialize clusters with random centroids
    let mut rng = rand::thread_rng();
    // Should probably choose highest covered
    let mut centroids: Vec<Point> = data.choose_multiple(&mut rng, k).cloned().collect();
    let mut clusters: Vec<Cluster> = centroids
        .iter()
        .map(|centroid| Cluster::new(centroid.clone()))
        .collect();

    loop {
        // Assign each point to the nearest cluster
        for (idx, point) in data.iter().enumerate() {
            let nearest_cluster_idx = clusters
                .iter()
                .enumerate()
                .min_by(|(_, c1), (_, c2)| {
                    distance(point, &c1.centroid)
                        .partial_cmp(&distance(point, &c2.centroid))
                        .unwrap_or(Ordering::Equal)
                })
                .map(|(idx, _)| idx)
                .unwrap();

            clusters[nearest_cluster_idx].points.push(point.clone());
            clusters[nearest_cluster_idx].points_idx.push(idx);
        }

        // Update cluster centroids
        let old_centroids: Vec<Centroid> = clusters.iter().map(|c| c.centroid.clone()).collect();
        for (cluster, _centroid) in clusters.iter_mut().zip(&old_centroids) {
            cluster.update_centroid();
        }

        // Check for convergence
        let new_centroids: Vec<Point> = clusters.iter().map(|c| c.centroid.clone()).collect();
        if new_centroids == old_centroids {
            break;
        } else {
            centroids = new_centroids;
            clusters = centroids
                .iter()
                .map(|centroid| Cluster::new(centroid.clone()))
                .collect();
        }
    }

    clusters
}
