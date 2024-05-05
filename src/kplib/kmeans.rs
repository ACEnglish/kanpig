// Written by chatgpt because the kmeans crate is broken because simd is??
use std::cmp::Ordering;

pub type Point = Vec<f32>;
pub type Centroid = Point;

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

fn choose_centroids(sorted_vec: &[Point], k: usize) -> Vec<Point> {
    let segment_size = sorted_vec.len() / k;
    let mut centroids = Vec::with_capacity(k);

    for i in 0..k {
        let centroid_index = i * segment_size + segment_size / 2;
        centroids.push(sorted_vec[centroid_index].clone());
    }

    centroids
}

pub fn kmeans(data: &[Point], k: usize) -> Vec<Cluster> {
    let mut centroids: Vec<Point> = choose_centroids(data, k); //data.iter().take(k).cloned().collect();
    let mut clusters: Vec<Cluster> = centroids
        .iter()
        .map(|centroid| Cluster::new(centroid.clone()))
        .collect();

    let mut ntries = 10;
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
        if new_centroids == old_centroids || ntries == 0 {
            break;
        } else {
            centroids = new_centroids;
            clusters = centroids
                .iter()
                .map(|centroid| Cluster::new(centroid.clone()))
                .collect();
        }
        ntries -= 1;
    }

    clusters
}
