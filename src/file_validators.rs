use crate::{commands::plup::PlupCommand, kplib::KDParams};
use rust_htslib::tbx::{self, Read as TbxRead};
use std::path::Path;

/// Helper function to validate a file's existence and type
pub fn validate_file(path: &Path, label: &str) -> bool {
    if !path.exists() {
        error!("{} does not exist", label);
        return false;
    }
    if !path.is_file() {
        error!("{} is not a file", label);
        return false;
    }
    true
}

pub fn validate_bam(file_path: &str) -> bool {
    let mut is_ok = true;
    if file_path.ends_with(".bam") || file_path.ends_with(".cram") {
        let index_extensions = [".bai", ".crai", ".csi"];
        let index_exists = index_extensions.iter().any(|ext| {
            let index_path = format!("{}{}", file_path, ext);
            let p = Path::new(&index_path);
            p.exists() & p.is_file()
        });

        if !index_exists {
            error!(
                "bam/cram index ({}) does not exist",
                index_extensions.join(", ")
            );
            is_ok = false;
        }
    } else {
        is_ok = false;
    }
    is_ok
}

pub fn validate_plup(file_path: &str, params: &KDParams) -> bool {
    let mut is_ok = true;
    if !file_path.ends_with(".plup.gz") {
        is_ok = false;
    } else {
        let tbi_path = format!("{}.tbi", file_path);
        if !validate_file(Path::new(&tbi_path), "plup index (.tbi)") {
            is_ok = false;
        } else {
            let tbx = tbx::Reader::from_path(file_path).expect("Failed to open TBX file");
            let header = tbx.header();
            if header.len() != 1 {
                error!("Malformed plup.gz header. Unable to validate parameters");
            } else {
                match serde_json::from_str::<PlupCommand>(&header[0][2..]) {
                    Ok(plup_args) => {
                        if plup_args.sizemin != params.sizemin {
                            warn!(
                                "plup created with --sizemin {} != gt --sizemin {}",
                                plup_args.sizemin, params.sizemin
                            );
                        }

                        if plup_args.sizemax != params.sizemax {
                            warn!(
                                "plup created with --sizemax {} != gt --sizemax {}",
                                plup_args.sizemax, params.sizemax
                            );
                        }

                        if plup_args.mapq != params.mapq {
                            warn!(
                                "plup created with --mapq {} != gt --mapq {}",
                                plup_args.mapq, params.mapq
                            );
                        }

                        if plup_args.mapflag != params.mapflag {
                            warn!(
                                "plup created with --mapflag {} != gt --mapflag {}",
                                plup_args.mapflag, params.mapflag
                            );
                        }
                    }
                    Err(e) => {
                        error!(
                            "Failed to parse plup.gz header for parameter validation: {}",
                            e
                        );
                    }
                }
            }
        }
    }
    is_ok
}
/// Helper function to validate reads (.bam, .cram, or .plup.gz)
pub fn validate_reads(reads: &Path, params: &KDParams) -> bool {
    let mut is_ok = validate_file(reads, "--reads");
    let file_path = reads.to_str().unwrap_or_default();
    let bam_ok = validate_bam(file_path);
    let plup_ok = validate_plup(file_path, params);
    if !(bam_ok || plup_ok) {
        error!("Unsupported file type: {}", file_path);
        is_ok = false;
    }
    is_ok
}

/// Checks reference and its .fai index
pub fn validate_reference(reference: &Path) -> bool {
    let mut is_ok = validate_file(reference, "--reference");

    let mut fai_path = reference.to_path_buf();
    fai_path.set_file_name(format!(
        "{}.fai",
        fai_path.file_name().unwrap().to_string_lossy()
    ));
    is_ok &= validate_file(&fai_path, "--reference index (.fai)");

    is_ok
}
