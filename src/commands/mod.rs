use clap::Subcommand;

/// Trait required for each of the kanpig commands
pub trait KanpigCommand: std::fmt::Debug {
    fn validate(&self) -> bool;
    fn debug(&self) -> bool;
    fn run(&mut self);
}

pub mod germ;
use crate::commands::germ::GermCommand;

pub mod plup;
use crate::commands::plup::PlupCommand;

// pub mod trio;
// use crate::commands::trio::TrioCommand;

// pub mod mosaic;
// use crate::commands::mosaic::MosaicCommand;

/// Set of commands
#[derive(Subcommand, Debug, Clone)]
pub enum Commands {
    #[command(about = "BAM/CRAM to Pileup Index")]
    Plup(PlupCommand),

    #[command(about = "Germline SV Genotyping")]
    GT(GermCommand),
    // #[command(about = "Trio SV Genotyping")]
    // Trio(TrioCommand),

    // #[command(about = "Mosaic SV Genotyping")]
    // Mosaic(MosaicCommand),
}
