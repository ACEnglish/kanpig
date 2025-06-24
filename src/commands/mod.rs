use clap::Subcommand;

/// Trait required for each of the kanpig commands
pub trait KanpigCommand: std::fmt::Debug {
    fn validate(&self) -> bool;
    fn debug(&self) -> bool;
    fn run(&mut self);
}

pub mod genotype;
use crate::commands::genotype::GTCommand;

pub mod plup;
use crate::commands::plup::PlupCommand;

pub mod trio;
use crate::commands::trio::TrioCommand;

/// Set of commands
#[derive(Subcommand, Debug, Clone)]
pub enum Commands {
    #[command(about = "Genotype SVs")]
    Gt(GTCommand),

    #[command(about = "BAM/CRAM to Pileup Index")]
    Plup(PlupCommand),

    #[command(about = "Trio Genotyping")]
    Trio(TrioCommand),
}
