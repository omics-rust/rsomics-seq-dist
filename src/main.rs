use clap::Parser;
use rsomics_common::Tool;
use rsomics_seq_dist::cli::Cli;

fn main() -> std::process::ExitCode {
    Cli::parse().run()
}
