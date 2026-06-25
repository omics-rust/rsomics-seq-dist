use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};

use crate::{Metric, run};

#[derive(Parser)]
#[command(
    name = "rsomics-seq-dist",
    version,
    about = "Pairwise nucleotide distances (Hamming, JC69, K2P) from an aligned FASTA"
)]
pub struct Cli {
    /// Aligned FASTA, or - for stdin.
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// Distance metric.
    #[arg(long, value_parser = ["hamming", "jc69", "k2p"], default_value = "hamming")]
    pub metric: String,

    /// Output distance-matrix TSV, or - for stdout.
    #[arg(short = 'o', long, default_value = "-")]
    pub output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        ToolMeta {
            name: "rsomics-seq-dist",
            version: env!("CARGO_PKG_VERSION"),
        }
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let metric: Metric = self.metric.parse()?;

        let mut reader: Box<dyn io::BufRead> = if self.input.as_os_str() == "-" {
            Box::new(BufReader::new(io::stdin().lock()))
        } else {
            Box::new(BufReader::new(File::open(&self.input).map_err(|e| {
                RsomicsError::InvalidInput(format!("{}: {e}", self.input.display()))
            })?))
        };

        let mut writer: Box<dyn io::Write> = if self.output == "-" {
            Box::new(BufWriter::new(io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };

        run(metric, &mut reader, &mut writer)
    }
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
