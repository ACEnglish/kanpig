use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use crate::kplib::{metrics::GTstate, GenotypeAnno};

use noodles_vcf::{
    self as vcf,
    header::record::value::map::format,
    header::record::value::Map,
    record::genotypes::{keys::Keys, sample::Value, Genotypes},
};

pub struct VcfWriter {
    writer: vcf::Writer<BufWriter<File>>,
    header: vcf::Header,
    keys: Keys,
    pub gtcounts: HashMap<GTstate, usize>,
}

impl VcfWriter {
    /// Given a path and a header, setup a new output vcf
    pub fn new(out_path: &PathBuf, mut header: vcf::Header, sample: &Option<String>) -> Self {
        // Ensure sample is correctly setup
        let sample_name = match sample {
            Some(name) => name.clone(),
            None => {
                if header.sample_names().is_empty() {
                    error!("--input contains no samples. --sample name must be provided");
                    std::process::exit(1);
                }
                let samp_name = header.sample_names()[0].clone();
                info!("setting sample to {}", samp_name);
                samp_name
            }
        };

        if header.sample_names().len() > 1 {
            warn!(
                "clearing {} sample columns in output",
                header.sample_names().len()
            );
            header.sample_names_mut().clear()
        }
        header.sample_names_mut().insert(sample_name);

        // Setup FORMAT header definitions
        // Overwrites existing definitions
        let all_formats = header.formats_mut();

        let keys: Keys = "GT:FT:SQ:GQ:PG:DP:AD:ZS:SS".parse().unwrap();

        // GT
        let gtid = keys[0].clone();
        let mut gtfmt = Map::<format::Format>::from(&gtid);
        *gtfmt.number_mut() = vcf::header::Number::Count(1);
        *gtfmt.type_mut() = format::Type::String;
        *gtfmt.description_mut() = "Kanplug genotype".to_string();
        all_formats.insert(gtid, gtfmt);

        // FT
        let ftid = keys[1].clone();
        let mut ftfmt = Map::<format::Format>::from(&ftid);
        *ftfmt.number_mut() = vcf::header::Number::Count(1);
        *ftfmt.type_mut() = format::Type::Integer;
        *ftfmt.description_mut() = "Kanpig filter".to_string();
        all_formats.insert(ftid, ftfmt);

        // SQ
        let sqid = keys[2].clone();
        let mut sqfmt = Map::<format::Format>::from(&sqid);
        *sqfmt.number_mut() = vcf::header::Number::Count(1);
        *sqfmt.type_mut() = format::Type::Integer;
        *sqfmt.description_mut() =
            "Phred scaled quality of sample being non-ref at this variant".to_string();
        all_formats.insert(sqid, sqfmt);

        // GQ
        let gqid = keys[3].clone();
        let mut gqfmt = Map::<format::Format>::from(&gqid);
        *gqfmt.number_mut() = vcf::header::Number::Count(1);
        *gqfmt.type_mut() = format::Type::Integer;
        *gqfmt.description_mut() = "Phred scaled quality of genotype".to_string();
        all_formats.insert(gqid, gqfmt);

        // PG
        let pgid = keys[4].clone();
        let mut pgfmt = Map::<format::Format>::from(&pgid);
        *pgfmt.number_mut() = vcf::header::Number::Count(1);
        *pgfmt.type_mut() = format::Type::Integer;
        *pgfmt.description_mut() = "Local phase group of entries".to_string();
        all_formats.insert(pgid, pgfmt);

        // DP
        let dpid = keys[5].clone();
        let mut dpfmt = Map::<format::Format>::from(&dpid);
        *dpfmt.number_mut() = vcf::header::Number::Count(1);
        *dpfmt.type_mut() = format::Type::Integer;
        *dpfmt.description_mut() = "Coverage over region".to_string();
        all_formats.insert(dpid, dpfmt);

        // AD
        let adid = keys[6].clone();
        let mut adfmt = Map::<format::Format>::from(&adid);
        *adfmt.number_mut() = vcf::header::Number::R;
        *adfmt.type_mut() = format::Type::Integer;
        *adfmt.description_mut() = "Coverage for reference and alternate alleles".to_string();
        all_formats.insert(adid, adfmt);

        // ZS
        let zsid = keys[7].clone();
        let mut zsfmt = Map::<format::Format>::from(&zsid);
        *zsfmt.number_mut() = vcf::header::Number::R;
        *zsfmt.type_mut() = format::Type::Integer;
        *zsfmt.description_mut() = "Size similarity of path to entry".to_string();
        all_formats.insert(zsid, zsfmt);

        // SS
        let ssid = keys[8].clone();
        let mut ssfmt = Map::<format::Format>::from(&ssid);
        *ssfmt.number_mut() = vcf::header::Number::R;
        *ssfmt.type_mut() = format::Type::Integer;
        *ssfmt.description_mut() = "Sequence similarity of path to entry".to_string();
        all_formats.insert(ssid, ssfmt);

        // Ready to make files
        let out_buf = BufWriter::with_capacity(
            page_size::get() * 1000,
            File::create(out_path).expect("Error Creating Output File"),
        );
        let mut writer = vcf::Writer::new(out_buf);
        let _ = writer.write_header(&header);

        Self {
            writer,
            header,
            keys,
            gtcounts: HashMap::new(),
        }
    }

    pub fn anno_write(&mut self, mut annot: GenotypeAnno, phase_group: i32) {
        *self.gtcounts.entry(annot.gt_state).or_insert(0) += 1;
        *annot.entry.genotypes_mut() =
            Genotypes::new(self.keys.clone(), vec![annot.make_fields(phase_group)]);

        let _result = self.writer.write_record(&self.header, &annot.entry);
    }

    pub fn write_entry(&mut self, mut entry: vcf::Record) {
        *entry.genotypes_mut() = Genotypes::new(
            "GT".parse().unwrap(),
            vec![vec![
                Some(Value::from("./.")), // GT
            ]],
        );
        *self.gtcounts.entry(GTstate::Non).or_insert(0) += 1;
        let _result = self.writer.write_record(&self.header, &entry);
    }
}
