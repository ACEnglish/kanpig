/*
* This needs the output pathbuf
* and a vcf header hold a copy since I think the vcfparser needs it, also
* Actually, def a copy because we're altering it
*
* Then we check if the sample name is in the header
* if yes, we want to preserve that format information
* if not, we just warn that they're getting a new column
* Lets be nice and warn that N samples will be dropped in the output vcf (for now because I don't
* want to do record keeping)
*
* Format fields we'll populate are
*  GT
*  PG : phase group of a chunk. this will just be an integer and for a single sample. The writer
*  will know what number its on
*  DP : total coverage depth over the analyzed region
*  AD : allele 1 and allele 2 depth
*  --- Advanced ---
*  GQ : (probably) Essentially copying what I did for
*  PL : (I think -- what was this other field?) geno quality a single number and then
*          per-combination, or something
*  FT : filter for the different cases I might want to look out for e.g. lowcoverage, low
*      similarity, low gq (which is just a non thresholded low coverage) I might make this a
*      bit flag just so they take up less space.
*
* has a method for write() which will do the annotation work
* before calling `writer.write_record`
*
*  "GT", genotype
   "PG", local phase group
   "DP", observed depth over region
   "AD", allele supporting depth (reference is assumed)
* Once this is setup, we'll work on the h1/h2 path finding where we stop trying to search for het
* paths (n=0) and just used the reported overall coverage to populate DP/AD
*/
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use petgraph::graph::NodeIndex;

use crate::kanpig::PathScore;

use noodles_vcf::{
    self as vcf,
    header::record::value::map::format,
    header::record::value::Map,
    record::genotypes::{keys::Keys, sample::Value, Genotypes},
};

pub struct VcfWriter {
    writer: vcf::Writer<BufWriter<File>>,
    header: vcf::Header,
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

        let keys: Keys = "GT:PG:DP:AD:ZS:SS".parse().unwrap(); // I shouldn't be making this each time?

        // GT
        let gtid = keys[0].clone();
        let mut gtfmt = Map::<format::Format>::from(&gtid);
        *gtfmt.number_mut() = vcf::header::Number::Count(1);
        *gtfmt.type_mut() = format::Type::String;
        *gtfmt.description_mut() = "Kanplug genotype".to_string();
        all_formats.insert(gtid, gtfmt);

        // PG
        let pgid = keys[1].clone();
        let mut pgfmt = Map::<format::Format>::from(&pgid);
        *pgfmt.number_mut() = vcf::header::Number::Count(1);
        *pgfmt.type_mut() = format::Type::Integer;
        *pgfmt.description_mut() = "Local phase group of entries".to_string();
        all_formats.insert(pgid, pgfmt);

        // DP
        let dpid = keys[2].clone();
        let mut dpfmt = Map::<format::Format>::from(&dpid);
        *dpfmt.number_mut() = vcf::header::Number::Count(1);
        *dpfmt.type_mut() = format::Type::Integer;
        *dpfmt.description_mut() = "Coverage over region".to_string();
        all_formats.insert(dpid, dpfmt);

        // AD
        let adid = keys[3].clone();
        let mut adfmt = Map::<format::Format>::from(&adid);
        *adfmt.number_mut() = vcf::header::Number::R;
        *adfmt.type_mut() = format::Type::Integer;
        *adfmt.description_mut() = "Coverage per allele".to_string();
        all_formats.insert(adid, adfmt);

        // ZS
        let zsid = keys[4].clone();
        let mut zsfmt = Map::<format::Format>::from(&zsid);
        *zsfmt.number_mut() = vcf::header::Number::R;
        *zsfmt.type_mut() = format::Type::Integer;
        *zsfmt.description_mut() = "Size similarity of path to entry".to_string();
        all_formats.insert(zsid, zsfmt);

        // SS
        let ssid = keys[5].clone();
        let mut ssfmt = Map::<format::Format>::from(&ssid);
        *ssfmt.number_mut() = vcf::header::Number::R;
        *ssfmt.type_mut() = format::Type::Integer;
        *ssfmt.description_mut() = "Sequence similarity of path to entry".to_string();
        all_formats.insert(ssid, ssfmt);

        // Ready to make files
        let out_buf = BufWriter::new(File::create(out_path).expect("Error Creating Output File"));
        let mut writer = vcf::Writer::new(out_buf);
        let _ = writer.write_header(&header);

        Self { writer, header }
    }

    pub fn anno_write(
        &mut self,
        mut entry: vcf::Record,
        var_idx: &NodeIndex,
        path1: &PathScore,
        path2: &PathScore,
        coverage: u64,
        phase_group: i32,
    ) {
        //let keys = "GT:PG:DP:AD".parse().unwrap(); // I shouldn't be making this each time?
        let keys = "GT:PG:DP:AD:ZS:SS".parse().unwrap(); // I shouldn't be making this each time?

        //https://docs.rs/noodles-vcf/0.49.0/noodles_vcf/record/genotypes/struct.Genotypes.html
        let gt = match (path1.path.contains(var_idx), path2.path.contains(var_idx)) {
            (true, true) => "1|1",
            (true, false) => "1|0",
            (false, true) => "0|1",
            (false, false) => "0|0",
        };
        // https://docs.rs/noodles-vcf/0.49.0/noodles_vcf/header/record/value/map/struct.Map.html
        let ad = Value::from(vec![
            Some(path1.coverage.unwrap() as i32),
            Some(path2.coverage.unwrap() as i32),
        ]);

        let zs = Value::from(vec![
            Some((path1.sizesim * 100.0) as i32),
            Some((path2.sizesim * 100.0) as i32),
        ]);

        let ss = Value::from(vec![
            Some((path1.seqsim * 100.0) as i32),
            Some((path2.seqsim * 100.0) as i32),
        ]);

        let genotypes = Genotypes::new(
            keys,
            vec![vec![
                Some(Value::from(gt)),              // GT
                Some(Value::from(phase_group)),     // PG
                Some(Value::from(coverage as i32)), // DP
                Some(ad),
                Some(zs),
                Some(ss),
            ]],
        );

        *entry.genotypes_mut() = genotypes.clone();
        let _result = self.writer.write_record(&self.header, &entry);
    }
}
