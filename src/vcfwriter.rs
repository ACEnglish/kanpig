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

use crate::PathScore;

use noodles_vcf::{
    self as vcf,
    record::genotypes::{sample::Value, Genotypes},
};

pub struct VcfWriter {
    writer: vcf::Writer<BufWriter<File>>,
    header: vcf::Header,
}

impl VcfWriter {
    pub fn new(out_path: &PathBuf, header: vcf::Header) -> Self {
        let out_buf = BufWriter::new(File::create(out_path).expect("Error Creating Output File"));
        let mut writer = vcf::Writer::new(out_buf);

        // Edit header here

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
    ) {
        //let keys = "GT:PG:DP:AD".parse().unwrap(); // I shouldn't be making this each time?
        let keys = "GT:DP:AD".parse().unwrap(); // I shouldn't be making this each time?

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
        let genotypes = Genotypes::new(
            keys,
            vec![vec![
                Some(Value::from(gt)), // GT
                //Some(Value::from(0)),    // PG
                Some(Value::from(coverage as i32)), // DP
                Some(ad),
            ]],
        );

        *entry.genotypes_mut() = genotypes.clone();
        let _result = self.writer.write_record(&self.header, &entry);
    }
}
