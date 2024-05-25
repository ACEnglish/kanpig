use crate::kplib::{metrics::GTstate, GenotypeAnno};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use noodles_vcf::{
    self as vcf,
    header::record::value::map::format,
    header::record::value::Map,
    variant::io::Write,
    variant::record_buf::samples::{keys::Keys, Samples},
};

pub struct VcfWriter {
    writer: vcf::io::Writer<BufWriter<File>>,
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

        if !header.sample_names().is_empty() {
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
        let new_fmts: Vec<String> = "GT:FT:SQ:GQ:PS:DP:AD:ZS:SS"
            .split(':')
            .map(String::from)
            .collect();
        let keys: Keys = Keys::from_iter(new_fmts);

        // GT
        let gtid = "GT";
        let mut gtfmt = Map::<format::Format>::from(gtid);
        *gtfmt.number_mut() = format::Number::Count(1);
        *gtfmt.type_mut() = format::Type::String;
        *gtfmt.description_mut() = "Kanplug genotype".to_string();
        all_formats.insert(gtid.to_string(), gtfmt);

        // FT
        let ftid = "FT";
        let mut ftfmt = Map::<format::Format>::from(ftid);
        *ftfmt.number_mut() = format::Number::Count(1);
        *ftfmt.type_mut() = format::Type::Integer;
        *ftfmt.description_mut() = "Kanpig filter".to_string();
        all_formats.insert(ftid.to_string(), ftfmt);

        // SQ
        let sqid = "SQ";
        let mut sqfmt = Map::<format::Format>::from(sqid);
        *sqfmt.number_mut() = format::Number::Count(1);
        *sqfmt.type_mut() = format::Type::Integer;
        *sqfmt.description_mut() =
            "Phred scaled quality of sample being non-ref at this variant".to_string();
        all_formats.insert(sqid.to_string(), sqfmt);

        // GQ
        let gqid = "GQ";
        let mut gqfmt = Map::<format::Format>::from(gqid);
        *gqfmt.number_mut() = format::Number::Count(1);
        *gqfmt.type_mut() = format::Type::Integer;
        *gqfmt.description_mut() = "Phred scaled quality of genotype".to_string();
        all_formats.insert(gqid.to_string(), gqfmt);

        // PS
        let psid = "PS";
        let mut psfmt = Map::<format::Format>::from(psid);
        *psfmt.number_mut() = format::Number::Count(1);
        *psfmt.type_mut() = format::Type::Integer;
        *psfmt.description_mut() = "Local phase group of entries".to_string();
        all_formats.insert(psid.to_string(), psfmt);

        // DP
        let dpid = "DP";
        let mut dpfmt = Map::<format::Format>::from(dpid);
        *dpfmt.number_mut() = format::Number::Count(1);
        *dpfmt.type_mut() = format::Type::Integer;
        *dpfmt.description_mut() = "Coverage over region".to_string();
        all_formats.insert(dpid.to_string(), dpfmt);

        // AD
        let adid = "AD";
        let mut adfmt = Map::<format::Format>::from(adid);
        *adfmt.number_mut() = format::Number::ReferenceAlternateBases;
        *adfmt.type_mut() = format::Type::Integer;
        *adfmt.description_mut() = "Coverage for reference and alternate alleles".to_string();
        all_formats.insert(adid.to_string(), adfmt);

        // ZS
        let zsid = "ZS";
        let mut zsfmt = Map::<format::Format>::from(zsid);
        *zsfmt.number_mut() = format::Number::ReferenceAlternateBases;
        *zsfmt.type_mut() = format::Type::Integer;
        *zsfmt.description_mut() = "Size similarity of path to entry".to_string();
        all_formats.insert(zsid.to_string(), zsfmt);

        // SS
        let ssid = "SS";
        let mut ssfmt = Map::<format::Format>::from(ssid);
        *ssfmt.number_mut() = format::Number::ReferenceAlternateBases;
        *ssfmt.type_mut() = format::Type::Integer;
        *ssfmt.description_mut() = "Sequence similarity of path to entry".to_string();
        all_formats.insert(ssid.to_string(), ssfmt);

        // Ready to make files
        let out_buf = BufWriter::with_capacity(
            page_size::get() * 500,
            File::create(out_path).expect("Error Creating Output File"),
        );
        let mut writer = vcf::io::Writer::new(out_buf);
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
        *annot.entry.samples_mut() =
            Samples::new(self.keys.clone(), vec![annot.make_fields(phase_group)]);
        //*annot.entry.reference_bases_mut() = replace_iupac(annot.entry.reference_bases());
        replace_iupac_inplace(annot.entry.reference_bases_mut());
        match self.writer.write_variant_record(&self.header, &annot.entry) {
            Ok(_) => {}
            Err(error) => panic!("Couldn't write record {:?}", error),
        }
    }
}

lazy_static::lazy_static! {
    static ref IUPAC: [u8; 128] = {
        let mut arr = [0u8; 128];
        arr[b'R' as usize] = b'A'; // A or G
        arr[b'Y' as usize] = b'C'; // C or T
        arr[b'S' as usize] = b'C'; // C or G
        arr[b'W' as usize] = b'A'; // A or T
        arr[b'K' as usize] = b'G'; // G or T
        arr[b'M' as usize] = b'A'; // A or C
        arr[b'B' as usize] = b'C'; // C, G, or T
        arr[b'D' as usize] = b'A'; // A, G, or T
        arr[b'H' as usize] = b'A'; // A, C, or T
        arr[b'V' as usize] = b'A'; // A, C, or G
        arr[b'N' as usize] = b'A'; // Any base (A, C, G, or T)
        arr[b'r' as usize] = b'a'; // a or g
        arr[b'y' as usize] = b'c'; // c or t
        arr[b's' as usize] = b'c'; // c or g
        arr[b'w' as usize] = b'a'; // a or t
        arr[b'k' as usize] = b'g'; // g or t
        arr[b'm' as usize] = b'a'; // a or c
        arr[b'b' as usize] = b'c'; // c, g, or t
        arr[b'd' as usize] = b'a'; // a, g, or t
        arr[b'h' as usize] = b'a'; // a, c, or t
        arr[b'v' as usize] = b'a'; // a, c, or g
        arr[b'n' as usize] = b'a'; // any base (a, c, g, or t)
        arr
    };
}

fn replace_iupac_inplace(sequence: &mut String) {
    unsafe {
        let bytes = sequence.as_bytes_mut();
        bytes.iter_mut().for_each(|b| {
            let t = IUPAC[*b as usize]; 
            if t != 0u8 {
                *b = t;
            }
        });
    }
}

/*fn replace_iupac(sequence: &str) -> String {
    let mut result = String::with_capacity(sequence.len());

    for c in sequence.chars() {
        result.push(REPLACEMENTS[c as usize]);
    }

    result
}*/
