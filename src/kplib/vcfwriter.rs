use crate::kplib::{metrics::GTstate, GenotypeAnno};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};

use noodles_vcf::{
    self as vcf,
    header::record::value::map::format,
    header::record::value::Map,
    variant::io::Write as vcfWrite,
    variant::record_buf::samples::{keys::Keys, Samples},
};

pub struct VcfWriter {
    writer: vcf::io::Writer<Box<dyn Write>>,
    header: vcf::Header,
    keys: Keys,
    pub gtcounts: HashMap<GTstate, usize>,
    pub iupac_fixed: bool,
    buf: Vec<u8>,
}

#[rustfmt::skip]
impl VcfWriter {
    /// Given a path and a header, setup a new output VCF
    pub fn new(
        out_path: &Option<PathBuf>,
        mut header: vcf::Header,
        sample: &Option<String>,
    ) -> Self {
        // Ensure sample is correctly set up
        let sample_name = match sample {
            Some(name) => name.clone(),
            None => {
                if header.sample_names().is_empty() {
                    error!("--input contains no samples. --sample name must be provided");
                    std::process::exit(1);
                }
                let samp_name = header.sample_names()[0].clone();
                info!("Setting sample to {}", samp_name);
                samp_name
            }
        };

        if !header.sample_names().is_empty() {
            warn!(
                "Clearing {} sample columns in output",
                header.sample_names().len()
            );
            header.sample_names_mut().clear();
        }
        header.sample_names_mut().insert(sample_name);

        // Setup FORMAT header definitions
        let all_formats = header.formats_mut();
        let num1 = format::Number::Count(1);
        // Edits to these must be sync'd with GenotypeAnno::make_fields
        let format_definitions = vec![
            ("GT", num1, format::Type::String, "Kanplug genotype"),
            ("FT", num1, format::Type::Integer, "Kanpig filter"),
            ("SQ", num1, format::Type::Integer, "Phred quality of being non-ref"),
            ("GQ", num1, format::Type::Integer, "Phred quality of genotype"),
            ("PS", num1, format::Type::Integer, "PhaseSet tag from reads"),
            ("NE", num1, format::Type::Integer, "Neighborhood phase set of entries"),
            ("DP", num1, format::Type::Integer, "Coverage over region"),
            ("AD", format::Number::ReferenceAlternateBases, format::Type::Integer, "Ref/Alt coverage"),
            ("KS", format::Number::Unknown, format::Type::Integer, "Kanpig score"),
        ];
        let new_fmts: Vec<String> = format_definitions.iter().map(|x| String::from(x.0)).collect();

        for (id, number, ty, desc) in format_definitions {
            all_formats.insert(id.to_string(), create_format(id, number, ty, desc));
        }

        // Prepare output
        let out_buf: Box<dyn Write> = match out_path {
            Some(ref path) => {
                let m_page = page_size::get() * 1000;
                let file = File::create(path).expect("Error creating output file");
                Box::new(BufWriter::with_capacity(m_page, file))
            }
            None => Box::new(BufWriter::new(std::io::stdout())),
        };
        let mut writer = vcf::io::Writer::new(out_buf);
        let _ = writer.write_header(&header);
        
        Self {
            writer,
            header,
            keys: Keys::from_iter(new_fmts),
            gtcounts: HashMap::new(),
            iupac_fixed: false,
            buf: vec![],
        }
    }

    pub fn anno_write(&mut self, mut annot: GenotypeAnno, neigh_id: i32) {
        *self.gtcounts.entry(annot.gt_state).or_insert(0) += 1;
        *annot.entry.samples_mut() =
            Samples::new(self.keys.clone(), vec![annot.make_fields(neigh_id)]);

        self.buf.clear();
        let mut tmp = vcf::io::Writer::new(&mut self.buf);
        if tmp.write_variant_record(&self.header, &annot.entry).is_err() {
            let changed = replace_iupac_inplace(annot.entry.reference_bases_mut());
            self.iupac_fixed |= changed;
            if let Err(error) = self.writer.write_variant_record(&self.header, &annot.entry) {
                panic!("Couldn't write record {:?}", error);
            }
        } else if let Err(error) = self.writer.get_mut().write_all(&self.buf) {
            panic!("Couldn't write record {:?}", error);
        }
    }
}

fn create_format(
    id: &str,
    number: format::Number,
    ty: format::Type,
    desc: &str,
) -> Map<format::Format> {
    let mut fmt = Map::<format::Format>::from(id);
    *fmt.number_mut() = number;
    *fmt.type_mut() = ty;
    *fmt.description_mut() = desc.to_string();
    fmt
}

lazy_static::lazy_static! {
    static ref IUPAC: [u8; 128] = {
        let mut arr = [0u8; 128];
        for &(iupac, replacement) in &[
            (b'R', b'A'), (b'Y', b'C'), (b'S', b'C'), (b'W', b'A'),
            (b'K', b'G'), (b'M', b'A'), (b'B', b'C'), (b'D', b'A'),
            (b'H', b'A'), (b'V', b'A'), (b'r', b'a'), (b'y', b'c'),
            (b's', b'c'), (b'w', b'a'), (b'k', b'g'), (b'm', b'a'),
            (b'b', b'c'), (b'd', b'a'), (b'h', b'a'), (b'v', b'a'),
        ] {
            arr[iupac as usize] = replacement;
        }
        arr
    };
}

fn replace_iupac_inplace(sequence: &mut str) -> bool {
    let mut any_change = false;
    unsafe {
        let bytes = sequence.as_bytes_mut();
        bytes.iter_mut().for_each(|b| {
            let t = IUPAC[*b as usize];
            if t != 0u8 {
                any_change = true;
                *b = t;
            }
        });
    }
    any_change
}
