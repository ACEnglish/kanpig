use crate::kplib::{metrics::GTstate, ChannelOutput, GenotypeAnno, Ploidy};
use crossbeam_channel::Receiver;
use indicatif::{ProgressBar, ProgressStyle};
use petgraph::graph::NodeIndex;
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
    sample_count: usize,
    header: vcf::Header,
    keys: Keys,
    pub gtcounts: HashMap<GTstate, usize>,
    pub iupac_fixed: bool,
    buf: Vec<u8>,
}

impl VcfWriter {
    /// Given a path and a header, setup a new output VCF
    pub fn new(
        out_path: &Option<PathBuf>,
        mut header: vcf::Header,
        sample_names: &Vec<String>,
    ) -> Self {
        if !header.sample_names().is_empty() {
            warn!(
                "Clearing {} sample columns in output",
                header.sample_names().len()
            );
            header.sample_names_mut().clear();
        }
        for i in sample_names {
            header.sample_names_mut().insert(i.to_owned());
        }

        // Setup FORMAT header definitions
        let all_formats = header.formats_mut();
        all_formats.clear();
        let num1 = format::Number::Count(1);
        // Edits to these must be sync'd with GenotypeAnno::make_fields
        #[rustfmt::skip]
        let format_definitions = vec![
            ("GT", num1, format::Type::String, "Kanpig genotype"),
            ("FT", num1, format::Type::Integer, "Kanpig filter"),
            ("SQ", num1, format::Type::Integer, "Phred quality of being non-ref"),
            ("GQ", num1, format::Type::Integer, "Phred quality of genotype"),
            ("PS", num1, format::Type::Integer, "PhaseSet tag from reads"),
            ("NE", num1, format::Type::Integer, "Neighborhood phase set of entries"),
            ("DP", num1, format::Type::Integer, "Coverage over region"),
            ("AD", format::Number::ReferenceAlternateBases, format::Type::Integer, "Ref/Alt coverage"),
            ("KS", format::Number::Unknown, format::Type::Integer, "Kanpig score"),
        ];
        let new_fmts: Vec<String> = format_definitions
            .iter()
            .map(|x| String::from(x.0))
            .collect();

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
            sample_count: sample_names.len(),
            keys: Keys::from_iter(new_fmts),
            gtcounts: HashMap::new(),
            iupac_fixed: false,
            buf: vec![],
        }
    }

    pub fn anno_write(
        &mut self,
        mut entry: vcf::variant::RecordBuf,
        mut annots: Vec<GenotypeAnno>,
    ) {
        // Some variants aren't annotated, we'll fill it in here
        if annots.is_empty() {
            annots.resize_with(self.sample_count, || {
                GenotypeAnno::new(&NodeIndex::new(0), &[], 0, &Ploidy::Zero, 0, 0)
            });
        }

        // TODO: This is broken... gtcounts will need to be done per-sample
        *self.gtcounts.entry(annots[0].gt_state).or_insert(0) += 1;

        let out_fields: Vec<_> = annots.iter().map(|a| a.make_fields()).collect();
        *entry.samples_mut() = Samples::new(self.keys.clone(), out_fields);

        self.buf.clear();
        let mut tmp = vcf::io::Writer::new(&mut self.buf);
        if tmp.write_variant_record(&self.header, &entry).is_err() {
            let changed = replace_iupac_inplace(entry.reference_bases_mut());
            self.iupac_fixed |= changed;
            if let Err(error) = self.writer.write_variant_record(&self.header, &entry) {
                panic!("Couldn't write record {:?}", error);
            }
        } else if let Err(error) = self.writer.get_mut().write_all(&self.buf) {
            panic!("Couldn't write record {:?}", error);
        }
    }
}

pub fn open_writer_thread(
    result_receiver: Receiver<ChannelOutput>,
    out_path: Option<PathBuf>,
    sample_names: Vec<String>,
    wt_header: vcf::Header,
    wt_num_variants: std::sync::Arc<std::sync::Mutex<u64>>,
) -> std::thread::JoinHandle<()> {
    std::thread::spawn(move || {
        let mut m_writer = VcfWriter::new(&out_path, wt_header.clone(), &sample_names);

        let mut pbar: Option<ProgressBar> = None;
        let sty = ProgressStyle::with_template(
            " [{elapsed_precise}] {bar:44.cyan/blue} > {pos} completed",
        )
        .unwrap()
        .progress_chars("・🐷🥫");

        let mut completed_variants: u64 = 0;
        loop {
            match result_receiver.recv() {
                Ok(None) | Err(_) => {
                    pbar.expect("I actually shouldn't be expecting the bar")
                        .finish();
                    break;
                }
                Ok(Some(result)) => {
                    let mut rsize: u64 = 0;
                    for (entry, annos) in result {
                        m_writer.anno_write(entry, annos);
                        rsize += 1;
                    }

                    if let Some(ref mut bar) = pbar {
                        bar.inc(rsize);
                    } else {
                        completed_variants += rsize;
                        // check if the reader is finished so we can setup the pbar
                        let value = *wt_num_variants.lock().unwrap();
                        if value != 0 {
                            let t_bar = ProgressBar::new(value).with_style(sty.clone());
                            t_bar.inc(completed_variants);
                            pbar = Some(t_bar);
                        }
                    }
                }
            }
        }
        if m_writer.iupac_fixed {
            warn!("Some IUPAC codes in REF sequences have been fixed in output");
        }
        info!("genotype counts: {:#?}", m_writer.gtcounts);
    })
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
