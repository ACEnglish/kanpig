use crate::kplib::{
    germ_genotyper::{GTstate, Genotyper},
    ChannelOutput, GenotypeAnno, Ploidy,
};
use crossbeam_channel::Receiver;
use indicatif::{ProgressBar, ProgressStyle};
use petgraph::graph::NodeIndex;
use std::{
    collections::HashMap,
    env,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};

use noodles_vcf::{
    self as vcf,
    header::record::{
        value::{map::format, Map},
        Value,
    },
    variant::io::Write as vcfWrite,
    variant::record::Ids,
    variant::record_buf::samples::{keys::Keys, Samples},
};

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

pub struct VcfWriter {
    writer: vcf::io::Writer<Box<dyn Write>>,
    sample_count: usize,
    header: vcf::Header,
    keys: Keys,
    pub gtcounts: Vec<HashMap<GTstate, usize>>,
    genotyper: Genotyper,
    rnames_writer: Option<Box<dyn Write>>,
}

impl VcfWriter {
    /// Given a path and a header, setup a new output VCF
    pub fn new(
        out_path: &Option<PathBuf>,
        rnames_path: &Option<PathBuf>,
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
        // Comment line for version/params
        let command = env::args().collect::<Vec<String>>().join(" ");
        let comment = format!("<version='v{}',command='{}'>", VERSION, command);
        let _ = header.insert(PKG_NAME.parse().expect("const"), Value::String(comment));

        // Setup FORMAT header definitions
        let all_formats = header.formats_mut();
        all_formats.clear();
        let format_definitions = GenotypeAnno::make_format();
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

        let genotyper = Genotyper::with_optional_config(None);

        let rnames_writer: Option<Box<dyn Write>> = match rnames_path {
            Some(ref path) => {
                let file = File::create(path).expect("Error creating output file");
                let m_page = page_size::get() * 1000;
                Some(Box::new(BufWriter::with_capacity(m_page, file)))
            }
            None => None,
        };

        Self {
            writer,
            header,
            sample_count: sample_names.len(),
            keys: Keys::from_iter(new_fmts),
            gtcounts: vec![HashMap::new(); sample_names.len()],
            genotyper,
            rnames_writer,
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
                GenotypeAnno::new(
                    &NodeIndex::new(0),
                    &[],
                    0,
                    &Ploidy::Zero,
                    0,
                    0,
                    &self.genotyper,
                )
            });
        }

        for (gtcount_map, annot) in self.gtcounts.iter_mut().zip(annots.iter()) {
            *gtcount_map.entry(annot.gt_state).or_insert(0) += 1;
        }

        let out_fields: Vec<_> = annots.iter().map(|a| a.make_fields()).collect();
        *entry.samples_mut() = Samples::new(self.keys.clone(), out_fields);

        if let Err(error) = self.writer.write_variant_record(&self.header, &entry) {
            panic!("Couldn't write record {:?}", error);
        }

        // rnames writing here
        if let Some(writer) = &mut self.rnames_writer {
            let m_id = entry.ids().iter().collect::<Vec<_>>().join(";");
            for annot in annots {
                for rname in &annot.rnames {
                    let _ = writeln!(writer, "{}\t{}", m_id, rname);
                }
            }
        }
    }
}

pub fn open_writer_thread(
    result_receiver: Receiver<ChannelOutput>,
    out_path: Option<PathBuf>,
    rnames_path: Option<PathBuf>,
    sample_names: Vec<String>,
    wt_header: vcf::Header,
    wt_num_variants: std::sync::Arc<std::sync::Mutex<u64>>,
) -> std::thread::JoinHandle<()> {
    std::thread::spawn(move || {
        let mut m_writer =
            VcfWriter::new(&out_path, &rnames_path, wt_header.clone(), &sample_names);

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
