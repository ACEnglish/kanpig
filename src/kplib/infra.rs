use crate::kplib::GenotypeAnno;
use noodles_vcf::variant::RecordBuf;
pub type ChannelInput = Option<Vec<RecordBuf>>;
pub type ChannelOutput = Option<Vec<(RecordBuf, Vec<GenotypeAnno>)>>;
