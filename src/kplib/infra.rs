use crate::kplib::GenotypeAnno;
use noodles_vcf::variant::RecordBuf;
pub type ChannelInput = Option<Vec<RecordBuf>>;
pub type ChannelOutput = Option<Vec<(RecordBuf, Vec<GenotypeAnno>)>>;
pub type KmerVec = Vec<(u64, f32)>;
