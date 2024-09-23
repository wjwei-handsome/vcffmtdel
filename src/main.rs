// receive two position arguments: 1. vcf path 2. fasta path
// output: a vcf file, but <DEL> symbol in REF column are replaced by the reference sequence in the fasta file

use noodles::{
    fasta,
    vcf::{
        self,
        variant::{io::Write, record::info::field::key, record_buf::AlternateBases},
    },
};
use std::{
    env,
    io::{self, BufWriter},
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input_vcf = env::args().nth(1).expect("missing input vcf");
    let mut vcf_reader = vcf::io::reader::Builder::default().build_from_path(input_vcf)?;
    let header = vcf_reader.read_header()?;

    let input_fasta = env::args().nth(2).expect("missing input fasta");
    let mut fasta_reader =
        fasta::io::indexed_reader::Builder::default().build_from_path(input_fasta)?;

    // let region = raw_region.parse()?;
    // let record = reader.query(&region)?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(BufWriter::new(stdout));

    writer.write_header(&header)?;

    for result in vcf_reader.record_bufs(&header) {
        let mut record = result?;
        let chro = record.reference_sequence_name();
        let info = record.info();
        let svtype_key = key::SV_TYPE;
        let svend_key = key::END_POSITION;
        let svtype = info
            .get(svtype_key)
            // .transpose()?
            .expect("missing SVTYPE field")
            .expect("missing SVTYPE value");
        let svtype = match svtype {
            vcf::variant::record_buf::info::field::Value::String(s) => s,
            _ => panic!("expected SVTYPE to be a string"),
        };
        if svtype == "DEL" {
            let del_start = record
                .variant_start()
                // .transpose()?
                .map(usize::from)
                .unwrap_or_default();
            let del_end = info
                .get(svend_key)
                // .transpose()?
                .expect("miss end field")
                .expect("miss end value");
            let del_end = match del_end {
                vcf::variant::record_buf::info::field::Value::Integer(n) => n,
                _ => panic!("expected END to be an integer"),
            };
            let del_end = *del_end as usize;
            let raw_region = format!("{}:{}-{}", chro, del_start, del_end);
            let query_region = raw_region.parse()?;
            let query_record = fasta_reader.query(&query_region)?;
            let query_seq = query_record.sequence().as_ref();
            let query_seq_str = std::str::from_utf8(query_seq)?;
            let alt_seq = &query_seq_str[0..1];
            // let mut record_buf = RecordBuf::try_from_variant_record(&header, &record)?;
            *record.reference_bases_mut() = String::from(query_seq_str);
            *record.alternate_bases_mut() = AlternateBases::from(vec![String::from(alt_seq)]);
        }
        let _ = writer.write_variant_record(&header, &record);

        // writer.write_record(&header, &record)?;
    }
    Ok(())
}
