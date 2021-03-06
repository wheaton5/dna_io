extern crate flate2;
extern crate rust_htslib;

use std::io::Error;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;

use std::io::BufReader;
use std::io::BufWriter;
use std::io::BufRead;
use std::io::Write;
use std::fs::File;

use rust_htslib::sam;
use rust_htslib::bam;
use rust_htslib::prelude::*;

#[derive(Debug)]
pub enum DnaFormat {
    Fastq,
    Fasta,
    Bam, 
    Sam,
    Cram, 
    TwoBit,
}
use DnaFormat::*;

#[derive(Debug,PartialEq,Clone)]
pub enum Compression {
    Gzipped,
    Uncompressed,
}
use Compression::*;

pub struct DnaRecord {
    pub seq: String,
    pub qual: Option<String>,
    pub name: String,
}

pub trait DnaRead {
    fn next(&mut self) -> Option<DnaRecord>;
    fn my_type(&self) -> DnaFormat;
    fn header(&self) -> Option<bam::Header>;
    fn extension(&self) -> String;
}

pub trait DnaWrite {
    fn write(&mut self, rec: &DnaRecord) -> Result<(), Error>;
}

pub struct DnaReader {
    pub reader: Box<DnaRead>,
}

pub struct DnaWriter {
    pub writer: Box<DnaWrite>,
}

fn check_extension(filename: &str) -> (DnaFormat, Compression) {
    let filetype = filename.split(".").collect::<Vec<&str>>();
    assert!(filetype.len() >= 2, "file {} has no extension", filename);
    match filetype[filetype.len()-1] {
        "gz" => {
            assert!(filetype.len() >= 3, "gz file {} has no format extension", filename);
            match filetype[filetype.len()-2] {
                "fa" | "fasta" => (Fasta, Gzipped),
                "fq" | "fastq" => (Fastq, Gzipped),
                _ => panic!("format of file {} not supported",filename),
            }
        },
        "fastq" | "fq" => (Fastq, Uncompressed),
        "fasta" | "fa" => (Fasta, Uncompressed),
        "sam" => (Sam, Uncompressed),
        "bam" => (Bam, Gzipped), // this isnt strictly true, can have uncompressed bam, but bam library will deal with this
        "cram" => (Cram, Gzipped), // same, also unimplemented
        "2Bit" => (TwoBit, Uncompressed), //unimplemented
        _ => panic!("format of file {} not supported ",filename),
    }
}

impl DnaReader {
    pub fn from_path(filename: &str) -> Self {
        let (file_fmt, compression) = check_extension(filename);
        let reader: Box<DnaRead> = match file_fmt {
            Fasta => Box::new(FastaReader::new(filename, compression)),
            Fastq => Box::new(FastqReader::new(filename, compression)),
            Bam => Box::new(BamReader::new(filename)),
            Sam => Box::new(SamReader::new(filename)),
            _ => panic!("file extension type {:?} not accepted.",file_fmt),
        };
        DnaReader{reader: reader}
    }
    pub fn header(&self) -> Option<bam::Header> { self.reader.header() }
    pub fn my_type(&self) -> DnaFormat { self.reader.my_type() }
    pub fn extension(&self) -> String { self.reader.extension() }
}

// for now we will assume that output is not gz'ed. they may want to stream to another program
impl DnaWriter {
    pub fn from_reader(filename: &str, reader: &DnaReader) -> Self {
        let writer: Box<DnaWrite> = match reader.my_type() {
            Fastq => Box::new(FastqWriter::new(filename, Uncompressed)),
            Fasta => Box::new(FastaWriter::new(filename, Uncompressed)),
            Sam | Bam | Cram => Box::new(SamWriter::new(filename, reader)),
            TwoBit => panic!("unimplemented"),
        };
        DnaWriter{ writer: writer }
    }
    pub fn from_path(filename: &str) -> Self {
        let (file_fmt, compression) = check_extension(filename);
        let writer: Box<DnaWrite> = match file_fmt {
            Fasta => Box::new(FastaWriter::new(filename, compression)),
            Fastq => Box::new(FastqWriter::new(filename, compression)),
            Sam => panic!("require from_reader for sam, I dont know how to make headers"),//Box::new(SamWriter(filename)),
            Bam => panic!("requires from_reader for bam writer, I dont know how to make headers"),//Box::new(BamWriter::new(filename, reader)),
            _ => panic!("file extension type {:?} not accepted.",file_fmt),
        };
        DnaWriter{ writer: writer }
    }
    pub fn write(&mut self, rec: &DnaRecord) -> Result<(), Error> { self.writer.write(rec) }
}
pub fn flush(writer: DnaWriter) {} // drop out of scope and flush/free automatically

impl Iterator for DnaReader {
    type Item = DnaRecord;
    fn next(&mut self) -> Option<DnaRecord> {
        self.reader.next()
    }
}

fn get_reader(filename: &str, compression: Compression) -> BufReader<Box<std::io::Read>> {
	let file = File::open(filename).expect("There was a problem opening the file");
    let reader: Box<std::io::Read> = match compression {
      	Gzipped => Box::new(GzDecoder::new(file)),
        Uncompressed => Box::new(file),
    };
    BufReader::new(reader)
}

fn get_writer(filename: &str, compression: Compression) -> BufWriter<Box<std::io::Write>> {
	let file = File::create(filename).expect("Unable to create file");
    let writer: Box<std::io::Write> = match compression {
        Gzipped => Box::new(GzEncoder::new(file, flate2::Compression::default())),
        Uncompressed => Box::new(file),
    };
    BufWriter::new(writer)
}

pub struct FastqReader {
    pub buf_reader: BufReader<Box<std::io::Read>>,
    compression: Compression,
}

pub struct FastqWriter {
    pub buf_writer: BufWriter<Box<std::io::Write>>,
}

impl FastqReader {
    fn new(filename: &str, compression: Compression) -> Self {
        FastqReader{ buf_reader: get_reader(filename, compression.clone()) , compression: compression}
    }
}

impl FastqWriter {
    fn new(filename: &str, compression: Compression) -> Self {
        FastqWriter{ buf_writer: get_writer(filename, compression) }
    }
}

impl DnaRead for FastqReader {
    fn next(&mut self) -> Option<DnaRecord> {     
        let mut name = String::new();
        let mut seq = String::new();
        let mut sep = String::new();
        let mut qual = String::new();
        match self.buf_reader.read_line(&mut name).expect("Could not read file") {0 => return None, _ => ()};
        name.pop();
        match self.buf_reader.read_line(&mut seq).expect("Could not read file") {0 => return None, _ => ()};
        seq.pop();
		match self.buf_reader.read_line(&mut sep).expect("Could not read file") {0 => return None, _ => ()};
		match self.buf_reader.read_line(&mut qual).expect("Could not read file") {0 => return None, _ => ()};
        qual.pop();
		Some(DnaRecord{ name: name, seq: seq, qual: Some(qual)})        
    }
	fn my_type(&self) -> DnaFormat {
		Fastq
	}
    fn header(&self) -> Option<bam::Header> {
        None
    }
    fn extension(&self) -> String {
        let mut to_ret = ".fastq".to_string();
        if self.compression == Gzipped { to_ret.push_str(".gz") }
        to_ret
    }
}

impl DnaWrite for FastqWriter {
    fn write(&mut self, rec: &DnaRecord) -> Result<(), Error> {
        let qual = match rec.qual {
            Some(ref x) => x,
            None => panic!("I have no qual i cant write fastq"),
        };
        let to_write = format!("{}\n{}\n+\n{}\n",rec.name,rec.seq,qual);
        self.buf_writer.write_all(&to_write.as_bytes())
    }
}

pub struct FastaReader {
    pub buf_reader: BufReader<Box<std::io::Read>>,
    pub last_name: Option<String>,
    compression: Compression,
}

pub struct FastaWriter {
    pub buf_writer: BufWriter<Box<std::io::Write>>,
}

impl FastaReader {
    fn new(filename: &str, compression: Compression) -> Self {
        FastaReader{ buf_reader: get_reader(filename, compression.clone()) , last_name: None, compression: compression}
    }
}

impl FastaWriter {
	fn new(filename: &str, compression: Compression) -> Self {
		FastaWriter{ buf_writer: get_writer(filename, compression) }
	}
}


impl DnaRead for FastaReader {
	fn next(&mut self) -> Option<DnaRecord> {
        let mut name = String::new();
        let mut seq = String::new();
        let mut next_name = String::new();

        match self.last_name {
            Some(ref my_name) => {
                println!("got here");
                name.push_str(my_name);
                'line_iter: loop {
                    let mut line = String::new();
                    match self.buf_reader.read_line(&mut line) {
                        Ok(_) => {
                            if line.starts_with(">") {
                                next_name.push_str(&line);
                                next_name.pop();
                                break 'line_iter;
                            }
                            else if line.is_empty() {
                                break 'line_iter; 
                            }
                            else {
                                seq.push_str(&line);
                                seq.pop();
                            }
                        },
                        Err(err) => {
                            panic!("{}",err);
                        },
                   }      
                }
            },
            None => {
                println!("first");
                let mut line = String::new();
                match self.buf_reader.read_line(&mut line) {
                    
                    Ok(_) => {
                        println!("{}",line);
                        if line.starts_with(">") {
                            name.push_str(&line);
                            name.pop();
                        } else if line.is_empty() {
                            println!("empty");
                            return None;
                        } else {
                            panic!("not fasta format?");
                        }
                    },
                    Err(err) => panic!("{}",err),
                }
                'line_iter2: loop {
                    let mut line = String::new();
                    match self.buf_reader.read_line(&mut line) {
                        Ok(x) => {
                            println!("here {} {}", line,x);
                            if line.starts_with(">") {
                                next_name.push_str(&line);
                                next_name.pop();
                                break 'line_iter2;
                            }
                            else if line.is_empty() { break 'line_iter2; }
                            else {
                                seq.push_str(&line);
                                seq.pop();
                            }
                        },
                        Err(err) => {
                            panic!("{}",err);
                        },
                    }
                }
            },
        }
        if !next_name.is_empty() {
            self.last_name = Some(next_name);
        } else {
            self.last_name = None;
        }
        Some(DnaRecord{ name: name, seq: seq, qual: None })
	}
    fn header(&self) -> Option<bam::Header> { None }
    fn my_type(&self) -> DnaFormat { Fasta }
    fn extension(&self) -> String {
        let mut to_ret = ".fasta".to_string();
        if self.compression == Gzipped { to_ret.push_str(".gz"); }
        to_ret
    }
}

impl DnaWrite for FastaWriter {
	fn write(&mut self, rec: &DnaRecord) -> Result<(), Error> {
        let to_write = format!("{}\n{}\n",rec.name, rec.seq);
        self.buf_writer.write_all(&to_write.as_bytes())
	}
}

pub struct BamReader {
    pub reader: bam::Reader,
}

pub struct BamWriter {
	pub writer: bam::Writer,
}

impl BamReader {
    fn new(filename: &str) -> Self {
        let bam = bam::Reader::from_path(filename).expect("could not open file for bam reading");
        BamReader { reader: bam }
    }
}

impl BamWriter {
	fn new(filename: &str, template: bam::Reader) -> Self {
		let header = bam::Header::from_template(template.header());
		let mut bam = bam::Writer::from_path(filename, &header).expect("could not open bam for writing");
		BamWriter{ writer: bam }
	}
}

impl DnaRead for BamReader {
    fn next(&mut self) -> Option<DnaRecord> {
        let mut record = bam::record::Record::new();
		match self.reader.read(&mut record) {
            Err(bam::ReadError::NoMoreRecord) => return None,
            Ok(_x) => (),//Some(Ok(x)),
            Err(_err) => panic!("bam error, im lazy and cant be bothered to make good error messages"),
        }
        Some(DnaRecord{ 
            name: String::from_utf8_lossy(record.qname()).to_string(), 
            qual: Some(String::from_utf8_lossy(record.qual()).to_string()),
            seq: String::from_utf8_lossy(&record.seq().as_bytes()).to_string(),
        })
            
    }
    fn my_type(&self) -> DnaFormat { Bam }
    fn header(&self) -> Option<bam::Header> { Some(bam::Header::from_template(self.reader.header())) }
    fn extension(&self) -> String { ".bam".to_string() }
}

impl DnaWrite for BamWriter {
	fn write(&mut self, rec: &DnaRecord) -> Result<(), Error> {
		panic!("unimplemented");
	}
}

pub struct SamReader {
    buf_reader: BufReader<Box<std::io::Read>>,
}

pub struct SamWriter {
	writer: sam::Writer,
}

impl SamReader {
    fn new(filename: &str) -> Self {
        let reader = get_reader(filename, Uncompressed);
        SamReader{ buf_reader: reader }
    }
}

impl SamWriter {
	fn new(filename: &str, template: &DnaReader) -> Self {
		let header = match template.header() {
            Some(x) => x,
            None => panic!("i have no header for template"),
        };
		let writer = sam::Writer::from_path(filename, &header).expect("could not open sam file for writing");
        SamWriter{ writer: writer }
	}
}

impl DnaRead for SamReader {
    fn next(&mut self) -> Option<DnaRecord> {
        let mut line = String::new();
        let mut to_ret: Option<DnaRecord> = None;
        while !self.buf_reader.read_line(&mut line).unwrap_or(0) > 0 && !line.starts_with("@") {
            let line: Vec<&str> = line.split_whitespace().collect();
            assert!( line.len() > 10, "is this sam format?, error parsing");
            to_ret = Some( DnaRecord{ name: line[0].to_string(), seq: line[9].to_string(), qual: Some(line[10].to_string()) } )
        };
        to_ret
    }
    fn my_type(&self) -> DnaFormat { Sam }
    fn header(&self) -> Option<bam::Header> { None }
    fn extension(&self) -> String { ".sam".to_string() }
}

impl DnaWrite for SamWriter {
	fn write(&mut self, rec: &DnaRecord) -> Result<(), Error> {
		panic!("unimplemented");
        let mut bam_rec = bam::Record::new();
        let qual = match rec.qual {
            Some(ref x) => x,
            None => panic!("what am i doing writing sam, i have no qual"),
        };
        bam_rec.set(&rec.name.as_bytes(), Some(&bam::record::CigarString::from_str("").unwrap()), &rec.seq.as_bytes(), &qual.as_bytes());
        match self.writer.write(&bam_rec) {
            Ok(_) => Ok(()),
            Err(_) => panic!("write error on sam write"),
        }
	}
}

mod tests {
    #[allow(unused_imports)]
    use DnaReader;
    #[allow(unused_imports)]
    use DnaRecord;
    #[allow(unused_imports)]
    use DnaWriter;
    use DnaWrite;
    use std::io::Read;
    use std::fs::File;
    use flush;

    #[test]
    fn test_fastq() {
        let mut reader = DnaReader::from_path("test/data/fastq.fastq");
        let rec = match reader.next() {
            Some(x) => x,
            None => panic!("no records"),
        };
        println!("{}",rec.seq);
        assert!("ACTGGTCA"==rec.seq);
        reader.next();
        match reader.next() {
            Some(_) => panic!("should not be any more records"),
            None => (),
        }

        let mut reader = DnaReader::from_path("test/data/fastq.fq");
        let rec = match reader.next() {
            Some(x) => x,
            None => panic!("no records"),
        };
        println!("{}",rec.seq);
        assert!("ACTGGTCA"==rec.seq);

        let mut reader = DnaReader::from_path("test/data/fastq.fastq.gz");
        let rec = match reader.next() {
            Some(x) => x,
            None => panic!("no records"),
        };
        println!("{}",rec.seq);
        assert!("ACTGGTCA"==rec.seq);

        let mut reader = DnaReader::from_path("test/data/fastq.fq.gz");
        let rec = match reader.next() {
            Some(x) => x,
            None => panic!("no records"),
        };
        println!("{}",rec.seq);
        assert!("ACTGGTCA"==rec.seq); 
    }

    #[test]
    fn test_fasta() {
        let mut reader = DnaReader::from_path("test/data/fasta.fasta");
        println!("about to");
        let rec = match reader.next() {
            Some(x) => x,
            None => panic!("no records"),
        };
        println!("{}",rec.seq);
        assert!("ACGTTTTTTTTTTTTTTACGT" == rec.seq);
        println!("done with firts");
        let rec = match reader.next() {
            Some(x) => x,
            None => panic!("no second record"),
        };
        println!("{}",rec.seq);
        assert!("GGGGGGGGGGGGGGGGGGGG" == rec.seq);
        match reader.next() {
            Some(_) => panic!("there shouldnt be any more records"),
            None => (),
        }
    }

    #[test]
    fn test_bam() {
        let mut reader = DnaReader::from_path("test/data/test.bam");
        let rec = match reader.next() {
            Some(x) => x,
            None => panic!("bam reader doesnt work"),
        };
        println!("{}",rec.seq);
        assert!("GTCCTAAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCTAACCATACCCATAACCCCAACCCTAACACTAACCCCAAACCCAACCATAACCAACACCCCACACCTA" == rec.seq);
    }

    #[test]
    fn test_write_fastq() {
        let mut reader = DnaReader::from_path("test/data/fastq.fastq");
        let mut writer = DnaWriter::from_reader("test/data/fastq_written.fastq",&reader);
        for rec in reader {
            writer.write(&rec).expect("failed to write fastq file in test");
        }
        flush(writer);
        let mut file = File::open("test/data/fastq.fastq").expect("test data not available");
        let mut contents = String::new();
        file.read_to_string(&mut contents).expect("cant read test data");
        let mut file2 = File::open("test/data/fastq_written.fastq").expect("written test data not available");
        let mut contents2 = String::new();
        file2.read_to_string(&mut contents2).expect("cant read written test data");
        assert!(contents == contents2);
    }
    
    #[test]
    fn test_write_fasta() {
        let mut reader = DnaReader::from_path("test/data/fasta.fasta");
        let mut writer = DnaWriter::from_reader("test/data/fasta_written.fasta",&reader);
        for rec in reader {
            writer.write(&rec).expect("failed to write fastq file in test");
        }
        flush(writer);
        let reader = DnaReader::from_path("test/data/fasta.fasta");
        let reader2 = DnaReader::from_path("test/data/fasta_written.fasta");
        for (rec1, rec2) in reader.zip(reader2) {
            assert!(rec1.seq == rec2.seq);
            assert!(rec1.name == rec2.name);
        }
    }

    #[test]
    fn test_big_fasta() {
        println!("fork me");
        let mut reader = DnaReader::from_path("blah.fasta");
        let mut writer = DnaWriter::from_reader("out2.fasta", &reader);
        for rec in reader {
            writer.write(&rec).expect("nope");
        }
        assert!(2==3);
        
    }

    


}
