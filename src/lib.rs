extern crate flate2;
use flate2::read::GzDecoder;

use std::io::BufReader;
use std::io::BufRead;
use std::fs::File;

use std::io;
use std::io::prelude::*;


pub struct DnaRecord {
    pub seq: String,
    pub qual: Option<String>,
    pub name: String,
}

pub trait DnaRead {
    fn next(&mut self) -> Option<DnaRecord>;
}

pub struct DnaReader {
    pub reader: Box<DnaRead>,
}

impl DnaReader {
    fn from_path(filename: &str) -> Self {
        let filetype = filename.split(".").collect::<Vec<&str>>();
        if filetype.len() < 2 {
            panic!("file {} has no extension");
        }
        let reader: Box<DnaRead> = match filetype[filetype.len()-1] {
            "gz" => match filetype[filetype.len()-2] {
                "fa" | "fasta" => Box::new(FastaReader::new(filename)),
                "fq" | "fastq" => Box::new(FastqReader::new(filename)),
                _ => panic!("file extension {}.{} not accepted",filetype[filetype.len()-1],filetype[filetype.len()-2])
            }
            "fastq" | "fq" => Box::new(FastqReader::new(filename)),
            "fasta" | "fa" => Box::new(FastaReader::new(filename)),
            "bam" => Box::new(BamReader::new(filename)),
            "sam" => Box::new(SamReader::new(filename)),
            _ => panic!("file extension type {} not accepted.",filetype[filetype.len()-1])
        };
        DnaReader{reader: reader}
    }
    fn next(&mut self) -> Option<DnaRecord> {
        self.reader.next()
    }
}

pub struct FastqReader {
    pub buf_reader: BufReader<Box<Read>>,
}

impl FastqReader {
    fn new(filename: &str) -> Self {
        let file = match File::open(filename.to_string()) {
            Ok(file) => file,
            Err(error) => {
                panic!("There was a problem opening the file {} with error {}", filename, error)
            },
        };
        let filetype = filename.split(".").collect::<Vec<&str>>();
        let filetype = filetype[filetype.len()-1];
        let reader: Box<Read> = match filetype {
            "gz" => Box::new(GzDecoder::new(file)),
            _ => Box::new(file),
        };
        let reader = BufReader::new(reader);
        FastqReader{ buf_reader: reader }
    }
}

impl DnaRead for FastqReader {
    fn next(&mut self) -> Option<DnaRecord> {     
        let mut name = String::new();
        let mut seq = String::new();
        let mut sep = String::new();
        let mut qual = String::new();
        match self.buf_reader.read_line(&mut name) {
            Ok(_) => (),
            Err(_) => return None,
        }
        name.pop();
        match self.buf_reader.read_line(&mut seq) {
       		Ok(_) => (),
			Err(_) => return None, 
        }
        seq.pop();
		match self.buf_reader.read_line(&mut sep) {
			Ok(_) => (),
			Err(_) => return None,
		}

		match self.buf_reader.read_line(&mut qual) {
			Ok(_) => (),
			Err(_) => return None,
		}
        qual.pop();
		Some(DnaRecord{ name: name, seq: seq, qual: Some(qual)})
        
    }
}

pub struct FastaReader {
    pub buf_reader: BufReader<Box<Read>>,
    pub last_name: Option<String>,
}

impl FastaReader {
    fn new(filename: &str) -> Self {
        let file = match File::open(filename.to_string()) {
            Ok(file) => file,
            Err(error) => {
                panic!("There was a problem opening the file {} with error {}", filename, error)
            },
        };
        let filetype = filename.split(".").collect::<Vec<&str>>();
        let filetype = filetype[filetype.len()-1];
        let reader: Box<Read> = match filetype {
            "gz" => Box::new(GzDecoder::new(file)),
            _ => Box::new(file),
        };
        let reader = BufReader::new(reader);
        FastaReader{ buf_reader: reader , last_name: None}
    }
}

impl DnaRead for FastaReader {
	fn next(&mut self) -> Option<DnaRecord> {
        let mut name = String::new();
        let mut seq = String::new();
        let mut next_name = String::new();

        match self.last_name {
            Some(ref my_name) => {
                name.push_str(my_name);
                'line_iter: loop {
                    let mut line = String::new();
                    match self.buf_reader.read_line(&mut line) {
                        Ok(_) => {
                            if line.starts_with(">") {
                                next_name.push_str(&line);
                                next_name.pop();
                                println!("trying to break");
                                break 'line_iter;
                            }
                            else if line.is_empty() {
                                break 'line_iter; 
                            }
                            else {
                                println!("pushing {}",line);
                                seq.push_str(&line);
                                seq.pop();
                            }
                        },
                        Err(_) => {
                            println!("trying to return");
                            return None;
                        },
                   }      
                }
            },
            None => {
                let mut line = String::new();
                match self.buf_reader.read_line(&mut line) {
                    Ok(_) => {
                        if line.starts_with(">") {
                            name.push_str(&line);
                            name.pop();
                        } else {
                            panic!("not fasta format?");
                        }
                    },
                    Err(_) => return None,
                }
                'line_iter: loop {
                    let mut line = String::new();
                    match self.buf_reader.read_line(&mut line) {
                        Ok(_) => {
                            if line.starts_with(">") {
                                next_name.push_str(&line);
                                next_name.pop();
                                println!("trying to break2");
                                break 'line_iter;
                            }
                            else if line.is_empty() { break 'line_iter; }
                            else {
                                println!("pushing2 {}",line);
                                seq.push_str(&line);
                                seq.pop();
                            }
                        },
                        Err(_) => {
                            println!("trying to return 2");
                            return None;
                        },
                    }
                }
            },
        }
        if !next_name.is_empty() {
            self.last_name = Some(next_name);
        }
        Some(DnaRecord{ name: name, seq: seq, qual: None })
	}
}

pub struct BamReader {
	//unimplemented
}

impl BamReader {
	//unimplemented
    fn new(filename: &str) -> Self {
        BamReader {}
    }
}

impl DnaRead for BamReader {
    fn next(&mut self) -> Option<DnaRecord> {
        //unimplemented
        None
    }
}

pub struct SamReader {
	//unimplemented
}

impl SamReader {
    //unimplemented
    fn new(filename: &str) -> Self {
        SamReader {}
    }
}

impl DnaRead for SamReader {
    fn next(&mut self) -> Option<DnaRecord> {
        //unimplemented
        None
    }
}

mod tests {
    use DnaReader;
    use DnaRecord;

    #[test]
    fn test_fastq() {
        let mut reader = DnaReader::from_path("test/data/fastq.fastq");
        let rec = match reader.next() {
            Some(x) => x,
            None => panic!("no records"),
        };
        println!("{}",rec.seq);
        assert!("ACTGGTCA"==rec.seq);

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
    }
}
