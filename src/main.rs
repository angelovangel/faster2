use std::{path::Path, process};
use kseq::parse_path;

extern crate clap;
use clap::{App, Arg, ArgGroup};
use indicatif::{HumanCount, ProgressBar};
use human_repr::HumanThroughput;
use prettytable::Slice;
use std::time::{Duration, Instant};
use colored::Colorize;
use serde_json::json;

use faster2;
use noodles_bam as bam;
use noodles_sam::alignment::record::QualityScores as _;

#[macro_use] extern crate prettytable;

struct SimpleRecord {
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl SimpleRecord {
    fn len(&self) -> usize {
        self.seq.len()
    }
    fn seq(&self) -> &[u8] {
        &self.seq
    }
    fn qual(&self) -> &[u8] {
        &self.qual
    }
}

trait ToBytes {
    fn as_bytes(&self) -> &[u8];
}

impl ToBytes for &[u8] {
    fn as_bytes(&self) -> &[u8] {
        self
    }
}

impl ToBytes for &Vec<u8> {
    fn as_bytes(&self) -> &[u8] {
        self.as_slice()
    }
}

impl ToBytes for String {
    fn as_bytes(&self) -> &[u8] {
        self.as_bytes()
    }
}

macro_rules! run_analysis {
    ($source:ident, $next_expr:expr, $infile:expr, $matches:expr, $start:expr) => {
        let inf = $infile;
        let mtch = $matches;
        let st = $start;
        
        if mtch.is_present("len") {
            while let Some(record) = $next_expr {
                println!( "{}", record.len() );
            }
            process::exit(0);

        } else if mtch.is_present("qual") {
            while let Some(record) = $next_expr {
                let mean_prob = faster2::qscore_probs( record.qual().as_bytes() ) / record.len() as f32;
                println!("{:.4}", -10.0 * mean_prob.log10()); 
            }
            process::exit(0);
        
        } else if mtch.is_present("gc") {
            while let Some(record) = $next_expr {
                let gc_content = faster2::get_gc_content(record.seq().as_bytes() );
                println!( "{:.4}", gc_content )
            }
            process::exit(0);
        } else if mtch.is_present("nx") {
            let nxvalue = mtch
                .value_of("nx")
                .unwrap()
                .parse::<f32>()
                .expect("Failed to parse nx value");
            match nxvalue {
                x if (0.0..=1.0).contains(&x) => {
                    let mut lengths: Vec<u64> = Vec::new();
                    while let Some(record) = $next_expr {
                        let len = record.len() as u64;
                        lengths.push(len);
                    }
                    let nx = faster2::get_nx(&mut lengths, 1.0 - nxvalue);
                    println!("{}", nx);
                }
                _=> {
                    eprintln!("The NX value should be between 0.0 and 1.0"); 
                    process::exit(0)
                }
            }

        } else if mtch.is_present("qyield") {
            let qvalue = mtch
                .value_of("qyield")
                .unwrap().parse::<u8>()
                .expect("Failed to parse Q score value");
            match qvalue {
                x if (1..=93).contains(&x) => {
                    let mut bases: i64 = 0;
                    let mut qualx: i64 = 0;
                    while let Some(record) = $next_expr {
                        bases += record.len() as i64;
                        qualx += faster2::get_qual_bases(record.qual().as_bytes(), qvalue + 33); 
                    }
                    let qx = qualx as f64 / bases as f64 * 100.0;
                    println!("Q{}\t{:.2}", qvalue, qx);
                    process::exit(0)
                }
                _=> {
                    eprintln!("The Q score value should be between 1 and 93"); 
                    process::exit(0)
                }
            }
    

        } else if mtch.is_present("table") {
            let filename = Path::new(inf).file_name().unwrap().to_str().unwrap();
            let mut reads: u64 = 0;
            let mut bases: u64 = 0;
            let mut num_n: u64 = 0;
            let mut qual20: i64 = 0;
            let mut gc_bases: u64 = 0;
            let mut len_vector: Vec<u64> = Vec::new();
            let pb = ProgressBar::new_spinner();
            pb.enable_steady_tick(Duration::from_millis(120));

            while let Some(record) = $next_expr {
                let len = record.len() as u64;
                len_vector.push(len);

                reads += 1;
                bases += len;
                gc_bases += faster2::get_gc_bases(record.seq().as_bytes());

                num_n += faster2::get_n_bases(record.seq().as_bytes() );
                qual20 += faster2::get_qual_bases(record.qual().as_bytes(), 53); 
                
                let speed = (reads as u128).checked_div(st.elapsed().as_millis()).ok_or(1);
                if reads % 10000 == 0 {
                    let message = format!(
                        "Processed reads: {} ({})", 
                        HumanCount(reads as u64).to_string().green(),
                        (speed.unwrap_or(0)*1000).human_throughput(" reads").to_string().red()
                    );
                    pb.set_message(message);
                }
            
            }

            let q20 = if bases > 0 { qual20 as f64 / bases as f64 * 100.0 } else { 0.0 };
            let n50 = faster2::get_nx(&mut len_vector, 0.5);
            let min_len = len_vector.iter().min().unwrap_or(&0);
            let max_len = len_vector.iter().max().unwrap_or(&0);
            let gc_content = if bases > 0 { gc_bases as f64/ bases as f64 * 100.0 } else { 0.0 };
    
            pb.finish_and_clear();

            if mtch.is_present("pretty") {
                let table = table!( 
                    ["file", "reads", "bases", "nbases", "min_len", "max_len", "N50", "GC_percent", "Q20_percent"],
                    [filename, HumanCount(reads), HumanCount(bases), 
                    HumanCount(num_n), HumanCount(*min_len), HumanCount(*max_len), HumanCount(n50), 
                    format!("{:.2}", gc_content), format!("{:.2}", q20)]);
                if mtch.is_present("skip_header") {
                    let slice = table.slice(1..);
                    slice.printstd();
                } else {
                    table.printstd();            
                }


            } else if mtch.is_present("json") {
                let json = json!({
                    "file": filename, "reads": reads, "bases": bases, "num_n": num_n, 
                    "min_len": min_len, "max_len": max_len, "n50": n50, "gc_percent": gc_content, "q20": q20
                });
                println!("{}", json.to_string())
            
            } else {
                if !mtch.is_present("skip_header") {
                    println!("file\treads\tbases\tn_bases\tmin_len\tmax_len\tN50\tGC_percent\tQ20_percent");
                }
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}", 
                    filename, reads, bases, num_n, min_len, max_len, n50, gc_content, q20
                );

            }
        }
    }
}

fn main(){

    let matches = App::new("faster2")
        .version("0.3.0")
        .author("https://github.com/angelovangel")
        .about("fast statistics of a fastx file")
        
        .arg(Arg::with_name("len")
            .long("len")
            .short('l')
            //.takes_value(true)
            .required(false)
            .help("Output read lengths"))

        .arg(Arg::with_name("gc")
            .long("gc")
            .short('g')
            .required(false)
            .help("output GC content per read"))
        
        .arg(Arg::with_name("qual")
            .long("qual")
            .short('q')
            .required(false)
            .takes_value(false)
            .help("Output 'average' q-score per read"))

        .arg(Arg::with_name("nx")
            .long("nx")
            .short('x')
            .takes_value(true)
            .help("Output NX value, provide the desired NX value as 0.5 for e.g. N50"))
        .arg(Arg::with_name("qyield")
            .long("qyield")
            .short('y')
            .takes_value(true)
            .help("Percent bases with Q score of x or higher (use range 1..93)"))
        .arg(Arg::with_name("table")
            .long("table")
            .short('t')
            .takes_value(false)
            .required(false)
            .help("Output table summary"))
        .arg(Arg::with_name("pretty")
            .long("pretty")
            .short('p')
            .takes_value(false)
            .required(false)
            .help("pretty print table summary, works only together with the --table flag"))
        .arg(Arg::with_name("skip_header")
            .long("skip_header")
            .short('s')
            .takes_value(false)
            .required(false)
            .help("Skip header in summary table, works only together with the --table flag"))
        .arg(Arg::with_name("json")
            .long("json")
            .short('j')
            .takes_value(false)
            .required(false)
            .help("output summary data as json, works only together with the --table flag"))
        .arg(Arg::with_name("INPUT")
            .help("Path to a fastq file")
            .required(true)
            .multiple_values(true)
            .min_values(1)
            .index(1))
        .group(ArgGroup::with_name("group").required(true).args(&["table", "len", "qual", "gc", "nx", "qyield"]))
        .get_matches();

    let infiles: Vec<&str> = matches.values_of("INPUT").unwrap().collect();
    for infile in infiles {
        let start = Instant::now();
        
        if infile.ends_with(".bam") {
            let mut reader = bam::io::reader::Builder::default()
                .build_from_path(infile)
                .expect("Failed to open BAM file");
            let _header = reader.read_header().expect("Failed to read BAM header");
            
            let mut records = reader.records().map(|r| {
                let record = r.expect("Failed to read BAM record");
                SimpleRecord {
                    seq: record.sequence().iter().map(|b| u8::from(b)).collect(),
                    qual: record.quality_scores().iter().map(|q| u8::from(q) + 33).collect(),
                }
            });
            run_analysis!(records, records.next(), infile, &matches, start);
        } else {
            let mut loader = parse_path(infile).unwrap();
            run_analysis!(loader, loader.iter_record().unwrap(), infile, &matches, start);
        }
    }
}