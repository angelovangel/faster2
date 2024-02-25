use std::process;
use kseq::parse_path;

extern crate clap;
use clap::{App, Arg, ArgGroup};
use indicatif::{HumanCount, HumanDuration, ProgressBar};
use human_repr::HumanThroughput;
use std::time::{Duration, Instant};
use owo_colors::OwoColorize;

use faster2;

#[macro_use] extern crate prettytable;

fn main(){

    let matches = App::new("faster2")
        .version("0.3.0")
        .author("https://github.com/angelovangel")
        .about("fast statistics and manipulation of a fastx file")
        
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
            .help("Output NX value, provide the desired NX value as 0.5 for e.g. N50 [numeric]"))
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
        .arg(Arg::with_name("INPUT")
            .help("Path to a fastq file")
            .required(true)
            .index(1))
        .group(ArgGroup::with_name("group").required(true).args(&["table", "len", "qual", "gc", "nx", "qyield"]))
        .get_matches();

	let start = Instant::now();
    let infile = matches.value_of("INPUT").unwrap();
	let mut records = parse_path(infile).unwrap();
	
    if matches.is_present("len") {
        while let Some(record) = records.iter_record().unwrap() {
            println!( "{}", record.len() );
        }
        process::exit(0);

    } else if matches.is_present("qual") {
        while let Some(record) = records.iter_record().unwrap() {
            let mean_prob = faster2::qscore_probs( record.qual().as_bytes() ) / record.len() as f32;
            println!("{:.4}", -10.0 * mean_prob.log10()); // convert to phred score
            
            //let mean_prob = faster2::phred_gm( record.qual().as_bytes() );
            //println!("{:.4}", mean_prob);
        }
        process::exit(0);

    } else if matches.is_present("gc") {
        while let Some(record) = records.iter_record().unwrap() {
            let gc_content = faster2::get_gc_content(record.seq().as_bytes() );
            println!( "{:.4}", gc_content )
        }
        process::exit(0);
    } else if matches.is_present("nx") {
        let nxvalue = matches
            .value_of("nx")
            .unwrap()
            .parse::<f32>()
            .expect("Failed to parse nx value");
        match nxvalue {
            x if (0.0..=1.0).contains(&x) => {
                let mut lengths: Vec<i64> = Vec::new();
                while let Some(record) = records.iter_record().unwrap() {
                    let len = record.len() as i64;
                    lengths.push(len);
                }
                let nx = faster2::get_nx(&mut lengths, 1.0 - nxvalue);
                
                println!("{}", nx);
                //println!("N{}\t{}", nx_value as i16, nx);
            }
            _=> {
                eprintln!("The NX value should be between 0.0 and 1.0"); 
                process::exit(0)
            }
        }

    } else if matches.is_present("qyield") {
        let qvalue = matches
            .value_of("qyield")
            .unwrap().parse::<u8>()
            .expect("Failed to parse Q score value");
        match qvalue {
            x if (1..=93).contains(&x) => {
                let mut bases: i64 = 0;
                let mut qualx: i64 = 0;
                while let Some(record) = records.iter_record().unwrap() {
                    bases += record.len() as i64;
                    qualx += faster2::get_qual_bases(record.qual().as_bytes(), qvalue + 33); // 33 offset   
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
    

    } else if matches.is_present("table") {
        let mut reads: i64 = 0;
        let mut bases: i64 = 0;
        let mut num_n: i64 = 0;
        let mut qual20: i64 = 0;
        //let mut qual30: i64 = 0;
        let mut len_vector: Vec<i64> = Vec::new();
        let pb = ProgressBar::new_spinner();
        pb.enable_steady_tick(Duration::from_millis(120));

        while let Some(record) = records.iter_record().unwrap() {
            let len = record.len() as i64;
            len_vector.push(len);
            reads += 1;
            bases += record.len() as i64;
            num_n += faster2::get_n_bases(record.seq().as_bytes() );
            qual20 += faster2::get_qual_bases(record.qual().as_bytes(), 53); // 33 offset
            //qual30 += get_qual_bases(record.qual().as_bytes(), 63);
            let speed = (reads as u128).checked_div(start.elapsed().as_millis()).ok_or(1);
            
            let message = format!(
                "Processed reads: {} ({})", 
                HumanCount(reads as u64).to_string().green(),
                (speed.unwrap_or(0)*1000).human_throughput(" reads")
                //HumanDuration(start.elapsed()) 
                //speed.unwrap_or(0)
            );
            pb.set_message(message);
        
        }

        let q20 = qual20 as f64 / bases as f64 * 100.0;
        //let q30 = qual30 as f64 / bases as f64 * 100.0;
        let n50 = faster2::get_nx(&mut len_vector, 0.5);
        let min_len = len_vector.iter().min().unwrap();
        let max_len = len_vector.iter().max().unwrap();
        pb.finish_and_clear();

        if matches.is_present("pretty") {
            let table = table!( 
                ["reads", "bases", "nbases", "min_len", "max_len", "N50", "Q20_percent"],
                [reads, bases, num_n, min_len, max_len, n50, format!("{:.2}", q20)]);
            table.printstd();
        } else {
            println!("reads\tbases\tn_bases\tmin_len\tmax_len\tN50\tQ20_percent");
            println!("{}\t{}\t{}\t{}\t{}\t{}\t{:.2}", reads, bases, num_n, min_len, max_len, n50, q20);
        }
    }
}