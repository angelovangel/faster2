use std::process;
use kseq::parse_path;

extern crate clap;
use clap::{App, Arg, ArgGroup};

use faster2;


fn main(){

    let matches = App::new("faster2")
        .version("0.1.0")
        .author("https://github.com/angelovangel")
        .about("fast statistics and manipulation of a fastx file")
        
        .arg(Arg::with_name("len")
            .long("len")
            .short("l")
            //.takes_value(true)
            .required(false)
            .help("Output read lengths"))

        .arg(Arg::with_name("gc")
            .long("gc")
            .short("g")
            .required(false)
            .help("output GC content per read"))
        
        .arg(Arg::with_name("qual")
            .long("qual")
            .short("q")
            .required(false)
            .takes_value(false)
            .help("Output 'average' q-score per read"))
        
        .arg(Arg::with_name("table")
            .long("table")
            .short("t")
            .takes_value(false)
            .required(false))
        
        .arg(Arg::with_name("INPUT")
            .help("Path to a fastq file")
            .required(true)
            .index(1))
        .group(ArgGroup::with_name("group").required(true).args(&["table", "len", "qual", "gc"]))
        .get_matches();
	
    let infile = matches.value_of("INPUT").unwrap().to_string();
	let mut records = parse_path(infile).unwrap();
	
    if matches.is_present("len") {
        while let Some(record) = records.iter_record().unwrap() {
            println!( "{}", record.len() );
        }
        process::exit(0);

    } else if matches.is_present("qual") {
        while let Some(record) = records.iter_record().unwrap() {
            //let mean_prob = faster2::qscore_probs( record.qual().as_bytes() ) / record.len() as f32;
            //println!("{:.4}", -10.0 * mean_prob.log10()); // convert to phred score
            let mean_prob = faster2::phred_gm( record.qual().as_bytes() );
            println!("{:.4}", mean_prob);
        }
        process::exit(0);

    } else if matches.is_present("gc") {
        while let Some(record) = records.iter_record().unwrap() {
            let gc_content = faster2::get_gc_content(record.seq().as_bytes() );
            println!( "{:.4}", gc_content )
        }
        process::exit(0);
    
    } else if matches.is_present("table") {
        let mut reads: i64 = 0;
        let mut bases: i64 = 0;
        let mut num_n: i64 = 0;
        let mut qual20: i64 = 0;
        //let mut qual30: i64 = 0;
        let mut len_vector: Vec<i64> = Vec::new();

        while let Some(record) = records.iter_record().unwrap() {
            let len = record.len() as i64;
            len_vector.push(len);
            reads += 1;
            bases += record.len() as i64;
            num_n += faster2::get_n_bases(record.seq().as_bytes() );
            qual20 += faster2::get_qual_bases(record.qual().as_bytes(), 53); // 33 offset
            //qual30 += get_qual_bases(record.qual().as_bytes(), 63);
        
        }

        let q20 = qual20 as f64 / bases as f64 * 100.0;
        //let q30 = qual30 as f64 / bases as f64 * 100.0;
        let n50 = faster2::get_nx(&mut len_vector, 0.5);
        let min_len = len_vector.iter().min().unwrap();
        let max_len = len_vector.iter().max().unwrap();

        println!("reads\tbases\tn_bases\tmin_len\tmax_len\tN50\tQ20_percent");
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{:.2}", reads, bases, num_n, min_len, max_len, n50, q20);
    }
}