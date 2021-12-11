use std::process;
use kseq::parse_path;

extern crate clap;
use clap::{App, Arg, ArgGroup};

fn get_n_bases(seq: &[u8]) -> i32 {
    let mut n = 0;
    for s in seq {
        if s == &78u8 || s == &110u8 {
            n += 1;
        }
    }
    n
}

fn get_qual_bases(q: &[u8], qx: u8) -> i64 {
    let mut n = 0;
    for item in q {
        if *item >= qx {
            n += 1
        }
    }
    n
}

// to get mean of q scores from a record - first convert to prob, calc mean, then back to phred
// this fn reads phred and converts to probs and returns their sum
//
fn qscore_probs(q: &[u8]) -> f32 {
    let mut qprob_sum = 0.0;
    for &item in q.iter() {
        let phred = *&item as f32 - 33.0;
        let prob = 10.0_f32.powf(-phred / 10.0);
        qprob_sum += prob
    }
    qprob_sum
}

fn get_nx(numbers: &mut [i64], fraction: f32) -> i64 {
    numbers.sort_unstable();

    // half of the bases
    let halfsum = numbers.iter().sum::<i64>() as f32 * fraction; // f32 * f32

    // cumsum of the sorted vector
    let cumsum = numbers
        .iter()
        .scan(0, |sum, i| {
            *sum += i;
            Some(*sum)
        })
        .collect::<Vec<_>>();
    let n50_index = cumsum.iter().position(|&x| x > halfsum as i64).unwrap();

    numbers[n50_index]
}

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
        .group(ArgGroup::with_name("group").required(true).args(&["table", "len", "qual"]))
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
            let qscore = qscore_probs( record.qual().as_bytes() ) / record.len() as f32;
            println!("{:.4}", -10.0 * qscore.log10());
        }
        process::exit(0);

    } else if matches.is_present("table") {
        let mut reads: i64 = 0;
        let mut bases: i64 = 0;
        let mut num_n: i32 = 0;
        let mut qual20: i64 = 0;
        //let mut qual30: i64 = 0;
        let mut len_vector: Vec<i64> = Vec::new();

        while let Some(record) = records.iter_record().unwrap() {
            let len = record.len() as i64;
            len_vector.push(len);
            reads += 1;
            bases += record.len() as i64;
            num_n += get_n_bases(record.seq().as_bytes() );
            qual20 += get_qual_bases(record.qual().as_bytes(), 53); // 33 offset
            //qual30 += get_qual_bases(record.qual().as_bytes(), 63);
        
        }

        let q20 = qual20 as f64 / bases as f64 * 100.0;
        //let q30 = qual30 as f64 / bases as f64 * 100.0;
        let n50 = get_nx(&mut len_vector, 0.5);
        let min_len = len_vector.iter().min().unwrap();
        let max_len = len_vector.iter().max().unwrap();

        println!("reads\tbases\tn_bases\tmin_len\tmax_len\tN50\tQ20_percent");
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}", reads, bases, num_n, min_len, max_len, n50, q20);
    }
}