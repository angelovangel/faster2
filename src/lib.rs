
pub fn get_n_bases(seq: &[u8]) -> i64 {
    let mut n = 0;
    for s in seq {
        if matches!(s, &78u8 | &110u8) { // N or n
            n += 1;
        }
    }
    n
}
// return per read content
pub fn get_gc_content(seq: &[u8]) -> f32 {
    
let mut n: i64 = 0;
    for s in seq {
        if matches!(s, &103u8 | &71u8 |  &99u8 | &67u8) { //G, g, C or c
            n += 1;
        }
    }
        n as f32 / seq.len() as f32
    }


pub fn get_qual_bases(q: &[u8], qx: u8) -> i64 {
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
pub fn qscore_probs(q: &[u8]) -> f32 {
    let mut qprob_sum = 0.0;
    for &item in q.iter() {
        let phred = *&item as f32 - 33.0;
        let prob = 10.0_f32.powf(-phred / 10.0);
        qprob_sum += prob
    }
    qprob_sum
}

// get geometric mean from phred scores
// not used, this is wrong
pub fn phred_gm(q: &[u8]) -> f64 {
    let mut phred_product = 0.;
    
    for &item in q {                       
    let phred = *&item as f64 - 33.0;

            phred_product += phred.ln();

    }
    let result = phred_product / q.len() as f64;
    result.exp()
}


pub fn get_nx(numbers: &mut [i64], fraction: f32) -> i64 {
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
