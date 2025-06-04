use strobemers_rs::{MinStrobes, Result};

fn main() -> Result<()> {
    // ---------------------------------------------------- //
    // 1. Input sequence and parameters
    // ---------------------------------------------------- //
    // Sequence to process (25 nucleotides)
    let seq = "ATCGTACGATGCATGCATGCTGACG";
    // Order of strobemer: 2 (dimer) or 3 (trimer)
    let n        = 2;
    // Length of each k-mer (strobe). Must be between 1 and 64.
    let l        = 6;
    // Minimum window offset for selecting the next strobe
    let w_min    = 4;
    // Maximum window offset for selecting the next strobe
    let w_max    = 12;

    // ---------------------------------------------------- //
    // 2. Instantiate a MinStrobes iterator
    // ---------------------------------------------------- //
    // Converts the sequence to a byte slice, validates parameters,
    // and precomputes k-mer hashes and sliding-window minima.
    let mut ms = MinStrobes::new(seq.as_bytes(), n, l, w_min, w_max)?;

    println!("# MinStrobes example");
    println!("sequence      : {}", seq);
    println!("order (n)     : {}", n);
    println!("strobe length : {}", l);
    println!("w_min, w_max  : {}, {}", w_min, w_max);
    println!();

    // ---------------------------------------------------- //
    // 3. Iterate over MinStrobes, printing indices and hash
    // ---------------------------------------------------- //
    // Header for the output table
    println!(" idx | m1  m2 (m3) | hash");
    println!("-----+-------------+-------------------------");

    // Loop until the iterator returns None
    while let Some(hash) = ms.next() {
        // Retrieve the indices of the most recent strobes:
        // m1 is the starting k-mer index, m2 (and m3 if n=3) are chosen next
        let [m1, m2, m3] = ms.indexes();

        // Format the output differently depending on the strobemer order
        match n {
            2 => {
                // For order 2, only m1 and m2 are used
                // Print m1 twice for clarity (start and first strobe),
                // then m2, followed by the hash in hexadecimal
                println!("{:4} | {:3} {:3}     | 0x{:016x}", m1, m1, m2, hash)
            }
            3 => {
                // For order 3, print m1, m2, and m3
                println!("{:4} | {:3} {:3} {:3} | 0x{:016x}", m1, m1, m2, m3, hash)
            }
            _ => unreachable!(), // Parameter validation ensures n is 2 or 3
        }
    }

    Ok(())
}
