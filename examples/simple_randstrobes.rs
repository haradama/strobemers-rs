use strobemers_rs::{RandStrobes, Result};

fn main() -> Result<()> {
    // ---------------------------------------------------- //
    // 1. Input sequence and strobemer parameters
    // ---------------------------------------------------- //
    // Sequence to process (25 bases)
    let seq = "ATCGTACGATGCATGCATGCTGACG";
    // Order of strobemer: 2 (dimer) or 3 (trimer)
    let n = 2;
    // Instantiate a RandStrobes iterator with:
    // - sequence as bytes
    // - order = 2
    // - k-mer length = 6
    // - minimum window offset = 4
    // - maximum window offset = 12
    let mut rs = RandStrobes::new(seq.as_bytes(), n, 6, 4, 12)?;

    // ---------------------------------------------------- //
    // 2. Iterate over RandStrobes, printing indices and hash
    // ---------------------------------------------------- //

    // Header for the output table
    println!(" idx | m1  m2 (m3) | hash");
    println!("-----+-------------+-------------------------");

    // Loop until the iterator returns None
    while let Some(hash) = rs.next() {
        // Retrieve the indices of the most recent strobes:
        // m1 is the starting k-mer index,
        // m2 (and m3 if n=3) are chosen next.
        let [m1, m2, m3] = rs.indexes();

        // Print differently depending on the strobemer order
        match n {
            2 => {
                // For order 2: print m1 twice (as start and first strobe), then m2,
                // followed by the combined hash in hexadecimal.
                println!("{:4} | {:3} {:3}     | 0x{:016x}", m1, m1, m2, hash);
            }
            3 => {
                // For order 3: print m1, m2, and m3, then the hash.
                println!("{:4} | {:3} {:3} {:3} | 0x{:016x}", m1, m1, m2, m3, hash);
            }
            _ => unreachable!(), // Parameter validation ensures n is 2 or 3
        }
    }

    Ok(())
}
