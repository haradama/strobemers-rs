use strobemers_rs::{
    MinStrobes, RandStrobes, KmerHasher
};

struct XorHasher;

impl KmerHasher for XorHasher {
    fn hash_all(&self, seq: &[u8], k: usize) -> strobemers_rs::Result<Vec<u64>> {
        if k == 0 || seq.len() < k {
            return Err(strobemers_rs::StrobeError::SequenceTooShort);
        }

        let mut hashes = Vec::with_capacity(seq.len() - k + 1);
        for window in seq.windows(k) {
            let h = window.iter().fold(0u64, |acc, &b| acc ^ (b as u64));
            hashes.push(h);
        }
        Ok(hashes)
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let seq = b"ACGATCTGGTACCTAG";
    let order = 2;
    let k = 3;
    let w_min = 2;
    let w_max = 5;

    let xor_hasher = XorHasher;

    println!("== RandStrobes using XorHasher ==");
    let mut rs = RandStrobes::with_hasher(seq, order, k, w_min, w_max, &xor_hasher)?;
    for (i, h) in rs.by_ref().take(5).enumerate() {
        println!("RandStrobe {i}: {:x}", h);
    }

    println!("\n== MinStrobes using XorHasher ==");
    let mut ms = MinStrobes::with_hasher(seq, order, k, w_min, w_max, &xor_hasher)?;
    for (i, h) in ms.by_ref().take(5).enumerate() {
        println!("MinStrobe {i}: {:x}", h);
    }

    Ok(())
}
