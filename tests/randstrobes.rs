//! These tests port the Go implementation’s RandStrobes tests into Rust,
//! ensuring that the Rust version produces at least one strobemer per order.
//! The tests validate both order-2 and order-3 RandStrobes over a fixed sequence.

use strobemers_rs::{RandStrobes, Result};

/// Fixed test sequence (ASCII bytes).
const SEQ: &[u8] = b"ACGATCTGGTACCTAG";
/// Length of each k-mer (strobe length).
const L: usize = 3;
/// Minimum window offset for selecting the next strobe.
const W_MIN: usize = 3;
/// Maximum window offset for selecting the next strobe.
const W_MAX: usize = 5;

/// Corresponds to the Go test `TestRandStrobesOrder2`.
///
/// Verifies that RandStrobes with order = 2 produces at least one strobemer
/// over the test sequence. In debug builds, it prints each strobe’s indices
/// and hash values for manual inspection.
///
/// # Returns
/// - `Ok(())` if at least one strobemer was generated.
/// - `Err(...)` if `RandStrobes::new` fails or no strobes are produced.
#[test]
fn randstrobes_order2() -> Result<()> {
    // Order-2 strobemers (dimers)
    let n = 2;
    // Instantiate a RandStrobes iterator with the fixed sequence and parameters.
    let mut rs = RandStrobes::new(SEQ, n, L, W_MIN, W_MAX)?;

    let mut iter_count = 0usize;
    while let Some(h) = rs.next() {
        iter_count += 1;

        // Only print details when running in debug mode (`cargo test -- --nocapture`)
        if cfg!(debug_assertions) {
            let [i1, i2, _] = rs.indexes(); // For order-2, the third index is unused
            // Print the full sequence
            println!("{}", std::str::from_utf8(SEQ).unwrap());
            // Print the first k-mer (m1) with its starting index
            println!(
                "{:>width$}{} i1:{}",
                "",
                std::str::from_utf8(&SEQ[i1..i1 + L]).unwrap(),
                i1,
                width = i1
            );
            // Print the second k-mer (m2) with its starting index
            println!(
                "{:>width$}{} i2:{}",
                "",
                std::str::from_utf8(&SEQ[i2..i2 + L]).unwrap(),
                i2,
                width = i2
            );
            // Print the combined hash value aligned to the right
            println!("{:>width$}{}\n", "", h, width = SEQ.len() + 1);
        }
    }

    // Ensure that at least one strobemer was generated, matching the Go reference.
    assert!(iter_count > 0, "iterator produced no items");
    Ok(())
}

/// Corresponds to the Go test `TestRandStrobesOrder3`.
///
/// Verifies that RandStrobes with order = 3 produces at least one strobemer
/// over the test sequence. In debug builds, it prints each strobe’s indices
/// and hash values for manual inspection.
///
/// # Returns
/// - `Ok(())` if at least one strobemer was generated.
/// - `Err(...)` if `RandStrobes::new` fails or no strobes are produced.
#[test]
fn randstrobes_order3() -> Result<()> {
    // Order-3 strobemers (trimers)
    let n = 3;
    // Instantiate a RandStrobes iterator with the fixed sequence and parameters.
    let mut rs = RandStrobes::new(SEQ, n, L, W_MIN, W_MAX)?;

    let mut iter_count = 0usize;
    while let Some(h) = rs.next() {
        iter_count += 1;

        // Only print details when running in debug mode (`cargo test -- --nocapture`)
        if cfg!(debug_assertions) {
            let [i1, i2, i3] = rs.indexes();
            // Print the full sequence
            println!("{}", std::str::from_utf8(SEQ).unwrap());
            // Print the first k-mer (m1) with its starting index
            println!(
                "{:>width$}{} i1:{}",
                "",
                std::str::from_utf8(&SEQ[i1..i1 + L]).unwrap(),
                i1,
                width = i1
            );
            // Print the second k-mer (m2) with its starting index
            println!(
                "{:>width$}{} i2:{}",
                "",
                std::str::from_utf8(&SEQ[i2..i2 + L]).unwrap(),
                i2,
                width = i2
            );
            // Print the third k-mer (m3) with its starting index
            println!(
                "{:>width$}{} i3:{}",
                "",
                std::str::from_utf8(&SEQ[i3..i3 + L]).unwrap(),
                i3,
                width = i3
            );
            // Print the combined hash value aligned to the right
            println!("{:>width$}{}\n", "", h, width = SEQ.len() + 1);
        }
    }

    // Ensure that at least one strobemer was generated, matching the Go reference.
    assert!(iter_count > 0, "iterator produced no items");
    Ok(())
}
