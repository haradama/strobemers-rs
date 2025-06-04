//! These tests port the Go implementation’s MinStrobes tests into Rust,
//! ensuring that the Rust version produces at least as many strobes as the Go reference.
//! The tests validate both order-2 and order-3 MinStrobes over a fixed sequence.

use strobemers_rs::{MinStrobes, Result};

/// Fixed test sequence (ASCII bytes).
const SEQ: &[u8] = b"ACGATCTGGTACCTAG";
/// Length of each k-mer (strobe length).
const L: usize = 3;
/// Minimum window offset for selecting the next strobe.
const W_MIN: usize = 3;
/// Maximum window offset for selecting the next strobe.
const W_MAX: usize = 5;

/// Corresponds to the Go test `TestMinStrobesOrders2`.
///
/// Verifies that MinStrobes with order = 2 produces at least one strobemer
/// over the test sequence. In debug builds, it also prints each strobe’s
/// indices and hash values for manual inspection.
///
/// # Returns
/// - `Ok(())` if at least one strobe was generated
/// - `Err(...)` if `MinStrobes::new` fails or no strobes are produced
#[test]
fn minstrobes_order2() -> Result<()> {
    // Order-2 strobemers (dimers)
    let n = 2;
    let mut ms = MinStrobes::new(SEQ, n, L, W_MIN, W_MAX)?;

    let mut iter_count = 0usize;
    while let Some(h) = ms.next() {
        iter_count += 1;

        // Only print details when running in debug mode
        if cfg!(debug_assertions) {
            let [i1, i2, _] = ms.indexes(); // For order-2, the third index is unused
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

    // The reference Go test considers it successful if at least one strobe is produced.
    assert!(iter_count > 0, "iterator produced no items");
    Ok(())
}

/// Corresponds to the Go test `TestMinStrobesOrder3`.
///
/// Verifies that MinStrobes with order = 3 produces at least one strobemer
/// over the test sequence. In debug builds, it also prints each strobe’s
/// indices and hash values for manual inspection.
///
/// # Returns
/// - `Ok(())` if at least one strobe was generated
/// - `Err(...)` if `MinStrobes::new` fails or no strobes are produced
#[test]
fn minstrobes_order3() -> Result<()> {
    // Order-3 strobemers (trimers)
    let n = 3;
    let mut ms = MinStrobes::new(SEQ, n, L, W_MIN, W_MAX)?;

    let mut iter_count = 0usize;
    while let Some(h) = ms.next() {
        iter_count += 1;

        // Only print details when running in debug mode
        if cfg!(debug_assertions) {
            let [i1, i2, i3] = ms.indexes();
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

    // The reference Go test considers it successful if at least one strobe is produced.
    assert!(iter_count > 0, "iterator produced no items");
    Ok(())
}
