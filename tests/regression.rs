//! Regression tests that verify the Rust implementation produces the same
//! strobemer hash sequences as a previous snapshot (from the Go reference).
//! If the strobemer algorithm is intentionally changed, update these expected
//! values and document a breaking change in the CHANGELOG.

use strobemers_rs::{MinStrobes, RandStrobes, Result};

// ==== Fixed parameters (identical to the original Go tests) ====
const SEQ: &[u8] = b"ACGATCTGGTACCTAG"; // Test sequence
const L: usize = 3; // k-mer (strobe) length
const W_MIN: usize = 3; // Minimum window offset
const W_MAX: usize = 5; // Maximum window offset

// ==== Snapshot of expected hash outputs (as of 2025-06-03 commit) ====
//
// * MinStrobes order-2 … 11 values
// * MinStrobes order-3 … 6 values
// * RandStrobes order-2 … 11 values
// * RandStrobes order-3 … 6 values
//
const MIN_O2: [u64; 11] = [
    5508583604130516576,
    7820137869046132365,
    5541303490076687811,
    5796921065369559009,
    7864972478291945971,
    6364449594620396814,
    4156992363689746675,
    5730802552933835827,
    8690393705976365196,
    11912708257446301134,
    8953117104403771765,
];

const MIN_O3: [u64; 6] = [
    5838247918869859075,
    5824753939158295439,
    4305531019845332403,
    4497244201314985802,
    7896767773547654737,
    6896419184433288632,
];

const RAND_O2: [u64; 11] = [
    6508932193244882681,
    8820486458160498470,
    5796921065369559009,
    8819188626893971357,
    7864972478291945971,
    8337510363315416394,
    6747875957559703611,
    8321686146803792763,
    8690393705976365196,
    11912708257446301134,
    8953117104403771765,
];

const RAND_O3: [u64; 6] = [
    7772345821922645402,
    9313381998533055928,
    4497244201314985802,
    6763944872458295062,
    7896767773547654737,
    8376214760954553316,
];

// ---------------------------------------------------------------------
//                         REGRESSION  TESTS
// ---------------------------------------------------------------------

/// Verifies that MinStrobes of order 2 produces exactly the expected sequence of 11 hash values.
/// Collects all outputs from the iterator and compares them to the hard-coded snapshot.
#[test]
fn regression_minstrobes_order2() -> Result<()> {
    // Collect all order-2 MinStrobe hashes into a Vec
    let v: Vec<u64> = MinStrobes::new(SEQ, 2, L, W_MIN, W_MAX)?.collect();
    // Compare to the expected snapshot
    assert_eq!(v, MIN_O2);
    Ok(())
}

/// Verifies that MinStrobes of order 3 produces exactly the expected sequence of 6 hash values.
/// Collects all outputs from the iterator and compares them to the hard-coded snapshot.
#[test]
fn regression_minstrobes_order3() -> Result<()> {
    let v: Vec<u64> = MinStrobes::new(SEQ, 3, L, W_MIN, W_MAX)?.collect();
    assert_eq!(v, MIN_O3);
    Ok(())
}

/// Verifies that RandStrobes of order 2 produces exactly the expected sequence of 11 hash values.
/// Collects all outputs from the iterator and compares them to the hard-coded snapshot.
#[test]
fn regression_randstrobes_order2() -> Result<()> {
    let v: Vec<u64> = RandStrobes::new(SEQ, 2, L, W_MIN, W_MAX)?.collect();
    assert_eq!(v, RAND_O2);
    Ok(())
}

/// Verifies that RandStrobes of order 3 produces exactly the expected sequence of 6 hash values.
/// Collects all outputs from the iterator and compares them to the hard-coded snapshot.
#[test]
fn regression_randstrobes_order3() -> Result<()> {
    let v: Vec<u64> = RandStrobes::new(SEQ, 3, L, W_MIN, W_MAX)?.collect();
    assert_eq!(v, RAND_O3);
    Ok(())
}
