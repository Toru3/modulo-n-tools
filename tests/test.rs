use modulo_n_tools::montgomery::*;
use modulo_n_tools::*;
use num::BigInt;
use rand::Rng;
use rug::Integer;
use std::convert::TryInto;
use test_case::test_case;

#[test_case(100; "small")]
#[test_case(10000; "medium")]
#[test_case(1000000; "big")]
fn add_mod_test(n: usize) {
    let mut rng = rand::thread_rng();
    for _ in 0..n {
        let m = rng.gen::<i64>().abs() / 2 + 1;
        let a = rng.gen::<i64>() % m;
        let b = rng.gen::<i64>() % m;
        assert_eq!(add_mod(&a, &b, &m), (a + b) % m);
    }
}

#[test_case(100; "small")]
#[test_case(10000; "medium")]
#[test_case(1000000; "big")]
fn sub_mod_test(n: usize) {
    let mut rng = rand::thread_rng();
    for _ in 0..n {
        let m = rng.gen::<i64>().abs() / 2 + 1;
        let a = rng.gen::<i64>() % m;
        let b = rng.gen::<i64>() % m;
        assert_eq!(sub_mod(&a, &b, &m), (a - b) % m);
    }
}

#[test_case(100; "small")]
#[test_case(10000; "medium")]
#[test_case(1000000; "big")]
fn mul_mod_test(n: usize) {
    let mut rng = rand::thread_rng();
    for _ in 0..n {
        let m = rng.gen::<i64>().abs() / 2 + 1;
        let a = i64::from(rng.gen::<i32>()) % m;
        let b = i64::from(rng.gen::<i32>()) % m;
        assert_eq!(mul_mod(&a, &b, &m), (a * b) % m);
    }
}

fn is_prime(n: i64) -> bool {
    if n == 2 {
        return true;
    }
    if n < 2 || n % 2 == 0 {
        return false;
    }
    let mut i = 3;
    while i * i <= n {
        if n % i == 0 {
            return false;
        }
        i += 2;
    }
    true
}

#[test]
fn pow_mod_test_small() {
    assert_eq!(pow_mod(2, 15, &30), 8);
    assert_eq!(pow_mod(3, 53, &31), 11);
    assert_eq!(pow_mod(6, 12, &77), 36);
    assert_eq!(pow_mod(2, 10, &16), 0);
    assert_eq!(pow_mod(6, 10, &12), 0);
}

#[test_case(100; "small")]
#[test_case(10000; "medium")]
#[test_case(100000; "big")]
fn pow_mod_test(n: usize) {
    let mut rng = rand::thread_rng();
    for _ in 0..n {
        let m = loop {
            let m = i64::from(rng.gen::<i32>().abs() / 2 + 1);
            if is_prime(m) {
                break m;
            }
        };
        let a = rng.gen::<i64>() % m;
        assert_eq!(pow_mod(a, m, &m), a);
    }
}

#[test]
fn montgomery32_test_small() {
    let m = Montgomery32::new(31);
    assert_eq!(m.powmod(3, 53), 11);
    assert_eq!(m.powmod(5, 9), 1);
    let m = Montgomery32::new(77);
    assert_eq!(m.powmod(6, 12), 36);
    assert_eq!(m.powmod(17, 9), 13);
}

#[test_case(100; "small")]
#[test_case(10000; "medium")]
#[test_case(100000; "big")]
fn montgomery32_test(n: usize) {
    let mut rng = rand::thread_rng();
    for _ in 0..n {
        let m = loop {
            let m = rng.gen::<u32>() / 2;
            if m != 2 && is_prime(m.into()) {
                break m;
            }
        };
        let mon = Montgomery32::new(m);
        let a = rng.gen::<u32>() % m;
        assert_eq!(mon.powmod(a, m), a);
    }
}

#[test]
fn montgomery64_test_small() {
    let m = Montgomery64::new(31);
    assert_eq!(m.powmod(3, 53), 11);
    assert_eq!(m.powmod(5, 9), 1);
    let m = Montgomery64::new(77);
    assert_eq!(m.powmod(6, 12), 36);
    assert_eq!(m.powmod(17, 9), 13);
}

#[test_case(100; "small")]
#[test_case(10000; "medium")]
#[test_case(100000; "big")]
fn montgomery64_test(n: usize) {
    let mut rng = rand::thread_rng();
    for _ in 0..n {
        let m = loop {
            let m = u64::from(rng.gen::<u32>());
            if m != 2 && is_prime(m.try_into().unwrap()) {
                break m;
            }
        };
        let mon = Montgomery64::new(m);
        let a = rng.gen::<u64>() % m;
        assert_eq!(mon.powmod(a, m), a);
    }
}

fn montgomery_num() {
    let m: BigInt = (BigInt::from(1) << 107) - 1;
    let mon = Montgomery::new(m.clone());
    let a = BigInt::from(2);
    let p = m - 1;
    let r = mon.powmod(a, p);
    assert_eq!(r, BigInt::from(1));
}

fn montgomery_rug() {
    let m = Integer::from((Integer::from(1) << 107) - 1);
    let mon = Montgomery::new(m.clone());
    let a = Integer::from(2);
    let p = m - 1;
    let r = mon.powmod(a, p);
    assert_eq!(r, Integer::from(1));
}
