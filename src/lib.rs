#![no_std]
//! modulo_tools
//! ```
//! use modulo_n_tools::*;
//! use modulo_n_tools::montgomery::*;
//! let a = add_mod(&3, &4, &5);
//! assert_eq!(a, 2);
//! let b = mul_mod(&3, &a, &5);
//! assert_eq!(b, 1);
//! let c = pow_mod(2, 6, &7);
//! assert_eq!(c, 1);
//! let m = Montgomery64::new(57);
//! let d = m.powmod(5, 42);
//! assert_eq!(d, 7);
//! ```
use core::ops::{Add, AddAssign, BitAnd, Mul, Neg, Rem, ShrAssign, Sub, SubAssign};
pub mod montgomery;

fn reduce<T>(mut a: T, modulo: &T) -> T
where
    T: Ord + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Neg<Output = T>,
{
    if &a >= modulo {
        a -= modulo;
    } else if a <= -modulo {
        a += modulo
    }
    a
}

/// $`a + b \bmod n`$
///
/// Input: $`-\text{modulo} \leq a,\, b \leq \text{modulo}`$  
/// Output: $`-\text{modulo} \leq x \leq \text{modulo}`$
/// ```
/// use modulo_n_tools::add_mod;
/// assert_eq!(add_mod(&3, &4, &5), 2);
/// assert_eq!(add_mod(&2, &5, &6), 1);
/// assert_eq!(add_mod(&-3, &-2, &4), -1);
/// assert_eq!(add_mod(&2, &3, &5), 0);
/// ```
pub fn add_mod<T>(a: &T, b: &T, modulo: &T) -> T
where
    T: Ord + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T>,
{
    let c = a + b;
    reduce(c, modulo)
}

/// $`a - b \bmod n`$
///
/// Input: $`-\text{modulo} \leq a,\, b \leq \text{modulo}`$  
/// Output: $`-\text{modulo} \leq x \leq \text{modulo}`$
/// ```
/// use modulo_n_tools::sub_mod;
/// assert_eq!(sub_mod(&3, &4, &5), -1);
/// assert_eq!(sub_mod(&2, &-5, &6), 1);
/// assert_eq!(sub_mod(&-2, &-3, &4), 1);
/// assert_eq!(sub_mod(&2, &2, &5), 0);
/// ```
pub fn sub_mod<T>(a: &T, b: &T, modulo: &T) -> T
where
    T: Ord + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Sub<Output = T> + Neg<Output = T>,
{
    let c = a - b;
    reduce(c, modulo)
}

/// $`ab \bmod n`$
///
/// Input: $`-\text{modulo} \leq a,\, b \leq \text{modulo}`$  
/// Output: $`-\text{modulo} \leq x \leq \text{modulo}`$
/// ```
/// use modulo_n_tools::mul_mod;
/// assert_eq!(mul_mod(&3, &4, &5), 2);
/// assert_eq!(mul_mod(&2, &5, &6), 4);
/// assert_eq!(mul_mod(&-2, &-3, &4), 2);
/// assert_eq!(mul_mod(&2, &3, &6), 0);
/// ```
pub fn mul_mod<T>(a: &T, b: &T, modulo: &T) -> T
where
    for<'x> &'x T: Mul<Output = T> + Rem<Output = T>,
{
    &(a * b) % modulo
}

/// $`a^b \bmod n`$
///
/// Input: $`-\text{modulo} \leq a \leq \text{modulo}`$,
/// b is non-negative integer.  
/// Output: $`-\text{modulo} \leq x \leq \text{modulo}`$
/// ```
/// use modulo_n_tools::pow_mod;
/// assert_eq!(pow_mod(3, 4, &5), 1);
/// assert_eq!(pow_mod(2, 5, &6), 2);
/// assert_eq!(pow_mod(-2, 3, &4), 0);
/// assert_eq!(pow_mod(2, 3, &7), 1);
/// ```
pub fn pow_mod<T, U>(a: T, mut b: U, modulo: &T) -> T
where
    T: From<u8>,
    for<'x> &'x T: Mul<Output = T> + Rem<Output = T>,
    U: Ord + ShrAssign<u8> + From<u8>,
    for<'x> &'x U: BitAnd<Output = U>,
{
    let c0 = U::from(0);
    let c1 = U::from(1);
    let mut x = a;
    let mut y = T::from(1);
    while b > c0 {
        if &b & &c1 != c0 {
            y = mul_mod(&x, &y, modulo);
        }
        x = mul_mod(&x, &x, modulo);
        b >>= 1;
    }
    y
}

/// $`a\cdot b^p \bmod n`$
///
/// Input: $`-\text{modulo} \leq b \leq \text{modulo}`$,
/// c is non-negative integer.  
/// Output: $`-\text{modulo} \leq x \leq \text{modulo}`$
/// ```
/// use modulo_n_tools::mul_pow_mod;
/// assert_eq!(mul_pow_mod(1, 3, 4, &5), 1);
/// assert_eq!(mul_pow_mod(1, 2, 5, &6), 2);
/// assert_eq!(mul_pow_mod(1, -2, 3, &4), 0);
/// assert_eq!(mul_pow_mod(1, 2, 3, &7), 1);
/// ```
pub fn mul_pow_mod<T, U>(a: T, base: T, mut power: U, modulo: &T) -> T
where
    for<'x> &'x T: Mul<Output = T> + Rem<Output = T>,
    U: Ord + ShrAssign<u8> + From<u8>,
    for<'x> &'x U: BitAnd<Output = U>,
{
    let c0 = U::from(0);
    let c1 = U::from(1);
    let mut x = base;
    let mut y = a;
    while power > c0 {
        if &power & &c1 != c0 {
            y = mul_mod(&x, &y, modulo);
        }
        x = mul_mod(&x, &x, modulo);
        power >>= 1;
    }
    y
}
