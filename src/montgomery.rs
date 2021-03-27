use core::convert::TryFrom;

/// Montgomery modular multiplication
///
/// ```
/// use modulo_n_tools::montgomery::{MontgomeryOperation, Montgomery64};
/// let m = Montgomery64::new(97);
/// let a = m.powmod(2, 77);
/// assert_eq!(a, 65);
/// ```
pub trait MontgomeryOperation<T, U> {
    /// $`NN' \equiv -1 \pmod R`$
    fn calc_n_prime(n: &T, s: u32) -> T;
    /// Constructor for $`{}\bmod N`$
    fn new(n: T) -> Self;
    /// Montgomery Reduction
    fn reduction(&self, x: U) -> T;
    /// $`x \mapsto xR \bmod N`$
    fn convert(&self, x: T) -> T;
    /// calcutate $`a^p \bmod N`$
    fn powmod<V>(&self, a: T, mut p: V) -> T
    where
        T: From<u8>,
        U: Clone + for<'x> core::ops::Mul<&'x U, Output = U> + From<T>,
        V: Clone
            + Ord
            + core::ops::ShrAssign<u32>
            + From<u8>
            + for<'x> core::ops::BitAnd<&'x V, Output = V>,
    {
        let c0 = V::from(0);
        let c1 = V::from(1);
        let mut x = self.convert(a);
        let mut y = self.convert(T::from(1));
        while p > c0 {
            let ux = U::from(x);
            if p.clone() & &c1 != c0 {
                let uy = U::from(y);
                y = self.reduction(ux.clone() * &uy);
            }
            x = self.reduction(ux.clone() * &ux);
            p >>= 1;
        }
        self.reduction(U::from(y))
    }
}

/// Montgomery $`R=2^{64}`$
///
/// ```
/// use modulo_n_tools::montgomery::{MontgomeryOperation, Montgomery64};
/// let m = Montgomery64::new(57);
/// let a = m.powmod(5, 42);
/// assert_eq!(a, 7);
/// ```
pub struct Montgomery64 {
    n: u64,
    np: u64,
    r2: u64,
}

impl MontgomeryOperation<u64, u128> for Montgomery64 {
    fn calc_n_prime(n: &u64, s: u32) -> u64 {
        let mut x = (!n + 1) % 8; // 3bit
        let mut b = 3;
        loop {
            let nx = n.wrapping_mul(x).wrapping_add(2).wrapping_mul(x);
            b *= 2;
            if b >= s {
                return nx;
            }
            x = nx;
        }
    }
    fn new(n: u64) -> Self {
        Montgomery64 {
            n,
            np: Self::calc_n_prime(&n, 64),
            r2: u64::try_from(u128::MAX % u128::from(n)).unwrap() + 1,
        }
    }
    fn reduction(&self, x: u128) -> u64 {
        let t = u64::try_from(x & u128::from(u64::MAX)).unwrap();
        let t = t.wrapping_mul(self.np);
        let t = x + u128::from(self.n) * u128::from(t);
        let t = u64::try_from(t >> 64).unwrap();
        if t >= self.n {
            t - self.n
        } else {
            t
        }
    }
    fn convert(&self, x: u64) -> u64 {
        self.reduction(u128::from(x) * u128::from(self.r2))
    }
}

/// Montgomery $`R=2^{32}`$
///
/// ```
/// use modulo_n_tools::montgomery::{MontgomeryOperation, Montgomery32};
/// let m = Montgomery32::new(89);
/// let a = m.powmod(3, 57);
/// assert_eq!(a, 23);
/// ```
pub struct Montgomery32 {
    n: u32,
    np: u32,
    r2: u32,
}

impl MontgomeryOperation<u32, u64> for Montgomery32 {
    fn calc_n_prime(n: &u32, s: u32) -> u32 {
        let mut x = (!n + 1) % 8; // 3bit
        let mut b = 3;
        loop {
            let nx = n.wrapping_mul(x).wrapping_add(2).wrapping_mul(x);
            b *= 2;
            if b >= s {
                return nx;
            }
            x = nx;
        }
    }
    fn new(n: u32) -> Self {
        Montgomery32 {
            n,
            np: Self::calc_n_prime(&n, 32),
            r2: u32::try_from(u64::MAX % u64::from(n)).unwrap() + 1,
        }
    }
    fn reduction(&self, x: u64) -> u32 {
        let t = u32::try_from(x & u64::from(u32::MAX)).unwrap();
        let t = t.wrapping_mul(self.np);
        let t = x + u64::from(self.n) * u64::from(t);
        let t = u32::try_from(t >> 32).unwrap();
        if t >= self.n {
            t - self.n
        } else {
            t
        }
    }
    fn convert(&self, x: u32) -> u32 {
        self.reduction(u64::from(x) * u64::from(self.r2))
    }
}

pub struct Montgomery<T> {
    n: T,
    np: T,
    s: u32,
    rm: T,
    r2: T,
}

fn bits<T>(mut n: T) -> u32
where
    T: Ord + core::ops::ShrAssign<u32> + core::ops::Shl<u32, Output = T> + From<u8>,
{
    let zero = T::from(0);
    let c64 = T::from(1) << 64;
    let mut b = 0;
    while n > c64 {
        n >>= 64;
        b += 64;
    }
    while n > zero {
        n >>= 1;
        b += 1;
    }
    b
}

impl<T> MontgomeryOperation<T, T> for Montgomery<T>
where
    T: Clone
        + Ord
        + for<'x> core::ops::BitAndAssign<&'x T>
        + for<'x> core::ops::AddAssign<&'x T>
        + for<'x> core::ops::SubAssign<&'x T>
        + for<'x> core::ops::MulAssign<&'x T>
        + for<'x> core::ops::BitAnd<&'x T, Output = T>
        + for<'x> core::ops::Sub<&'x T, Output = T>
        + for<'x> core::ops::Rem<&'x T, Output = T>
        + core::ops::Neg<Output = T>
        + core::ops::Shl<u32, Output = T>
        + core::ops::ShrAssign<u32>
        + core::ops::ShlAssign<u32>
        + From<u8>,
{
    fn calc_n_prime(n: &T, s: u32) -> T {
        let mut x = (-n.clone()) & &T::from(7);
        let mut b = 3;
        let rm = {
            let mut t = T::from(1);
            t <<= s;
            t -= &T::from(1);
            t
        };
        let two = T::from(2);
        loop {
            let mut nx = x.clone();
            nx *= n;
            nx += &two;
            nx *= &x;
            nx &= &rm;
            b *= 2;
            if b >= s {
                return nx;
            }
            x = nx;
        }
    }
    fn new(n: T) -> Self {
        let s = bits(n.clone());
        let rm = {
            let mut t = T::from(1);
            t <<= s;
            t -= &T::from(1);
            t
        };
        Montgomery::<T> {
            r2: (T::from(1) << (2 * s)) % &n,
            np: Self::calc_n_prime(&n, s),
            n,
            s,
            rm,
        }
    }
    fn reduction(&self, x: T) -> T {
        let mut t = x.clone();
        t &= &self.rm;
        t *= &self.np;
        t &= &self.rm;
        t *= &self.n;
        t += &x;
        t >>= self.s;
        if t >= self.n {
            t - &self.n
        } else {
            t
        }
    }
    fn convert(&self, mut x: T) -> T {
        x *= &self.r2;
        self.reduction(x)
    }
}
