# modulo_tools
add, sub, mul, pow in modulo, Montgomery multiplication.
```rust
use modulo_n_tools::*;
use modulo_n_tools::montgomery::*;
let a = add_mod(&3, &4, &5);
assert_eq!(a, 2);
let b = mul_mod(&3, &a, &5);
assert_eq!(b, 1);
let c = pow_mod(2, 6, &7);
assert_eq!(c, 1);
let m = Montgomery64::new(57);
let d = m.powmod(5, 42);
assert_eq!(d, 7);
```
