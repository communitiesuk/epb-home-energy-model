# Using fuzzing

Fuzzing is a form of testing whereby a program is given large numbers of random inputs (either random bytes or in a specific format) to find whether any possible inputs might break the program. If a break is found, the input is shown so that a failing unit test can be created for that case. Once a failing case is found, it is also usually possible for a fuzzer to reduce the case down to the simplest form that will still cause the same break.

Property testing and fuzzing as concepts are similar, overlapping practices. Generally property testing will test whether code, given random inputs, has a particular abstract property, whereas fuzzing will just attempt to cause a break. However, it is possible to use a property testing approach with a fuzzer as well.

## Running the fuzzer on the HEM code

You will need to perform the [installation steps](https://rust-fuzz.github.io/book/cargo-fuzz/setup.html) on your machine to run fuzzing, which includes installing a nightly Rust compiler. (NB. you may find from time to time that using this nightly compiler causes weird crashes - in such circumstances you may wish to install a nightly from a previous and recent date using the command `rustup install nightly-YYYY-MM-DD`, with YYYY-MM-DD replaced by the particular date in question.)

To list the available fuzz targets (once cargo-fuzz has been installed):

```sh
cargo fuzz list
```

To run a target:

```sh
cargo fuzz run {TARGET_NAME}
```
