# BORIS - Borel Resummation with conformal mapping

## Compilation and running

Clone the repository and compile with

```
git clone https://github.com/apik/boris.git
cd boris
cargo build --release
```

To check installation we provided number of internal tests available from

```
cargo test
```

To install prodeced binaries to working directory use

```
cargo install --path . --root INSTALL_DIR
```

Two binaries will be available `INSTALL_DIR/bin/point` and `INSTALL_DIR/bin/scan`

## Description of implemented operation modes

### Single parameter mode

`cargo run --release --bin point -- input/point.toml`

or with

`INSTALL_DIR/bin/point input/point.toml`

Example config file should contain information about original series coefficients, target resummation point(`ep`), resummation algorithm parameters(`a,b,lam,q`) 

```toml
# info string
info = "Single point mode"
# # series
f = [ 2.0, -0.25, -0.0859375, 0.11442530392291644, -0.2875129509289455, 0.9561331447210607, -3.855754549203757 ]
# ep
ep=1
# a
a=0.375
# b
b   = 10.0
# lambda
lam = 1.42
# q
q=0.04
```
`info` is arbitrry text string passed directly to the output

### Scan parameter mode

`cargo run --release --bin scan -- -w [WORKERS] input/scan.toml` 

or with 

`INSTALL_DIR/bin/scan -w [WORKERS] input/scan.toml`

In the scan mode parallel execution implemented with OpenMP is supported. Parameter `[WORKERS]` sets number of parallel threads and default value is equal to the number of CPU cores.

Sample config file is the following:
```toml
###
# Info string to be printed during program run
info = "nu^-1(n=0) = 1.70251"

###
# Series to be resumed, starts with constant term
f   = [2.0, -0.25, -0.0859375, 0.11442530392291644, -0.2875129509289455, 0.9561331447210607, -3.855754549203757]

### Fixed resummation parameters
# ep - final point
ep = 1
#  a - general fixed by symmetry considerations
a   = 0.375

### Resummation parameters to be scanned
# b [min, max, subdivisions]
b   = [0.0, 40.0, 80]
# lambda
lam = [0.0, 4.5, 225]
# q = q/2
q   = [0.0, 0.4, 40]


### Scan plateaus parameters in number of subdivisions units 
deltab   = 40
deltalam = 50
deltaq   = 1


### If provided results are saved in <dump> as .tsv file: b,lam,q,f(l),f(l-1),f(l-2)
# dump = gridNU

### Save final results
# output file name
out = "nuInvN0.m"


### Integration order
# set order of untegration rules, default 20
#n = 40
```

Now resummation parameters `b,lam,q` are defined on a grid with syntax `[min,max,divisions]`.

After resummation in all grid points total error estimated and point minimaly sensitive to the resummation parameter detected. 
Parameters controling are `deltab,deltalam,deltaq` and are in units of grid steps.

Resummation results on all grid points can be dumped with `dump= FILENAME` and final result is saved in `out=OUTFILE` in the form `{res_min, res_mean, sigma_mean}`




