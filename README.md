# tqecd

The `tqecd` package is a spin-off from the [`tqec`](https://github.com/tqec/tqec) package that allows the automatic search of detectors in structured quantum error correction circuits.

This package was initially part of `tqec` and has been outsourced to make it accessible to anyone wanting to use it outside of the `tqec` library.

## Documentation

Documentation is available at [https://tqec.github.io/tqecd/index.html](https://tqec.github.io/tqecd/index.html).

## Installation

Currently, `tqecd` needs to be installed from source using

```sh
python -m pip install git+https://github.com/tqec/tqecd.git
```

The `tqecd` package has some dependencies that might be harder to install than a simple `pip install`. If you have any issues with the simple installation method above, please look at the [full installation page](https://tqec.github.io/tqecd/user_guide/installation.html).

## Basic usage

```py
import stim

from tqecd import annotate_detectors_automatically

# Example circuit that contains a detector
XXXX = stim.Circuit("""\
QUBIT_COORDS(1, 1) 0
QUBIT_COORDS(0, 0) 1
QUBIT_COORDS(2, 0) 2
QUBIT_COORDS(0, 2) 3
QUBIT_COORDS(2, 2) 4
RX 0 1 2 3 4
TICK
CX 0 1
TICK
CX 0 2
TICK
CX 0 3
TICK
CX 0 4
TICK
MX 0
""")
assert XXXX.num_detectors == 0

annotated_circuit = annotate_detectors_automatically(XXXX)
print(annotated_circuit)
```

should output the following quantum circuit:

```text
QUBIT_COORDS(1, 1) 0
QUBIT_COORDS(0, 0) 1
QUBIT_COORDS(2, 0) 2
QUBIT_COORDS(0, 2) 3
QUBIT_COORDS(2, 2) 4
RX 0 1 2 3 4
TICK
CX 0 1
TICK
CX 0 2
TICK
CX 0 3
TICK
CX 0 4
TICK
MX 0
DETECTOR(1, 1, 0) rec[-1]
```

## Contributing

Pull requests and issues are more than welcomed!

See the [contributing page](https://tqec.github.io/tqecd/contributor_guide.html) for specific contributing instructions.
