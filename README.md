# center: C Program For Centering Gromacs Simulations

Centers a selected group of atoms in a simulation box. Simpler, faster, and more robust than `gmx trjconv`.

Uses state-of-the-art algorithm for the calculation of center of geometry in periodic systems developed by Linge Bai & David Breen (https://doi.org/10.1080/2151237X.2008.10129266). Therefore, `center` can center any selection of atoms (that is not completely homogeneously distributed) into the center of the simulation box, no matter whether any of the selected atoms crosses box boundaries.

## Dependencies

`center` requires you to have groan library installed. You can get groan from [here](https://github.com/Ladme/groan). See also the [installation instructions](https://github.com/Ladme/groan#installing) for groan.

## Installation

1) Run `make groan=PATH_TO_GROAN` to create a binary file `center` that you can place wherever you want. `PATH_TO_GROAN` is a path to the directory containing groan library (containing `groan.h` and `libgroan.a`).
2) (Optional) Run `make install` to copy the the binary file `center` into `${HOME}/.local/bin`.

## Options

```
Usage: center -c GRO_FILE -o OUTPUT_FILE [OPTION]...

OPTIONS
-h               print this message and exit
-c STRING        gro file to read
-f STRING        xtc file to read (optional)
-n STRING        ndx file to read (optional, default: index.ndx)
-o STRING        output file name
-r STRING        selection of atoms centered (default: Protein)
-s INTEGER       only center every Nth frame (default: 1)
-x/-y/-z         center in individual x/y/z dimensions (default: center in xyz)
```

You can specify any selection of atoms for centering using the flag `-r` and the [groan selection language](https://github.com/Ladme/groan#groan-selection-language). 
Note that if the query contains multiple words, they must be enclosed into quotation marks or apostrophes:
```
center -c system.gro -r "name CA and resname LEU" -o system_centered.gro
```

## Example usage

```
center -c md.gro -f md.xtc -o md_centered.xtc -r Backbone -x -y -s 5
```

This command will read every 5th frame of the 'md.xtc' trajectory and in each such frame center a selection of atoms corresponding to ndx group 'Backbone' to the center of the simulation box in the _xy_ plane. The result will be written into 'md_centered.xtc'.

The position of atoms on the _z_ axis will not be changed. `center` will read ndx file named `index.ndx` (the default option) searching for the ndx group 'Backbone'. If the ndx file does not exist or if the ndx group 'Backbone' does not exist, `center` will complain that no atoms have been selected for centering and will exit.

Note that if an `xtc` file is supplied, atom coordinates from the `gro` file are not used at all.

## Limitations

Assumes that the simulation box is rectangular and that periodic boundary conditions are applied in all three dimensions.

Note that the `center` uses center of _geometry_, not center of _mass_, for centering.

Only tested on Linux. Probably will not work on anything that is not UNIX-like.