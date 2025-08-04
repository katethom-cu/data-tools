# Data Tools

A collection of python scripts I wrote for data analysis:


## Files
- `src/cavity_pert.py` - cavity perturbation dielectric measurements. designed
to operate on tdms files output from that one labview program on that computer
in C3.06
- `src/xrd.py` - analyses xrd based on `.acs` and `.xlsx` files exported from
higschore in a specific format.
- `src/hpr_ppr.py` - HPR and PPR dielectric measurement analysis. Extracts real
and imaginary permittivity from 
- `src/catpuccin.py` - colour scheme!
- `src/style.py` - specific matplotlib styling things.
- `templates/` templates showing examples of each type of analysis.
- `link.sh` install files as a symlink to `~/.local/lib/util/`

## Install
I like to simlink my files to `~/.local/lib/util/` and then i include that in
my pypath for the virtual environment common to all the data anlysis stuff i
do.
