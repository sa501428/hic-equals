# Comparing the content of .hic files

This is a quick script that can validate comparisons
of different .hic files (e.g. different versions or
minor floating point differences)

## Usage:
```
python3 main.py file1.hic file2.hic <NORMS>
```
NORMS: which normalizations to test (e.g. KR, SCALE, etc)
they should be comma separated, but with no spaces
(e.g. `SCALE,GW_SCALE`)

## Example:
```
python3 main.py gm12878_v1.hic gm12878_v2.hic SCALE,GW_SCALE
```

## Requirements:
- numpy
- scipy
- strawC>=0.0.9
