# The Farthest Line-Segment Voronoi Diagram - UROP 2017

Summer internship (UROP: Undergraduate Research Opportunity Program) at
[USI Lugano](http://inf.usi.ch): Implementation of Geometrical Algorithms -
The Farthest Line-Segment Voronoi Diagram (supervised by Prof. Evanthia
Papadopoulou, Ioannis Mantas, Martin Suderland).


## Description

## Setup

### Prerequisites
### Installation



--------------------------------------------------------------------------------
## Notes

### Compilation

To compile an own non-shipped program creating `executable` (adapt according to
CGAL version):
```shell
cd /path/to/program
cgal_create_CMakeLists -s executable
cmake -DCGAL_DIR=$HOME/CGAL-4.9.1 .
make
```

To compile an Ipelet copy the file `CMakeLists.txt` from
`.../CGAL-X.X/demo/CGAL_Ipelets/CMakeLists.txt`, then after
`set(CGAL_IPELETS ${CGAL_IPELETS})` delete the other Ipelets and add own. Don't
forget to create a `lua/` folder where to put the `libCGAL_ipeletname.lua`
files (as described on the
[Ipelets manual](http://ipe.otfried.org/manual/cpp-ipelets.html)).

### CGAL

To see prerequisites: `brew info cgal`
To install: `brew install cgal`
Then download the source code from the official website (get the same version
as in `brew`), extract it where you need it (for example your home directory),
and:
```shell
cd CGAL-X.Y # e.g. CGAL-4.9.1
cmake .
make
```

### Ipe and Ipelets


### C++
