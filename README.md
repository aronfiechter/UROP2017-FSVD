# The Farthest Line-Segment Voronoi Diagram - UROP 2017

Summer internship (UROP: Undergraduate Research Opportunity Program) at
[USI Lugano](http://inf.usi.ch): Implementation of Geometrical Algorithms -
The Farthest Line-Segment Voronoi Diagram (supervised by Prof. Evanthia
Papadopoulou, Ioannis Mantas, Martin Suderland).


## Description

## Setup

### Prerequisites
### Installation



## TODO list

- [x] Read David Mount notes on Computational geometry:
  - [x] Introduction
  - [x] Convex hulls
  - [x] Line segment intersections
  - [x] Polygon triangulation (partially)
  - [x] Halfplane intersection and point-line duality
  - [ ] Linear Programming
  - [x] Voronoi diagrams and Fortune's algorithm
- [x] Read chapter 7 of Computational Geometry book
  - [x] Definitions and properties
  - [x] Computing the Voronoi diagram
  - [ ] Voronoi diagram of line segments (TODO: read again)
  - [x] Farthest-point Voronoi diagram
- [x] Read and understand Aurenhammer 2006 paper
- [ ] Understand how Ipelets work
  - [x] walk through the build process
  - [ ] read documentation
  - [x] compile copy of an ipelet myself
  - [ ] implement simple test ipelet
- [ ] Start implementing algorithm from Aurenhammer paper


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
files.

### CGAL

To see prerequisites: `brew info cgal`
To install: `brew install cgal`
Then download the source code from the official website (get the same version
as in `brew`), extract it where you need it, and:
```shell
cd CGAL-X.Y # version
cmake .
make
```

### Ipe and Ipelets


### C++


--------------------------------------------------------------------------------
## Temporal log

### 2017-07-03 to 2017-07-07

Read lecture notes, book chapter and paper

### 2017-07-10

Copied entire folder `OurGroupIpelets` to `old/OurGroupIpelets_SAVE`. Then
modified original folder; removed old files from `FCVD_HVD_Bisectors`.
Files removed:
  * `bisectors 2.cpp`
  * `bisectors very old.cpp`
  * `bisectors1.cpp`
  * `points-2.ipe`
  * `temp`
  * `inputs/*`

Created new folder TestIpelet to try creting an Ipelet from scratch. The only
needed files are:
  * test.cpp
  * lua/libCGAL_testIpelet.lua
  * CMakeLists.txt (copy of the file in `FCVD_HVD_Bisectors`)

### 2017-07-11
Copied example `simple_triangulation.cpp` ipelet from documentation and tried to
compile it and add it to Ipe. It works.
Took `min_cricle.cpp` example from CGAL documentation, converted it into Ipelet.
Both Ipelets are in `Test_Ipelets/`.
