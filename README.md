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
  - [ ] implement simple test ipelet



--------------------------------------------------------------------------------
## Temporal log

### 2017-07-03
### 2017-07-04
### 2017-07-05
### 2017-07-06
### 2017-07-07

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
