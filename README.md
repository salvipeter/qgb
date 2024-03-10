# Quadratic GB patch

C0 GB patch with quadratic Bézier boundary curves and a settable midpoint.
Uses a regular domain; for four-sided surfaces it is the same as a quadratic
tensor-product Bézier patch.

The example program reads files of the following format:
```
<# of sides>
<1st side 1st CP>
<1st side 2nd CP>
<1st side 3rd CP>
<2st side 1st CP = 1st side 3rd CP>
<2st side 2nd CP>
...
<nth side 3rd CP = 1st side 1st CP>
[central CP]
```
(cf. `test.qgb`)

The default position of the central CP is the mass center of the 2nd control points.

As a library it should be used as follows:

1. Call constructor `QGB(size_t)` with the number of sides
1. Set all boundary curves with `void setBoundary(size_t)`
1. Set midpoint with either `void setMidpoint(const Point3D &)` or `void resetMidpoint()`
1. Evaluate at a given parameter with `Point3D eval(const Point2D &)`
   or at a given resolution with `TriMesh eval(size_t)`
