This is an experimental Minknowski Sum/Difference computation implemenetation based on JTS functionality.

Current implementation supports Minkowski sums between polygon/linestring/multipolygon/multilinestring/GeometryCollection and polygon/linestring, and minkowski difference between polygon/multipolygon/GeometryCollection and polygon/linestring 

It doesn't require polygons to be convex. The "source" polygon may contain holes. 

Any holes in "reference" polygon are ignored (in most cases it doesn't make much practical sense anyway).

Dut to change of package names in JTS 1.15, two versions of the code are provided for JTS 1.14 and 1.15
