metric visualization

currently uses spherical metric and displays basis along theta and phi

soon to display additional information such as 
- better lighting / models?
- better display of coordinate system ... embedding coordinates? popup text? 
- some better way to show connection coefficients ... how when they're zero for the associated basis, the geodesic and partial line up
- maybe something to do with riemann metric tensor and geodesic deviation
- rendering of the mesh space and its dual space? multiply all points by g_ij?
- more fomulas
- custom formulas - either by rewriting my CAS in javascript (or using a pre-existing one) 
- - ... or send the formula to the server for lua to run it against my own lua CAS here: https://github.com/thenumbernine/symmath-lua
- - ... or emulate that lua CAS in javascript using Emscripten
