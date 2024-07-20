[![Donate via Stripe](https://img.shields.io/badge/Donate-Stripe-green.svg)](https://buy.stripe.com/00gbJZ0OdcNs9zi288)<br>

## Metric Visualization

### [Launch](https://thenumbernine.github.io/metric/)

Currently uses spherical metric and displays basis along theta and phi.

Be sure to check out the [Lua version](https://github.com/thenumbernine/lua-metric).

soon to display additional information such as
- better lighting / models?
- better display of coordinate system ... embedding coordinates? popup text?
- some better way to show connection coefficients ... how when they're zero for the associated basis, the geodesic and partial line up
- maybe something to do with riemann metric tensor and geodesic deviation
- rendering of the mesh space and its dual space? multiply all points by $g\_{ij}$?
- more fomulas
- custom formulas - either by rewriting my CAS in javascript (or using a pre-existing one)
- - ... or send the formula to the server for lua to run it against my own lua CAS here: https://github.com/thenumbernine/symmath-lua
- - ... or emulate that lua CAS in javascript using Emscripten
