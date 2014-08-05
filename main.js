var canvas;
var gl;
var renderer;
var mouse;
var sceneObj;
var meshObj;

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	renderer.resize();

	var info = $('#info');
	var width = window.innerWidth 
		- parseInt(info.css('padding-left'))
		- parseInt(info.css('padding-right'));
	info.width(width);
	var height = window.innerHeight
		- parseInt(info.css('padding-top'))
		- parseInt(info.css('padding-bottom'));
	info.height(height - 32);
	$('#panel').height(height - 16);
}

function update() {
	//update
	renderer.draw();
	requestAnimFrame(update);
};

/*

mapping:
	y_i(x_j)
	coord = x_j
	returns y_i

diff:
	dy_i/dx_j
	 dim = j
	 coord = x_j
	 return dy_i/dx_j

unitDiff: diff, normalized algebraically to remove singularities


*/
var coordCharts = {
	/*
	metric:
		g_theta_theta = r^2
		g_phi_phi = r^2 sin^2 theta
	partials:
		g_phi_phi,theta = 2 r^2 sin(theta) cos(theta)
	connections:
		conn_phi_phi_theta = 1/2 g_phi_phi,theta = r^2 sin(theta) cos(theta)
		conn_phi_theta_phi = 1/2 g_phi_phi,theta = r^2 sin(theta) cos(theta)
		conn_theta_phi_phi = -1/2 g_phi_phi,theta = -r^2 sin(theta) cos(theta)
	2nd kind:
		conn^phi_phi_theta = conn^phi_theta_phi = cos(theta)/sin(theta)
		conn^theta_phi_phi = -sin(theta) cos(theta)
	*/
	Spherical : {
		constants : {r : 1},
		parameters : ['theta', 'phi'],
		coordinateMin : [0, -Math.PI],
		coordinateMax : [Math.PI, Math.PI],
		equations : [
			'r * sin(theta) * cos(phi)',
			'r * sin(theta) * sin(phi)',
			'r * cos(theta)'
		],
		initialCoord : [0.26623666555845704, 1.8215957403167709],
		initialDirection : [0, 1]
	},
	Polar : {
		initialCoord : [0.26623666555845704, 1.8215957403167709],
		initialDirection : [0, 1],
		coordinateMin : [0, -Math.PI],
		coordinateMax : [2, Math.PI],
		parameters : ['r', 'phi'],
		equations : [
			'r * cos(phi)',
			'r * sin(phi)'
		],
		initialCoord : [1, 0],
		initialDirection : [0, 1]
	},
	Torus : {
		constants : {
			r : .25,
			R : 1
		},
		parameters : ['theta', 'phi'],
		coordinateMin : [-Math.PI, -Math.PI],
		coordinateMax : [Math.PI, Math.PI],
		equations : [
			'(r * sin(theta) + R) * cos(phi)',
			'(r * sin(theta) + R) * sin(phi)',
			'r * cos(theta)'
		],
		initialCoord : [0.26623666555845704, 1.8215957403167709],
		initialDirection : [0, 1]
	},
	Paraboloid : {
		coordinateMin : [1,1],
		coordinateMax : [-1,-1],
		parameters : ['u', 'v'],
		equations : [
			'u', 'v', '-u * u - v * v'
		],
		initialCoord : [.25,.25],
		initialDirection : [0,1]
	},
	'Hyperbolic Paraboloid' : {
		parameters : ['u', 'v'],
		coordinateMin : [-1,-1],
		coordinateMax : [1,1],
		equations : [ 'u', 'v', 'u * u - v * v' ],
		initialCoord : [.25,.25],
		initialDirection : [0,1]
	}
};
var currentCoordChart;

var intDivs = [120,120];

var intCoordToCoord = function(intCoord) {
	var coord = [];
	var coordMin = currentCoordChart.coordinateMin;
	var coordMax = currentCoordChart.coordinateMax;
	for (var k = 0; k < intCoord.length; ++k) {
		coord[k] = intCoord[k] / (intDivs[k]-1) * (coordMax[k] - coordMin[k]) + coordMin[k];
	}
	return coord;
};

function reset() {
	//generate metric geometry
	var thetaDiv = intDivs[0];
	var phiDiv = intDivs[1];
	for (var itheta = 0; itheta < thetaDiv; ++itheta) {
		for (var iphi = 0; iphi < phiDiv; ++iphi) {
			var intCoord = [itheta, iphi];
			var coord = intCoordToCoord(intCoord);
			var mappedCoord = currentCoordChart.mapping(coord);
			var normal = currentCoordChart.unitDiff(coord, 2);
			for (var k = 0; k < 2; ++k) {
				meshObj.attrs.intCoord.data[k + 2 * (iphi + phiDiv * itheta)] = intCoord[k];
			}
			for (var k = 0; k < 2; ++k) {
				meshObj.attrs.coord.data[k + 2 * (iphi + phiDiv * itheta)] = coord[k];
			}
			for (var k = 0; k < 3; ++k) {
				meshObj.attrs.vertex.data[k + 3 * (iphi + phiDiv * itheta)] = mappedCoord[k];
			}
			for (var k = 0; k < 3; ++k) {
				meshObj.attrs.normal.data[k + 3 * (iphi + phiDiv * itheta)] = normal[k];
			}
		}
	}
	meshObj.attrs.intCoord.updateData();
	meshObj.attrs.coord.updateData();
	meshObj.attrs.vertex.updateData();
	meshObj.attrs.normal.updateData();
	
	currentCoord = currentCoordChart.initialCoord.clone();
	currentDirection = currentCoordChart.initialDirection.clone();
	//TODO mesh regen as well since we will want varying resolutions ... or not?
	selectCoord(currentCoord);
}

var selectionObj;
var selectionRes = 100;
var basisObjs = [];
//var connObjs = [];
var partialPathObj;
var geodesicPathObj;
var pathLength = 5000;
var currentCoord = [0.26623666555845704, 1.8215957403167709];
var currentDirection = [0, 1];
function selectCoord(coord) {
	currentCoord = coord.clone();
	var basis0 = currentCoordChart.diff(coord, 0);
	var basis1 = currentCoordChart.diff(coord, 1);
	var basis = [basis0, basis1];
	
	var unitBasis0 = currentCoordChart.unitDiff(coord, 0);
	var unitBasis1 = currentCoordChart.unitDiff(coord, 1);
	var normal = currentCoordChart.unitDiff(coord,2);
	var unitBasis = [unitBasis0, unitBasis1];
	
	var mappedCoord = currentCoordChart.mapping(coord);

	var selectionRadius = .1;
	var modv = vec3.create();
	for (var i = 0; i < selectionRes; ++i) {
		var psi = i/selectionRes * Math.PI * 2;
		vec3.scaleAndAdd(modv, mappedCoord, unitBasis0, selectionRadius * Math.cos(psi));
		vec3.scaleAndAdd(modv, modv, unitBasis1, selectionRadius * Math.sin(psi));
		for (var j = 0; j < 3; ++j) {
			selectionObj.attrs.vertex.data[3 * i + j] = modv[j];
		}
	}
	selectionObj.attrs.vertex.updateData();

	var basisLength = .4;
	for (var i = 0; i < 2; ++i) {
		var basisObj = basisObjs[i];
		for (var k = 0; k < 3; ++k) {
			basisObj.attrs.vertex.data[k] = mappedCoord[k];
			basisObj.attrs.vertex.data[3 + k] = mappedCoord[k] + basisLength * basis[i][k] + .01 * normal[k];
		}
		basisObj.attrs.vertex.updateData();
		/* need a better way to visualize this ...
		for (var j = 0; j < 2; ++j) {
			var conn = currentCoordChart.connection(coord, i, j);
			var connObj = connObjs[i][j];
			for (var k = 0; k < 3; ++k) {
				connObj.attrs.vertex.data[k] = mappedCoord[k] + basisLength * basis[i][k];
				connObj.attrs.vertex.data[3+k] = mappedCoord[k] + basisLength * (basis[i][k] + conn[0] * basis[0][k] + conn[1] * basis[1][k]);
			}
			connObj.attrs.vertex.updateData();
		}
		*/
	}

	chooseDirection(currentDirection);
}

function chooseDirection(direction) {
	currentDirection = direction.clone();
	var partialCoord = currentCoord.clone();
	var partialDir = direction.clone();
	var geoCoord = currentCoord.clone();
	var geoDir = direction.clone();
	var step = .001;
	for (var i = 0; i < pathLength; ++i) {
		partialCoord[0] += step * partialDir[0];
		partialCoord[1] += step * partialDir[1];
		var mapped = currentCoordChart.mapping(partialCoord);
		var normal = currentCoordChart.unitDiff(partialCoord,2);
		partialPathObj.attrs.vertex.data[3*i+0] = mapped[0] + .01 * normal[0];
		partialPathObj.attrs.vertex.data[3*i+1] = mapped[1] + .01 * normal[1];
		partialPathObj.attrs.vertex.data[3*i+2] = mapped[2] + .01 * normal[2];
	}
	partialPathObj.attrs.vertex.updateData();
	for (var i = 0; i < pathLength; ++i) {
		geoCoord[0] += step * geoDir[0];
		geoCoord[1] += step * geoDir[1];
		var geoAccel = [0,0];
		for (var j = 0; j < 2; ++j) {
			for (var k = 0; k < 2; ++k) {
				for (var l = 0; l < 2; ++l) {
					geoAccel[j] -= currentCoordChart.connection(geoCoord,k,l)[j] * geoDir[k] * geoDir[l];
				}
			}
		}
		geoDir[0] += step * geoAccel[0];
		geoDir[1] += step * geoAccel[1];
		var mapped = currentCoordChart.mapping(geoCoord);
		var normal = currentCoordChart.unitDiff(geoCoord,2);
		geodesicPathObj.attrs.vertex.data[3*i+0] = mapped[0] + .01 * normal[0];
		geodesicPathObj.attrs.vertex.data[3*i+1] = mapped[1] + .01 * normal[1];
		geodesicPathObj.attrs.vertex.data[3*i+2] = mapped[2] + .01 * normal[2];
	}
	geodesicPathObj.attrs.vertex.updateData();
}

var findClickedCoord;
(function(){
	var tmp = vec3.create();
	var tmp2 = vec3.create();
	findClickedCoord = function() {
		var mouseDir = vec3.create();
		renderer.mouseDir(mouseDir, mouse.xf, mouse.yf);
		vec3.normalize(mouseDir, mouseDir);
		//ray intersection test with coordinate chart ...
		// for now just search for closest point in geometry?
		var bestDist = undefined;
		var bestCoord = undefined;
		vec3.quatZAxis(tmp2, renderer.view.angle);
		var considered = 0;
		var pt = vec3.create();
		for (var i = 0; i < meshObj.attrs.vertex.data.length/3; ++i) {
			var normal = [
				meshObj.attrs.normal.data[3*i+0],
				meshObj.attrs.normal.data[3*i+1],
				meshObj.attrs.normal.data[3*i+2]
			];
			vec3.transformQuat(normal, normal, sceneObj.angle);
			if (vec3.dot(normal, tmp2) < 0) continue;	//forward dot normal > 0 means the surface is back-facing
			
			var vertex = [
				meshObj.attrs.vertex.data[3*i+0],
				meshObj.attrs.vertex.data[3*i+1],
				meshObj.attrs.vertex.data[3*i+2]
			];
			vec3.transformQuat(pt, vertex, sceneObj.angle);

			
			//make sure we're on the right side of the view plane
			vec3.sub(tmp, pt, renderer.view.pos);
			if (vec3.dot(tmp, tmp2) > 0) continue;	//fwd dot delta > 0 means we're good, so -fwd dot delta < 0 means we're good, so -fwd dot delta > 0 means we're bad
		
			considered++;

			//ray/point distance from view pos / mouse line
			vec3.sub(tmp, pt, renderer.view.pos);
			vec3.cross(tmp, tmp, mouseDir);
			var dist = vec3.length(tmp);
			
			if (bestDist === undefined || dist < bestDist) {
				bestDist = dist;
				bestCoord = [
					meshObj.attrs.coord.data[2*i+0],
					meshObj.attrs.coord.data[2*i+1]
				];
			}
		}
		return bestCoord;
	}
})();

$(document).ready(function() {
	$('#panelButton').click(function() {
		$('#panel').show();	
		$('#info').hide();
	});
	$('#infoButton').click(function() {
		$('#info').show();
		$('#panel').hide();
	});
	
	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(canvas).disableSelection()

	try {
		renderer = new GL.CanvasRenderer({canvas:canvas});
		gl = renderer.context;
	} catch (e) {
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}
	$('#menu').show();
	$('#panel').show();

	if ($.url().param('info')) {
		$('#info').show();
		$('#panel').hide();
	}

	window.inputState = 'rotate';

	$('#tools_rotate').click(function() { inputState = 'rotate'; });
	$('#tools_select').click(function() { inputState = 'select'; });
	$('#tools_direction').click(function() { inputState = 'direction'; });

	currentCoordChart = coordCharts.Spherical;

	//init lua
	console.log('loading lua');
	var luaDoneLoading = false;
	var lua = new EmbeddedLuaInterpreter({
		packages : ['ext', 'symmath'],
		packageTests : ['symmath'],
		done : function() {
			console.log('loaded lua');
			Lua.execute("package.path = package.path .. ';./?/init.lua'");
			this.executeAndPrint("require 'symmath'");
			this.executeAndPrint("for k,v in pairs(symmath) do _G[k] = v end");
			luaDoneLoading = true;
		},
		autoLaunch : true
	});
		
	var coordLabels = ['x', 'y', 'z'];

	var allInputs = coordLabels.map(function(x) { return 'equation'+x.toUpperCase(); }).concat(['parameters', 'constants']).map(function(id) { return $('#'+id); });
	
	//args: done: what to execute next
	var updateEquations = function() {
		var doUpdateEquations = function() {
			//args: callback = what to execute, output = where to redirect output, error = where to redirect errors
			var capture = function(args) {
				//now cycle through coordinates, evaluate data points, and get the data back into JS
				//push module output and redirect to a buffer of my own
				var oldPrint = lua.print;
				var oldError = lua.printErr;
				if (args.output !== undefined) lua.print = args.output;
				if (args.error !== undefined) lua.printErr = args.error;
				args.callback();
				lua.print = oldPrint;
				lua.printErr = oldError;
			};
			
			//declare parameter variables
			var parameters = $('#parameters').val().split(',').map(function(s) { return s.trim(); });
			$.each(parameters, function(j,param) {
				lua.executeAndPrint(param+" = Variable('"+param+"')");
			});
			
			$('#constants').val().split(',').map(function(s) { 
				s = s.trim();
				var eqs = s.split('=').map(function(side) { return side.trim(); });
				var param = eqs[0];
				var value = eqs[1];
				lua.executeAndPrint(param+" = Constant("+value+")");
			});		
			
			var propertiesHTML = '';
			
			var luaEqnToJSFunc = function(eqn, params, texLabel) {
				var failed = false;
				//here I'm using output for errors
				//since directing error doesn't work -- all errors result in stdout printing "ERROR attempt to call string" 
				capture({
					callback : function() {
						Lua.execute("eqn = simplify("+eqn+")");
						Lua.execute("if type(eqn) == 'number' then eqn = Constant(eqn) end");
					},
					output : function(s) {
						//don't throw -- Lua.execute will catch it.
						console.log('Lua error! '+s);
						failed = true;
					}
				});
				if (failed) throw 'Lua error!';	//'return true' means break in $.each	
			
				//here I'm using output for capturing the compiled lua code
				// I'm recording errors if the captured code fails to compile in JavaScript
				var resultFunction = undefined;
				capture({
					callback : function() {
						var luaCmd = "print(symmath.ToJavaScriptCode:compile(eqn, {"+params+"}))";
						console.log('executing lua '+luaCmd);
						//print commands are going to the old output ...
						Lua.execute(luaCmd);
					},
					output : function(jsCmd) {
						console.log('got JS output '+jsCmd);
						try {
							resultFunction = eval(jsCmd);
						} catch (e) {
							failed = true;
						}
					}
				});
				if (failed) throw 'Lua error!';

				//while we're here, let's store the LaTex generated from the equations ...
				capture({
					callback : function() {
						var luaCmd = "print(symmath.ToLaTeX(eqn))"
						console.log('executing lua '+luaCmd);
						Lua.execute(luaCmd);
					},
					output : function(TeX) {
						console.log('got TeX output '+TeX);
						if (texLabel !== undefined) TeX = texLabel + ' = ' + TeX;
						propertiesHTML += '\\( ' + TeX + ' \\)<br>\n';
					}
				});

				return resultFunction;
			};

			//compile equations here
			var coordCallbacks = [];	//coordCallbacks[xyz]
			var coordDerivs = [];		//coordDerivs[xyz][uv|normal]
			var metric = [];			//metric[uv][uv] = g_ij = e_i dot e_j
			var metricEqn = [];
			var invMetricEqn = [];		//invMetricEqn[uv][uv] = g^ij
			var diffMetricEqn = [];		//diffMetricEqn[uv][uv][uv] = g_ij,k
			var conn1stEqn = [];		//conn1stEqn[uv][uv][uv] = conn_ijk = 1/2 (g_ij,k + g_ik,j - g_jk,i)
			var conn2nd = [];			//conn2nd[uv][uv][uv] = conn^i_jk = g^il * conn_ljk
			var conn2ndEqn = [];
			console.log('generating coordinate chart functions...');
			var failed = false;
			var equations = [];
			for (var i = 0; i < coordLabels.length; ++i) {
				equations[i] = $('#equation'+coordLabels[i].toUpperCase()).val();
			}
			try {
				//get functions
				$.each(equations, function(i,eqn) {
					coordCallbacks[i] = luaEqnToJSFunc(eqn, parameters.join(', '), coordLabels[i]);
					coordDerivs[i] = [];
				});
				$.each(equations, function(i,eqn) {
					$.each(parameters, function(j,param) {
						coordDerivs[i][j] = luaEqnToJSFunc('diff('+eqn+', '+param+')', parameters.join(', '), '{{d'+coordLabels[i]+'} \\over {d'+param+'}}');
					});
				});

				//get metric and its derivative 
				$.each(parameters, function(i,u) {
					metric[i] = [];
					metricEqn[i] = [];
					diffMetricEqn[i] = [];
					$.each(parameters, function(j,v) {
						var sum = [];
						$.each(equations, function(k,eqn) {
							sum.push('diff('+eqn+', '+u+') * diff('+eqn+', '+v+')');
						});
						metricEqn[i][j] = sum.join(' + ');
						metric[i][j] = luaEqnToJSFunc(metricEqn[i][j], parameters.join(', '), 'g_{'+u+v+'}');
						diffMetricEqn[i][j] = [];
						$.each(parameters, function(k,w) {
							diffMetricEqn[i][j][k] = 'diff('+metricEqn[i][j]+', '+w+')';
						});
					});
				});

				//inverse metric
				if (parameters.length != 2) throw 'only works with two parameters';
				var metricDetEqn = '('+metricEqn[0][0]+') * ('+metricEqn[1][1]+') - ('+metricEqn[1][0]+') * ('+metricEqn[0][1]+')';
				var invMetricEqn = [
					[	'('+metricEqn[1][1]+') / ('+metricDetEqn+')',	//00
						'-('+metricEqn[0][1]+') / ('+metricDetEqn+')'],	//01
					[	'-('+metricEqn[1][0]+') / ('+metricDetEqn+')',	//10
						'('+metricEqn[0][0]+') / ('+metricDetEqn+')']
				];
	
				//1st conns
				$.each(parameters, function(i,u) {
					conn1stEqn[i] = [];
					$.each(parameters, function(j,v) {
						conn1stEqn[i][j] = [];
						$.each(parameters, function(k,w) {
							conn1stEqn[i][j][k] = '1/2 * ('+diffMetricEqn[i][j][k]+' + '+diffMetricEqn[i][k][j]+' - '+diffMetricEqn[j][k][i]+')';
						});
					});
				});

				//2nd conns
				$.each(parameters, function(i,u) {
					conn2nd[i] = [];
					conn2ndEqn[i] = [];
					$.each(parameters, function(j,v) {
						conn2nd[i][j] = [];
						conn2ndEqn[i][j] = [];
						$.each(parameters, function(k,w) {
							var sum = [];
							for (var l = 0; l < parameters.length; ++l) {
								sum.push('('+invMetricEqn[i][l]+') * ('+conn1stEqn[l][j][k]+')');
							}
							conn2ndEqn[i][j][k] = sum.join(' + ');
							conn2nd[i][j][k] = luaEqnToJSFunc(conn2ndEqn[i][j][k], parameters.join(', '), '{Gamma^'+u+'}_{'+v+w+'}');
						});
					});
				});
				
				//get normal
				$.each(equations, function(i,eqn) {
					var j = (i+1)%equations.length;
					var k = (j+1)%equations.length;
					var dj_du = coordDerivs[j][0];
					var dk_du = coordDerivs[k][0];
					var dj_dv = coordDerivs[j][1];
					var dk_dv = coordDerivs[k][1];
					coordDerivs[i][2] = function() { 
						return dj_du.apply(undefined, arguments) 
							* dk_dv.apply(undefined, arguments) 
							- dk_du.apply(undefined,arguments) 
							* dj_dv.apply(undefined,arguments);
					};
				});
			} catch (e) {
				failed = true;
			}
			if (failed) {
				//color the failed input textarea red
				$.each(allInputs, function(i,input) {
					input.css({background:'rgb(255,127,127)'});
				});
				return;
			} else {
				$.each(allInputs, function(i,input) {
					input.css({background:'white'});
				});
			}
			currentCoordChart.mapping = function(coord) {
				var result = [];
				for (var k = 0; k < coordLabels.length; ++k) {
					result[k] = coordCallbacks[k].apply(undefined, coord);
				}
				return result;
			};
			currentCoordChart.diff = function(coord, dim) {
				var result = [];
				for (var k = 0; k < coordLabels.length; ++k) {
					result[k] = coordDerivs[k][dim].apply(undefined, coord);
				}
				return result;
			};
			currentCoordChart.unitDiff = function(coord, dim) {
				var diff = this.diff(coord, dim);
				var len = vec3.length(diff);
				return [diff[0] / len, diff[1] / len, diff[2] / len];
			};
			currentCoordChart.connection = function(coord, k, j) {
				var result = [];
				for (var i = 0; i < parameters.length; ++i) {
					result[i] = conn2nd[i][j][k].apply(undefined, coord);
				}
				return result;
			};
			
			//process TeX
			var symbols = ['theta', 'phi', 'rho', 'gamma', 'Gamma'];
			$.each(symbols, function(i,sym) {
				propertiesHTML = propertiesHTML.replace(new RegExp(sym, 'g'), '\\'+sym);
			});
			$('#properties').html(propertiesHTML);
			MathJax.Hub.Queue(["Typeset", MathJax.Hub, "properties"]);
			
			reset();	//regenerate mesh
		};
		//wait for done lua loading
		var loadingInterval = setInterval(function() {
			if (luaDoneLoading) {
				clearInterval(loadingInterval);
				doUpdateEquations();
			}
		}, 100);
		
	};

	$.each(allInputs, function(i,input) {
		input
			.on('change', updateEquations)
			.on('keypress', updateEquations)
			.on('paste', updateEquations)
			.on('input', updateEquations);
	});

	$.each(coordCharts, function(name,coordChart) {
		var option = $('<option>', {
			text : name,
			value : name
		}).appendTo($('#coordinateSystem'));
		if (coordChart[name] == currentCoordChart) option.attr('selected', 'selected');
	});
	$('#coordinateSystem').change(function() {
		selectCoordChart($(this).val());
	});
	function selectCoordChart(name) {
		currentCoordChart = coordCharts[name];
		$.each(coordLabels, function(i,x) { $('#equation'+x.toUpperCase()).val('0'); });
		if (currentCoordChart.equations !== undefined) {
			$.each(currentCoordChart.equations, function(i,equation) {
				$('#equation'+coordLabels[i].toUpperCase()).val(equation);
			});
		}
		$('#parameters').val(currentCoordChart.parameters.join(', '));
		$('#constants').val('');
		if (currentCoordChart.constants !== undefined) {
			var constantText = [];
			$.each(currentCoordChart.constants, function(k,v) {
				constantText.push(k+'='+v);
			});
			if (constantText.length > 0) {
				$('#constants').val(constantText.concat());
			}
		}
		updateEquations();
	}

	var tmpQ = quat.create();	
	mouse = new Mouse3D({
		pressObj : canvas,
		move : function(dx,dy) {
			if (inputState == 'select') {
				var bestCoord = findClickedCoord();
			
				if (bestCoord !== undefined) {
					selectCoord(bestCoord);
				}
			}
			
			if (inputState == 'rotate') {
				var rotAngle = Math.PI / 180 * .03 * Math.sqrt(dx*dx + dy*dy);
				quat.setAxisAngle(tmpQ, [dy, dx, 0], rotAngle);

				quat.mul(sceneObj.angle, tmpQ, sceneObj.angle);
				quat.normalize(sceneObj.angle, sceneObj.angle);
			}

			if (inputState == 'direction') {
				//re-integrate geodesic and uncorrected coordinates along surface
				var clickedCoord = findClickedCoord();
				if (clickedCoord !== undefined) {
					var mappedCurrentCoord = currentCoordChart.mapping(currentCoord);
					//TODO
					//rotate view ray and origin by inverse angle
					//intersect view ray with plane of surface normal and mapped coordinate
					
					
					var mappedClickedCoord = currentCoordChart.mapping(clickedCoord);
					var delta = [
						mappedClickedCoord[0] - mappedCurrentCoord[0],
						mappedClickedCoord[1] - mappedCurrentCoord[1],
						mappedClickedCoord[2] - mappedCurrentCoord[2]
					];
					var invAngle = quat.create();
					quat.conjugate(invAngle, renderer.view.angle);
					vec3.transformQuat(delta, delta, invAngle);
					var basis0 = currentCoordChart.unitDiff(currentCoord, 0);
					var basis1 = currentCoordChart.unitDiff(currentCoord, 1);
					var d0 = vec3.dot(basis0, delta);
					var d1 = vec3.dot(basis1, delta);
					var dl = vec2.length([d0, d1]); 
					chooseDirection([d0/dl, d1/dl]);
				}
			}
		},
		zoom : function(dz) {
			renderer.view.pos[2] += .0001 * dz;
		}
	});

	var offset = [ [0,0], [1,0], [1,1], [1,1], [0,1], [0,0] ];
	var indexes = [];
	for (var itheta = 0; itheta < intDivs[0]-1; ++itheta) {
		for (var iphi = 0; iphi < intDivs[1]-1; ++iphi) {
			for (var k = 0; k < offset.length; ++k) {
				indexes.push(
					iphi + offset[k][0] + (itheta + offset[k][1]) * intDivs[1] 
				);
			}
		}
	}

	renderer.view.pos[2] = 2;
	renderer.view.zFar = 100;
	renderer.view.zNear = .1;
	var plainShader = new GL.ShaderProgram({
		context : gl,
		vertexPrecision : 'best',
		vertexCode : mlstr(function(){/*
attribute vec3 vertex;
uniform mat4 projMat;
uniform mat4 mvMat;
void main() {
	gl_Position = projMat * mvMat * vec4(vertex, 1.);
}
*/}),
		fragmentPrecision : 'best',
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
void main() {
	gl_FragColor = color;
}
*/})
	});
	var meshShader = new GL.ShaderProgram({
		context : gl,
		vertexPrecision : 'best',
		vertexCode : mlstr(function(){/*
attribute vec3 vertex;
attribute vec2 intCoord;
attribute vec3 normal;
uniform mat4 projMat;
uniform mat4 mvMat;
varying vec2 intCoordV;
varying vec3 normalV;
varying vec3 vertexV;
void main() {
	normalV = (mvMat * vec4(normal, 0.)).xyz;
	vec4 mvtx = mvMat * vec4(vertex, 1.);
	vertexV = mvtx.xyz;
	gl_Position = projMat * mvtx;
	intCoordV = intCoord;
}
*/}),
		fragmentPrecision : 'best',
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
varying vec2 intCoordV;
varying vec3 normalV;
varying vec3 vertexV;
void main() {
	vec3 n = normalize(normalV);
	if (n.z < 0.) n = -n;	//backface lighting
	
	vec2 fc = mod(intCoordV.xy / 10., 1.); //grid
	float i = 1. - 8. * fc.x * fc.y * (1. - fc.x) * (1. - fc.y);
	i = pow(i, 50.);
	
	gl_FragColor = vec4(.25, .5, .5, 1.);
	gl_FragColor.rgb *= 1. - i;
	vec3 u = normalize(vertexV);
	float l = dot(n, u);
	gl_FragColor.rgb *= max(abs(l), .3);
}
*/})
	});
	sceneObj = new GL.SceneObject({
		context : gl,
		scene : renderer.scene,
		static : false
	});

	meshObj = new GL.SceneObject({
		context : gl,
		scene : renderer.scene,
		parent : sceneObj,
		mode : gl.TRIANGLES,
		attrs : {
			vertex : new GL.ArrayBuffer({
				context : gl,
				count : intDivs[0] * intDivs[1],
				usgae : gl.DYNAMIC_DRAW
			}),
			coord : new GL.ArrayBuffer({
				context : gl,
				dim : 2,
				count : intDivs[0] * intDivs[1]
			}),
			intCoord : new GL.ArrayBuffer({
				context : gl,
				dim : 2,
				count : intDivs[0] * intDivs[1],
				usage : gl.DYNAMIC_DRAW
			}),
			normal : new GL.ArrayBuffer({
				context : gl,
				count : intDivs[0] * intDivs[1],
				usage : gl.DYNAMIC_DRAW
			})
		},
		indexes : new GL.ElementArrayBuffer({
			context : gl,
			data : new Uint16Array(indexes)
		}),
		uniforms : {
			color : [1,1,1,1]
		},
		shader : meshShader,
		static : false
	});
	
	selectionObj = new GL.SceneObject({
		context : gl,
		scene : renderer.scene,
		parent : sceneObj,
		mode : gl.LINE_LOOP,
		attrs : {
			vertex : new GL.ArrayBuffer({
				context : gl,
				usage : gl.DYNAMIC_DRAW,
				count : selectionRes
			})
		},
		uniforms : {
			color : [1,1,1,1]
		},
		shader : plainShader,
		static : false
	});

	var basisColors = [
		[1,0,0,1],
		[0,1,0,1]
	];
	var connColors = [
		[
			[1,1,0,1],
			[0,1,1,1]
		],
		[
			[1,0,1,1],
			[1,.5,0,1]
		]
	];
	for (var i = 0; i < 2; ++i) {
		basisObjs[i] = new GL.SceneObject({
			context : gl,
			scene : renderer.scene,
			parent : sceneObj,
			mode : gl.LINES,
			attrs : {
				vertex : new GL.ArrayBuffer({
					context : gl,
					count : 2,
					usage : gl.DYNAMIC_DRAW
				})
			},
			uniforms : {
				color : basisColors[i]
			},
			shader : plainShader,
			static : false
		});
		/*
		connObjs[i] = [];
		for (var j = 0; j < 2; ++j) {
			connObjs[i][j] = new GL.SceneObject({
				context : gl,
				scene : renderer.scene,
				parent : sceneObj,
				mode : gl.LINES,
				attrs : {
					vertex : new GL.ArrayBuffer({
						context : gl,
						count : 2,
						usage : gl.DYNAMIC_DRAW
					})
				},
				uniforms : {
					color : connColors[i][j]
				},
				shader : plainShader,
				static : false
			});
		}
		*/
	}

	partialPathObj = new GL.SceneObject({
		context : gl,
		scene : renderer.scene,
		parent : sceneObj,
		mode : gl.LINE_STRIP,
		attrs : {
			vertex : new GL.ArrayBuffer({
				context : gl,
				count : pathLength,
				usage : gl.DYNAMIC_DRAW
			})
		},
		uniforms : {
			color : [1,1,0,1]
		},
		shader : plainShader,
		static : false
	});
	geodesicPathObj = new GL.SceneObject({
		context : gl,
		scene : renderer.scene,
		parent : sceneObj,
		mode : gl.LINE_STRIP,
		attrs : {
			vertex : new GL.ArrayBuffer({
				context : gl,
				count : pathLength,
				usage : gl.DYNAMIC_DRAW
			})
		},
		uniforms : {
			color : [1,0,1,1]
		},
		shader : plainShader,
		static : false
	});

	gl.enable(gl.DEPTH_TEST);

	selectCoordChart('Spherical');	
	
	$(window).resize(resize);
	resize();
	
	update();
});

