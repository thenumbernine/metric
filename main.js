var canvas;
var gl;
var mouse;
var sceneObj;

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	GL.resize();

	var info = $('#info');
	var width = window.innerWidth 
		- parseInt(info.css('padding-left'))
		- parseInt(info.css('padding-right'));
	info.width(width);
	var height = window.innerHeight
		- parseInt(info.css('padding-top'))
		- parseInt(info.css('padding-bottom'));
	info.height(height - 32);
}

function update() {
	//update
	GL.draw();
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
var r = 1;
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
		mapping : function(coord) {
			var theta = coord[0];
			var phi = coord[1];
			var x = r * Math.sin(theta) * Math.cos(phi);
			var y = r * Math.sin(theta) * Math.sin(phi);
			var z = r * Math.cos(theta);
			return [x, y, z];
		},
		diff : function(coord, dim) {
			var theta = coord[0];
			var phi = coord[1];
			switch (dim) {
			case 0:
				return [
					r * Math.cos(theta) * Math.cos(phi),
					r * Math.cos(theta) * Math.sin(phi),
					-r * Math.sin(theta)];
			case 1:
				return [
					-r * Math.sin(theta) * Math.sin(phi),
					r * Math.sin(theta) * Math.cos(phi),
					0];
			case 2:	//normal
				return [
					Math.sin(theta) * Math.cos(phi),
					Math.sin(theta) * Math.sin(phi),
					Math.cos(theta)];
			}
		},
		unitDiff : function(coord, dim) {
			var theta = coord[0];
			var phi = coord[1];
			switch (dim) {
			case 0:
				return [
					Math.cos(theta) * Math.cos(phi),
					Math.cos(theta) * Math.sin(phi),
					-Math.sin(theta)];
			case 1:
				return [
					-Math.sin(phi),
					Math.cos(phi),
					0];
			case 2:
				return [
					Math.sin(theta) * Math.cos(phi),
					Math.sin(theta) * Math.sin(phi),
					Math.cos(theta)];
			}
		},
		// delta_k e_j = Gamma^i_jk e_i
		// maybe I should add the extrinsic components ... 
		//connection(coord, k,j)[i] == Gamma^i_jk
		connection : function(coord, basisIndex, wrtIndex) {
			var theta = coord[0];
			var phi = coord[1];
			return [
				//deriv of e_theta
				[
					//wrt e_theta
					[
						0, //conn^theta_theta_theta,
						0, //conn^phi_theta_theta
					],
					//wrt e_phi
					[
						0, //conn^theta_theta_phi
						1/Math.tan(theta), //conn^phi_theta_phi
					],
				],
				//deriv of e_phi
				[
					//wrt e_theta
					[
						0, //conn^theta_phi_theta
						1/Math.tan(theta), //conn^phi_phi_theta
					],
					//wrt e_phi
					[
						-.5 * Math.sin(2 * theta), //conn^theta_phi_phi
						0, //conn^phi_phi_phi
					],
				],
			][basisIndex][wrtIndex];	//one more reason why I hate javascript: array construction and dereferencing syntax
		}
	},
	Polar : {
		mapping : function(coord) {
			var r = coord[0];
			var phi = coord[1];
			return [r * Math.cos(phi), r * Math.sin(phi), 0];
		},
		diff : function(coord, dim) {
			var r = coord[0];
			var phi = coord[1];
			switch (dim) {
			case 0: return [Math.cos(phi), Math.sin(phi), 0];
			case 1: return [-r * Math.sin(phi), r * Math.cos(phi), 0];
			case 2: return [0,0,1];
			}
		},
		unitDiff : function(coord, dim) {
			var r = coord[0];
			var phi = coord[1];
			switch (dim) {
			case 0: return [Math.cos(phi), Math.sin(phi), 0];
			case 1: return [-Math.sin(phi), Math.cos(phi), 0];
			case 2: return [0,0,1];
			}
		}
	}
};
var currentCoordChart;

var selectionObj;
var selectionRes = 20;
var basisObjs = [];
//var connObjs = [];
var partialPathObj;
var geodesicPathObj;
var pathLength = 200;
var currentCoord = [0.26623666555845704, 1.8215957403167709];
var currentDirection = [0, 1];
function selectCoord(coord) {
	currentCoord = coord.clone();
	var basis0 = currentCoordChart.diff(coord, 0);
	var basis1 = currentCoordChart.diff(coord, 1);
	var basis = [basis0, basis1];
	
	var unitBasis0 = currentCoordChart.unitDiff(coord, 0);
	var unitBasis1 = currentCoordChart.unitDiff(coord, 1);
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
			basisObj.attrs.vertex.data[3 + k] = mappedCoord[k] + basisLength * basis[i][k];
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
	var step = .025;
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

var vertexArray = [];
var intCoordArray = [];
var coordArray = [];
var normalArray = [];


var findClickedCoord;
(function(){
	var tmp = vec3.create();
	var tmp2 = vec3.create();
	findClickedCoord = function() {
		var mouseDir = vec3.create();
		GL.mouseDir(mouseDir, mouse.xf, mouse.yf);
		vec3.normalize(mouseDir, mouseDir);
		//ray intersection test with coordinate chart ...
		// for now just search for closest point in geometry?
		var bestDist = undefined;
		var bestCoord = undefined;
		vec3.quatZAxis(tmp2, GL.view.angle);
		var considered = 0;
		var pt = vec3.create();
		for (var i = 0; i < vertexArray.length/3; ++i) {
			var normal = [
				normalArray[3*i+0],
				normalArray[3*i+1],
				normalArray[3*i+2]
			];
			vec3.transformQuat(normal, normal, sceneObj.angle);
			if (vec3.dot(normal, tmp2) < 0) continue;	//forward dot normal > 0 means the surface is back-facing
			
			var vertex = [
				vertexArray[3*i+0],
				vertexArray[3*i+1],
				vertexArray[3*i+2]
			];
			vec3.transformQuat(pt, vertex, sceneObj.angle);

			
			//make sure we're on the right side of the view plane
			vec3.sub(tmp, pt, GL.view.pos);
			if (vec3.dot(tmp, tmp2) > 0) continue;	//fwd dot delta > 0 means we're good, so -fwd dot delta < 0 means we're good, so -fwd dot delta > 0 means we're bad
		
			considered++;

			//ray/point distance from view pos / mouse line
			// to meshMapping[i].dst
			vec3.sub(tmp, pt, GL.view.pos);
			vec3.cross(tmp, tmp, mouseDir);
			var dist = vec3.length(tmp);
			
			if (bestDist === undefined || dist < bestDist) {
				bestDist = dist;
				bestCoord = [
					coordArray[3*i+0],
					coordArray[3*i+1],
					coordArray[3*i+2]
				];
			}
		}
		return bestCoord;
	}
})();

$(document).ready(function() {
	$('#panelButton').click(function() {
		var panel = $('#panel');	
		if (panel.css('display') == 'none') {
			panel.show();
			$('#info').hide();
		} else {
			panel.hide();
		}
	});
	$('#infoButton').click(function() {
		var info = $('#info');
		if (info.css('display') == 'none') {
			info.show();
			$('#panel').hide();
		} else {
			info.hide();
		}
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
		gl = GL.init(canvas);
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
	$('#coord_Spherical').attr('checked', 'checked');
	
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
					var mappedClickedCoord = currentCoordChart.mapping(clickedCoord);
					var mappedCurrentCoord = currentCoordChart.mapping(currentCoord);
					var delta = [
						mappedClickedCoord[0] - mappedCurrentCoord[0],
						mappedClickedCoord[1] - mappedCurrentCoord[1],
						mappedClickedCoord[2] - mappedCurrentCoord[2]
					];
					var invAngle = quat.create();
					quat.conjugate(invAngle, GL.view.angle);
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
			GL.view.fovY *= Math.exp(-.0003 * dz);
			GL.view.fovY = Math.clamp(GL.view.fovY, 1, 179);
			GL.updateProjection();
		}
	});

	//generate metric geometry
	var thetaMin = 0;
	var thetaMax = Math.PI;
	var thetaDiv = 60;
	var phiMin = -Math.PI;
	var phiMax = Math.PI;
	var phiDiv = 120;
	var lines = [];
	var intCoordToCoord = function(intCoord) {
		//intCoord
		var itheta = intCoord[0];
		var iphi = intCoord[1];
		//converted to coordinate chart coordinates
		var theta = (itheta / (thetaDiv-1)) * (thetaMax - thetaMin) + thetaMin;
		var phi = (iphi / (phiDiv-1)) * (phiMax - phiMin) + phiMin;
		return [theta, phi];
	};
	meshMapping = [];
	var intDivs = [thetaDiv, phiDiv];
	for (var itheta = 0; itheta < thetaDiv; ++itheta) {
		for (var iphi = 0; iphi < phiDiv; ++iphi) {
			var intCoord = [itheta, iphi];
			var coord = intCoordToCoord(intCoord);
			var mappedCoord = currentCoordChart.mapping(coord);
			var normal = currentCoordChart.unitDiff(coord, 2);
			for (var k = 0; k < 2; ++k) {
				intCoordArray.push(intCoord[k]);
			}
			for (var k = 0; k < 3; ++k) {
				coordArray.push(coord[k]);
			}
			for (var k = 0; k < 3; ++k) {
				vertexArray.push(mappedCoord[k]);
			}
			for (var k = 0; k < 3; ++k) {
				normalArray.push(normal[k]);
			}
		}
	}
	var offset = [ [0,0], [1,0], [1,1], [1,1], [0,1], [0,0] ];
	var indexes = [];
	for (var itheta = 0; itheta < thetaDiv-1; ++itheta) {
		for (var iphi = 0; iphi < phiDiv-1; ++iphi) {
			for (var k = 0; k < offset.length; ++k) {
				indexes.push(
					iphi + offset[k][0] + (itheta + offset[k][1]) * phiDiv 
				);
			}
		}
	}

	GL.view.pos[2] = 2;
	GL.view.zFar = 100;
	GL.view.zNear = .1;
	var plainShader = new GL.ShaderProgram({
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
		vertexPrecision : 'best',
		vertexCode : mlstr(function(){/*
attribute vec3 vertex;
attribute vec2 intCoord;
attribute vec3 normal;
uniform mat4 projMat;
uniform mat4 mvMat;
varying vec2 intCoordV;
varying vec3 normalV;
void main() {
	normalV = (mvMat * vec4(normal, 0.)).xyz;
	gl_Position = projMat * mvMat * vec4(vertex, 1.);
	intCoordV = intCoord;
}
*/}),
		fragmentPrecision : 'best',
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
varying vec2 intCoordV;
varying vec3 normalV;
void main() {
	vec3 n = normalize(normalV);
	if (n.z < 0.) n = -n;	//backface
	vec2 fc = mod(intCoordV.xy / 10., 1.); 
	float i = 1. - 8. * fc.x * fc.y * (1. - fc.x) * (1. - fc.y);
	i = pow(i, 50.);
	gl_FragColor = vec4(.25, .5, .5, 1.);
	gl_FragColor.rgb *= 1. - i;
	gl_FragColor.rgb *= n.z;
}
*/})
	});
	sceneObj = new GL.SceneObject({
		static : false
	});

	var meshObj = new GL.SceneObject({
		parent : sceneObj,
		mode : gl.TRIANGLES,
		attrs : {
			vertex : new GL.ArrayBuffer({
				data : new Float32Array(vertexArray)
			}),
			intCoord : new GL.ArrayBuffer({
				dim : 2,
				data : new Float32Array(intCoordArray)
			}),
			normal : new GL.ArrayBuffer({
				data : new Float32Array(normalArray)
			})
		},
		indexes : new GL.ElementArrayBuffer({
			data : new Uint16Array(indexes)
		}),
		uniforms : {
			color : [1,1,1,1]
		},
		shader : meshShader,
		static : false
	});

	selectionObj = new GL.SceneObject({
		parent : sceneObj,
		mode : gl.LINE_LOOP,
		attrs : {
			vertex : new GL.ArrayBuffer({
				usage : gl.DYNAMIC_DRAW,
				count : selectionRes,
				keep : true
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
			parent : sceneObj,
			mode : gl.LINES,
			attrs : {
				vertex : new GL.ArrayBuffer({
					count : 2,
					usage : gl.DYNAMIC_DRAW,
					keep : true
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
				parent : sceneObj,
				mode : gl.LINES,
				attrs : {
					vertex : new GL.ArrayBuffer({
						count : 2,
						usage : gl.DYNAMIC_DRAW,
						keep : true
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
		parent : sceneObj,
		mode : gl.LINE_STRIP,
		attrs : {
			vertex : new GL.ArrayBuffer({
				count : pathLength,
				usage : gl.DYNAMIC_DRAW,
				keep : true
			})
		},
		uniforms : {
			color : [1,1,0,1]
		},
		shader : plainShader,
		static : false
	});
	geodesicPathObj = new GL.SceneObject({
		parent : sceneObj,
		mode : gl.LINE_STRIP,
		attrs : {
			vertex : new GL.ArrayBuffer({
				count : pathLength,
				usage : gl.DYNAMIC_DRAW,
				keep : true
			})
		},
		uniforms : {
			color : [1,0,1,1]
		},
		shader : plainShader,
		static : false
	});

	gl.enable(gl.DEPTH_TEST);

	selectCoord(currentCoord);
	chooseDirection(currentDirection);

	$(window).resize(resize);
	resize();
	
	update();
});

