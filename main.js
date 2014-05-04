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
	coords = x_j
	returns y_i

diff:
	dy_i/dx_j
	 dim = j
	 coords = x_j
	 return dy_i/dx_j

unitDiff: diff, normalized algebraically to remove singularities
*/
var coordCharts = {
	Spherical : {
		mapping : function(coords) {
			var r = coords[0];
			var theta = coords[1];
			var phi = coords[2];
			var x = r * Math.sin(theta) * Math.cos(phi);
			var y = r * Math.sin(theta) * Math.sin(phi);
			var z = r * Math.cos(theta);
			return [x, y, z];
		},
		diff : function(coords, dim) {
			var r = coords[0];
			var theta = coords[1];
			var phi = coords[2];
			switch (dim) {
			case 0:
				return [
					Math.sin(theta) * Math.cos(phi),
					Math.sin(theta) * Math.sin(phi),
					Math.cos(theta)];
			case 1:
				return [
					r * Math.cos(theta) * Math.cos(phi),
					r * Math.cos(theta) * Math.sin(phi),
					-r * Math.sin(theta)];
			case 2:
				return [
					-r * Math.sin(theta) * Math.sin(phi),
					r * Math.sin(theta) * Math.cos(phi),
					0];
			}
		},
		unitDiff : function(coords, dim) {
			var r = coords[0];
			var theta = coords[1];
			var phi = coords[2];
			switch (dim) {
			case 0:
				return [
					Math.sin(theta) * Math.cos(phi),
					Math.sin(theta) * Math.sin(phi),
					Math.cos(theta)];
			case 1:
				return [
					Math.cos(theta) * Math.cos(phi),
					Math.cos(theta) * Math.sin(phi),
					-Math.sin(theta)];
			case 2:
				return [
					-Math.sin(phi),
					Math.cos(phi),
					0];
			}
		}
	}
};
var currentCoordChart;

var selectionObj;
var selectionRes = 20;
var basisObjs = [];
function selectCoord(coords) {
	var basis0 = currentCoordChart.unitDiff(coords, 1);
	var basis1 = currentCoordChart.unitDiff(coords, 2);
	var basis = [basis0, basis1];
	var mappedCoords = currentCoordChart.mapping(coords);

	var radius = .1;
	var modv = vec3.create();
	for (var i = 0; i < selectionRes; ++i) {
		var psi = i/selectionRes * Math.PI * 2;
		vec3.scaleAndAdd(modv, mappedCoords, basis0, radius * Math.cos(psi));
		vec3.scaleAndAdd(modv, modv, basis1, radius * Math.sin(psi));
		for (var j = 0; j < 3; ++j) {
			selectionObj.attrs.vertex.data[3 * i + j] = modv[j];
		}
	}
	selectionObj.attrs.vertex.updateData();

	$.each(basisObjs, function(i,basisObj) {
		var e = basis[i];
		for (var j = 0; j < 3; ++j) {
			basisObj.attrs.vertex.data[j] = mappedCoords[j];
			basisObj.attrs.vertex.data[3 + j] = mappedCoords[j] + .4 * e[j];
		}
		basisObj.attrs.vertex.updateData();
	});
}

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

	$('#tools_rotate').click(function() { dragging = true; });
	$('#tools_select').click(function() { dragging = false; });

	currentCoordChart = coordCharts.Spherical;
	$('#coord_Spherical').attr('checked', 'checked');

	window.dragging = true;	//GUI-ize me plz
	var tmp = vec3.create();
	var tmp2 = vec3.create();
	var tmpQ = quat.create();	
	mouse = new Mouse3D({
		pressObj : canvas,
		move : function(dx,dy) {
			if (!dragging) {
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
				for (var i = 0; i < meshMapping.length; ++i) {
					vec3.transformQuat(pt, meshMapping[i].dst, sceneObj.angle);
					
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
						bestCoord = meshMapping[i].src;
					}
				}
				
				if (bestDist !== undefined) {
					selectCoord(bestCoord);
				}
			}
			
			if (dragging) {
				var rotAngle = Math.PI / 180 * .03 * Math.sqrt(dx*dx + dy*dy);
				quat.setAxisAngle(tmpQ, [dy, dx, 0], rotAngle);

				quat.mul(sceneObj.angle, tmpQ, sceneObj.angle);
				quat.normalize(sceneObj.angle, sceneObj.angle);
			}
		},
		zoom : function(dz) {
			GL.view.fovY *= Math.exp(-.0003 * dz);
			GL.view.fovY = Math.clamp(GL.view.fovY, 1, 179);
			GL.updateProjection();
		}
	});

	//generate metric geometry
	var rMin = 1;
	var rMax = 1;
	var rDiv = 1;
	var thetaMin = 0;
	var thetaMax = Math.PI;
	var thetaDiv = 6;
	var phiMin = -Math.PI;
	var phiMax = Math.PI;
	var phiDiv = 12;
	var lines = [];
	var intCoordsToCoords = function(intCoords) {
		//intCoords
		var ir = intCoords[0];
		var itheta = intCoords[1];
		var iphi = intCoords[2];
		//converted to coordinate chart coordinates
		var r = (ir / rDiv) * (rMax - rMin) + rMin;
		var theta = (itheta / thetaDiv) * (thetaMax - thetaMin) + thetaMin;
		var phi = (iphi / phiDiv) * (phiMax - phiMin) + phiMin;
		return [r, theta, phi];
	};
	meshMapping = [];
	vertexes = [];
	var indexes = [];
	//TODO line strip collections
	var addLine = function(ia, ib) {
		var res = 100;
		for (var j = 0; j < res; ++j) { 
			var f = j / res;
			var nf = (j + 1) / res;
			var iaf = vec3.create();
			vec3.lerp(iaf, ia, ib, f);
			var ibf = vec3.create();
			vec3.lerp(ibf, ia, ib, nf);
			var ca = intCoordsToCoords(iaf);
			var cb = intCoordsToCoords(ibf);
			var va = currentCoordChart.mapping(ca);
			var vb = currentCoordChart.mapping(cb);
			for (var i = 0; i < 3; ++i) {
				vertexes.push(va[i]);
			}
			meshMapping.push({
				src : ca.clone(),
				dst : va.clone()
			});
			indexes.push(vertexes.length/3-1);
			for (var i = 0; i < 3; ++i) {
				vertexes.push(vb[i]);
			}
			indexes.push(vertexes.length/3-1);
			meshMapping.push({
				src : cb.clone(),
				dst : vb.clone()
			});
		}
	};
	var intDivs = [rDiv, thetaDiv, phiDiv];
	for (ir = 0; ir < rDiv; ++ir) {
		for (var itheta = 0; itheta < thetaDiv; ++itheta) {
			for (var iphi = 0; iphi < phiDiv; ++iphi) {
				var intCoords = [ir, itheta, iphi];
				for (var i = 0; i < 3; ++i) {
					if (!(i == 0 && intCoords[i] >= intDivs[i]-1)) {
						var ia = [ir, itheta, iphi];
						var ib = ia.clone();
						++ib[i];
						addLine(ia, ib);
					}
				}
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

	sceneObj = new GL.SceneObject({
		static : false
	});

	var meshObj = new GL.SceneObject({
		parent : sceneObj,
		mode : gl.LINES,
		attrs : {
			vertex : new GL.ArrayBuffer({
				data : new Float32Array(vertexes)
			})
		},
		//indexes : indexes,
		uniforms : {
			color : [1,1,1,1]
		},
		shader : plainShader,
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
			color : [0,1,1,1]
		},
		shader : plainShader,
		static : false
	});

	var basisColors = [
		[1,0,0,1],
		[0,1,0,1]
	];
	for (var i = 0; i < 2; ++i) {
		basisObjs.push(new GL.SceneObject({
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
		}));
	}

	selectCoord([1,0,0]);

	$(window).resize(resize);
	resize();
	
	update();
});

