import {vec2, vec3, quat} from '/js/gl-matrix-3.4.1/index.js';
import {DOM, getIDs, removeFromParent, show, hide, arrayClone} from '/js/util.js';
import {EmbeddedLuaInterpreter} from '/js/lua.vm-util.js';
import {GLUtil, quatZAxis} from '/js/gl-util.js';
import {Mouse3D} from '/js/mouse3d.js';
import {makeGradient} from '/js/gl-util-Gradient.js';
import {makeUnitQuad} from '/js/gl-util-UnitQuad.js';
import {makeFloatTexture2D} from '/js/gl-util-FloatTexture2D.js';

const ids = getIDs();
window.ids = ids;
const urlparams = new URLSearchParams(window.location.search);

let canvas;
let gl;
let glutil;
let mouse;
let sceneObj;
let meshObj;

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	glutil.resize();

	const width = window.innerWidth
		- parseInt(ids.info.style.paddingLeft)
		- parseInt(ids.info.style.paddingRight);
	ids.info.style.width = width+'px';
	const height = window.innerHeight
		- parseInt(ids.info.style.paddingTop)
		- parseInt(ids.info.style.paddingBottom);
	ids.info.style.height = (height - 32)+'px';
	ids.panel.style.height = (height - 16)+'px';
}

function update() {
	//update
	glutil.draw();
	requestAnimationFrame(update);
};

// https://docs.mathjax.org/en/latest/web/typeset.html#typeset-async
// new MathJax is a bit more restrictive of how to handle concurrent rendering ...
function typeset(code) {
	MathJax.startup.promise = MathJax.startup.promise
		.then(() => MathJax.typesetPromise(code()))
		.catch((err) => console.log('Typeset failed: ' + err.message));
	return MathJax.startup.promise;
}

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
let coordCharts = {
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
let currentCoordChart;

let intDivs = [120,120];

let intCoordToCoord = function(intCoord) {
	let coord = [];
	let coordMin = currentCoordChart.coordinateMin;
	let coordMax = currentCoordChart.coordinateMax;
	for (let k = 0; k < intCoord.length; ++k) {
		coord[k] = intCoord[k] / (intDivs[k]-1) * (coordMax[k] - coordMin[k]) + coordMin[k];
	}
	return coord;
};

function reset() {
	//generate metric geometry
	let thetaDiv = intDivs[0];
	let phiDiv = intDivs[1];
	for (let itheta = 0; itheta < thetaDiv; ++itheta) {
		for (let iphi = 0; iphi < phiDiv; ++iphi) {
			let intCoord = [itheta, iphi];
			let coord = intCoordToCoord(intCoord);
			let mappedCoord = currentCoordChart.mapping(coord);
			let normal = currentCoordChart.unitDiff(coord, 2);
			for (let k = 0; k < 2; ++k) {
				meshObj.attrs.intCoord.buffer.data[k + 2 * (iphi + phiDiv * itheta)] = intCoord[k];
			}
			for (let k = 0; k < 2; ++k) {
				meshObj.attrs.coord.buffer.data[k + 2 * (iphi + phiDiv * itheta)] = coord[k];
			}
			for (let k = 0; k < 3; ++k) {
				meshObj.attrs.vertex.buffer.data[k + 3 * (iphi + phiDiv * itheta)] = mappedCoord[k];
			}
			for (let k = 0; k < 3; ++k) {
				meshObj.attrs.normal.buffer.data[k + 3 * (iphi + phiDiv * itheta)] = normal[k];
			}
		}
	}
	meshObj.attrs.intCoord.buffer.updateData();
	meshObj.attrs.coord.buffer.updateData();
	meshObj.attrs.vertex.buffer.updateData();
	meshObj.attrs.normal.buffer.updateData();

	currentCoord = arrayClone(currentCoordChart.initialCoord);
	currentDirection = arrayClone(currentCoordChart.initialDirection);
	//TODO mesh regen as well since we will want varying resolutions ... or not?
	selectCoord(currentCoord);
}

let selectionObj;
let selectionRes = 100;
let basisObjs = [];
//let connObjs = [];
let partialPathObj;
let geodesicPathObj;
let pathLength = 5000;
let currentCoord = [0.26623666555845704, 1.8215957403167709];
let currentDirection = [0, 1];
function selectCoord(coord) {
	currentCoord = arrayClone(coord);
	let basis0 = currentCoordChart.diff(coord, 0);
	let basis1 = currentCoordChart.diff(coord, 1);
	let basis = [basis0, basis1];

	let unitBasis0 = currentCoordChart.unitDiff(coord, 0);
	let unitBasis1 = currentCoordChart.unitDiff(coord, 1);
	let normal = currentCoordChart.unitDiff(coord,2);
	let unitBasis = [unitBasis0, unitBasis1];

	let mappedCoord = currentCoordChart.mapping(coord);

	let selectionRadius = .1;
	let modv = vec3.create();
	for (let i = 0; i < selectionRes; ++i) {
		let psi = i/selectionRes * Math.PI * 2;
		vec3.scaleAndAdd(modv, mappedCoord, unitBasis0, selectionRadius * Math.cos(psi));
		vec3.scaleAndAdd(modv, modv, unitBasis1, selectionRadius * Math.sin(psi));
		for (let j = 0; j < 3; ++j) {
			selectionObj.attrs.vertex.buffer.data[3 * i + j] = modv[j];
		}
	}
	selectionObj.attrs.vertex.buffer.updateData();

	let basisLength = .4;
	for (let i = 0; i < 2; ++i) {
		let basisObj = basisObjs[i];
		for (let k = 0; k < 3; ++k) {
			basisObj.attrs.vertex.buffer.data[k] = mappedCoord[k];
			basisObj.attrs.vertex.buffer.data[3 + k] = mappedCoord[k] + basisLength * basis[i][k] + .01 * normal[k];
		}
		basisObj.attrs.vertex.buffer.updateData();
		/* need a better way to visualize this ...
		for (let j = 0; j < 2; ++j) {
			let conn = currentCoordChart.connection(coord, i, j);
			let connObj = connObjs[i][j];
			for (let k = 0; k < 3; ++k) {
				connObj.attrs.vertex.buffer.data[k] = mappedCoord[k] + basisLength * basis[i][k];
				connObj.attrs.vertex.buffer.data[3+k] = mappedCoord[k] + basisLength * (basis[i][k] + conn[0] * basis[0][k] + conn[1] * basis[1][k]);
			}
			connObj.attrs.vertex.buffer.updateData();
		}
		*/
	}

	chooseDirection(currentDirection);
}

function chooseDirection(direction) {
	currentDirection = arrayClone(direction);
	let partialCoord = arrayClone(currentCoord);
	let partialDir = arrayClone(direction);
	let geoCoord = arrayClone(currentCoord);
	let geoDir = arrayClone(direction);
	let step = .001;
	for (let i = 0; i < pathLength; ++i) {
		partialCoord[0] += step * partialDir[0];
		partialCoord[1] += step * partialDir[1];
		let mapped = currentCoordChart.mapping(partialCoord);
		let normal = currentCoordChart.unitDiff(partialCoord,2);
		partialPathObj.attrs.vertex.buffer.data[3*i+0] = mapped[0] + .01 * normal[0];
		partialPathObj.attrs.vertex.buffer.data[3*i+1] = mapped[1] + .01 * normal[1];
		partialPathObj.attrs.vertex.buffer.data[3*i+2] = mapped[2] + .01 * normal[2];
	}
	partialPathObj.attrs.vertex.buffer.updateData();
	for (let i = 0; i < pathLength; ++i) {
		geoCoord[0] += step * geoDir[0];
		geoCoord[1] += step * geoDir[1];
		let geoAccel = [0,0];
		for (let j = 0; j < 2; ++j) {
			for (let k = 0; k < 2; ++k) {
				for (let l = 0; l < 2; ++l) {
					geoAccel[j] -= currentCoordChart.connection(geoCoord,k,l)[j] * geoDir[k] * geoDir[l];
				}
			}
		}
		geoDir[0] += step * geoAccel[0];
		geoDir[1] += step * geoAccel[1];
		let mapped = currentCoordChart.mapping(geoCoord);
		let normal = currentCoordChart.unitDiff(geoCoord,2);
		geodesicPathObj.attrs.vertex.buffer.data[3*i+0] = mapped[0] + .01 * normal[0];
		geodesicPathObj.attrs.vertex.buffer.data[3*i+1] = mapped[1] + .01 * normal[1];
		geodesicPathObj.attrs.vertex.buffer.data[3*i+2] = mapped[2] + .01 * normal[2];
	}
	geodesicPathObj.attrs.vertex.buffer.updateData();
}

let findClickedCoord;
{
	const tmp = vec3.create();
	const tmp2 = vec3.create();
	findClickedCoord = function() {
		const mouseDir = vec3.create();
		glutil.mouseDir(mouseDir, mouse.xf, mouse.yf);
		vec3.normalize(mouseDir, mouseDir);
		//ray intersection test with coordinate chart ...
		// for now just search for closest point in geometry?
		let bestDist = undefined;
		let bestCoord = undefined;
		quatZAxis(tmp2, glutil.view.angle);
		let considered = 0;
		const pt = vec3.create();
		for (let i = 0; i < meshObj.attrs.vertex.buffer.data.length/3; ++i) {
			const normal = [
				meshObj.attrs.normal.buffer.data[3*i+0],
				meshObj.attrs.normal.buffer.data[3*i+1],
				meshObj.attrs.normal.buffer.data[3*i+2]
			];
			vec3.transformQuat(normal, normal, sceneObj.angle);
			if (vec3.dot(normal, tmp2) < 0) continue;	//forward dot normal > 0 means the surface is back-facing

			const vertex = [
				meshObj.attrs.vertex.buffer.data[3*i+0],
				meshObj.attrs.vertex.buffer.data[3*i+1],
				meshObj.attrs.vertex.buffer.data[3*i+2]
			];
			vec3.transformQuat(pt, vertex, sceneObj.angle);


			//make sure we're on the right side of the view plane
			vec3.sub(tmp, pt, glutil.view.pos);
			if (vec3.dot(tmp, tmp2) > 0) continue;	//fwd dot delta > 0 means we're good, so -fwd dot delta < 0 means we're good, so -fwd dot delta > 0 means we're bad

			considered++;

			//ray/point distance from view pos / mouse line
			vec3.sub(tmp, pt, glutil.view.pos);
			vec3.cross(tmp, tmp, mouseDir);
			const dist = vec3.length(tmp);

			if (bestDist === undefined || dist < bestDist) {
				bestDist = dist;
				bestCoord = [
					meshObj.attrs.coord.buffer.data[2*i+0],
					meshObj.attrs.coord.buffer.data[2*i+1]
				];
			}
		}
		return bestCoord;
	}
}

ids.panelButton.addEventListener('click', e => {
	show(ids.panel);
	hide(ids.info);
});
ids.infoButton.addEventListener('click', e => {
	show(ids.info);
	hide(ids.panel);
});

canvas = DOM('canvas', {
	css : {
		left : 0,
		top : 0,
		position : 'absolute',
		userSelect : 'none',
	},
	prependTo : document.body,
});
window.canvas = canvas;

try {
	glutil = new GLUtil({canvas:canvas});
	gl = glutil.context;
} catch (e) {
	console.log('glutil failed', e);
	removeFromParent(canvas);
	show(ids.webglfail);
	throw e;
}
show(ids.menu);
show(ids.panel);

if (urlparams.get('info')) {
	show(ids.info);
	hide(ids.panel);
}

window.inputState = 'rotate';

ids.tools_rotate.addEventListener('click', e => { inputState = 'rotate'; });
ids.tools_select.addEventListener('click', e => { inputState = 'select'; });
ids.tools_direction.addEventListener('click', e => { inputState = 'direction'; });

currentCoordChart = coordCharts.Spherical;




//init lua
let luaDoneLoading = false;
const lua = new EmbeddedLuaInterpreter({
	packages : ['ext', 'symmath', 'complex'],
	packageTests : ['symmath'],
	done : function() {
console.log('loaded lua');
		const thiz = this;
		this.capture({
			callback : () => {
				thiz.executeAndPrint(`
local symmath = require 'symmath'
symmath.setup()
local LaTeX = symmath.export.LaTeX
symmath.tostring = LaTeX
LaTeX.openSymbol = ''
LaTeX.closeSymbol = ''
`);
console.log('initialized symmath');
			},
			output : s => {
console.log(s);
			}
		});
		luaDoneLoading = true;
	},
	autoLaunch : true
});
window.lua = lua;
const coordLabels = ['x', 'y', 'z'];

const allInputs = coordLabels
.map(x => { return 'equation'+x.toUpperCase(); })
.concat(['parameters', 'constants'])
.map(id => { return ids[id]; });

//args: done: what to execute next
const updateEquations = () => {
	const doUpdateEquations = () => {
		//declare parameter variables
		const parameters = ids.parameters.value.split(',').map(s => { return s.trim(); });
		parameters.forEach(param => {
			lua.executeAndPrint(param+" = Variable('"+param+"')");
		});

		ids.constants.value
		.split(',')
		.map(s => {
			s = s.trim();
			const eqs = s
			.split('=')
			.map(side => { return side.trim(); });
			const param = eqs[0];
			const value = eqs[1];
			lua.executeAndPrint(param+" = Constant("+value+")");
		});

		let propertiesHTML = '';

		const luaEqnToJSFunc = (eqn, params, texLabel) => {
			let failed = false;
			//here I'm using output for errors
			//since directing error doesn't work -- all errors result in stdout printing "ERROR attempt to call string"
			lua.capture({
				callback : function() {
					lua.execute("eqn = simplify("+eqn+")");
					lua.execute("if type(eqn) == 'number' then eqn = Constant(eqn) end");
				},
				output : function(s) {
					//don't throw -- lua.execute will catch it.
console.log('Lua error!', s);
					failed = true;
				}
			});
			if (failed) throw 'Lua error!';

			//here I'm using output for capturing the compiled lua code
			// I'm recording errors if the captured code fails to compile in JavaScript
			let resultFunction = undefined;
			lua.capture({
				callback : () => {
					//execute it as a single line, so output() could capture it all at once (because output() seems to be called line-by-line)
					let luaCmd = `
print((
	require 'symmath.export.JavaScript'
	:toFuncCode{
		output = {eqn},
		input = {`+params+`},
	}:gsub('\\n', ' ')
))
`;
console.log('in-callback prints', lua.LuaModule.print, lua.LuaModule.printErr);
console.log('executing lua', luaCmd);
					//print commands are going to the old output ...
					lua.execute(luaCmd);
					//TODO if lua has a syntax error, I just get "attempt to call a string value"
					// and I think this is going on inside of lua.vm.js ... time to replace it yet?
				},
				output : jsCmd => {
console.log('got JS output', jsCmd);
					try {
						//seems I used to be able to eval("function(){}") and get a function value back ... not anymore?
						//resultFunction = eval(jsCmd);
						// did eval() just get worse with ES6?
						// now I can't find any way for eval() to return a function, except by assigning it to a global (window)
						//ironic if ES6 is phasing out 'window' access in modules.
						//javascript is such a shitty language ...
						eval('window.__tmpf = '+jsCmd);
						resultFunction = window.__tmpf;
						window.__tmpf = undefined;
					} catch (e) {
console.log("Lua error!", e);
						failed = true;
					}
				},
			});
			if (failed) {
				throw 'Lua error!';
			}

			//while we're here, let's store the LaTex generated from the equations ...
			lua.capture({
				callback : () => {
					let luaCmd = "print((require 'symmath.export.LaTeX'(eqn)))"
console.log('executing lua '+luaCmd);
					lua.execute(luaCmd);
				},
				output : TeX => {
console.log('got TeX output '+TeX);
					if (texLabel !== undefined) TeX = texLabel + ' = {' + TeX + '}';
					propertiesHTML += '$ ' + TeX + ' $<br>\n';
				},
			});

			if (!resultFunction) throw "failed to create symmath->javascript function!";
			return resultFunction;
		};

		//compile equations here
		const coordCallbacks = [];	//coordCallbacks[xyz]
		const coordDerivs = [];		//coordDerivs[xyz][uv|normal]
		const metric = [];			//metric[uv][uv] = g_ij = e_i dot e_j
		const metricEqn = [];
		const invMetricEqn = [];		//invMetricEqn[uv][uv] = g^ij
		const diffMetricEqn = [];		//diffMetricEqn[uv][uv][uv] = g_ij,k
		const conn1stEqn = [];		//conn1stEqn[uv][uv][uv] = conn_ijk = 1/2 (g_ij,k + g_ik,j - g_jk,i)
		const conn2nd = [];			//conn2nd[uv][uv][uv] = conn^i_jk = g^il * conn_ljk
		const conn2ndEqn = [];
console.log('generating coordinate chart functions...');
		const equations = [];
		let failed = false;
		for (let i = 0; i < coordLabels.length; ++i) {
			equations[i] = ids['equation'+coordLabels[i].toUpperCase()].value;
		}
		try {
			//get functions
			equations.forEach((eqn,i) => {
				coordCallbacks[i] = luaEqnToJSFunc(eqn, parameters.join(', '), coordLabels[i]);
				coordDerivs[i] = [];
			});
			equations.forEach((eqn,i) => {
				parameters.forEach((param,j) => {
					coordDerivs[i][j] = luaEqnToJSFunc('diff('+eqn+', '+param+')', parameters.join(', '), '{{d'+coordLabels[i]+'} \\over {d'+param+'}}');
				});
			});

			//get metric and its derivative
			parameters.forEach((u,i) => {
				metric[i] = [];
				metricEqn[i] = [];
				diffMetricEqn[i] = [];
				parameters.forEach((v,j) => {
					let sum = [];
					equations.forEach((eqn,k) => {
						sum.push('diff('+eqn+', '+u+') * diff('+eqn+', '+v+')');
					});
					metricEqn[i][j] = sum.join(' + ');
					metric[i][j] = luaEqnToJSFunc(metricEqn[i][j], parameters.join(', '), 'g_{'+u+' '+v+'}');
					diffMetricEqn[i][j] = [];
					parameters.forEach((w,k) => {
						diffMetricEqn[i][j][k] = 'diff('+metricEqn[i][j]+', '+w+')';
					});
				});
			});

			//inverse metric
			if (parameters.length != 2) throw 'only works with two parameters';
			const metricDetEqn = '('+metricEqn[0][0]+') * ('+metricEqn[1][1]+') - ('+metricEqn[1][0]+') * ('+metricEqn[0][1]+')';
			const invMetricEqn = [
				[	'('+metricEqn[1][1]+') / ('+metricDetEqn+')',	//00
					'-('+metricEqn[0][1]+') / ('+metricDetEqn+')'],	//01
				[	'-('+metricEqn[1][0]+') / ('+metricDetEqn+')',	//10
					'('+metricEqn[0][0]+') / ('+metricDetEqn+')']
			];

			//1st conns
			parameters.forEach((u,i) => {
				conn1stEqn[i] = [];
				parameters.forEach((v,j) => {
					conn1stEqn[i][j] = [];
					parameters.forEach((w,k) => {
						conn1stEqn[i][j][k] = '1/2 * ('+diffMetricEqn[i][j][k]+' + '+diffMetricEqn[i][k][j]+' - '+diffMetricEqn[j][k][i]+')';
					});
				});
			});

			//2nd conns
			parameters.forEach((u,i) => {
				conn2nd[i] = [];
				conn2ndEqn[i] = [];
				parameters.forEach((v,j) => {
					conn2nd[i][j] = [];
					conn2ndEqn[i][j] = [];
					parameters.forEach((w,k) => {
						const sum = [];
						for (let l = 0; l < parameters.length; ++l) {
							sum.push('('+invMetricEqn[i][l]+') * ('+conn1stEqn[l][j][k]+')');
						}
						conn2ndEqn[i][j][k] = sum.join(' + ');
						conn2nd[i][j][k] = luaEqnToJSFunc(conn2ndEqn[i][j][k], parameters.join(', '), '{Gamma^'+u+'}_{'+v+' '+w+'}');
					});
				});
			});

			//get normal
			equations.forEach((eqn,i) => {
				const j = (i+1)%equations.length;
				const k = (j+1)%equations.length;
				const dj_du = coordDerivs[j][0];
				const dk_du = coordDerivs[k][0];
				const dj_dv = coordDerivs[j][1];
				const dk_dv = coordDerivs[k][1];
				coordDerivs[i][2] = function() {
					return dj_du.apply(undefined, arguments)
						* dk_dv.apply(undefined, arguments)
						- dk_du.apply(undefined,arguments)
						* dj_dv.apply(undefined,arguments);
				};
			});
		} catch (e) {
			console.log("failed: ", e);
			failed = true;
		}
		if (failed) {
			//color the failed input textarea red
			allInputs.forEach(input => {
				input.style.background = 'rgb(255,127,127)';
			});
			return;
		} else {
			allInputs.forEach(input => {
				input.style.background = 'white';
			});
		}
		currentCoordChart.mapping = function(coord) {
			const result = [];
			for (let k = 0; k < coordLabels.length; ++k) {
				result[k] = coordCallbacks[k].apply(undefined, coord);
			}
			return result;
		};
		currentCoordChart.diff = function(coord, dim) {
			const result = [];
			for (let k = 0; k < coordLabels.length; ++k) {
				result[k] = coordDerivs[k][dim].apply(undefined, coord);
			}
			return result;
		};
		currentCoordChart.unitDiff = function(coord, dim) {
			const diff = this.diff(coord, dim);
			const len = vec3.length(diff);
			return [diff[0] / len, diff[1] / len, diff[2] / len];
		};
		currentCoordChart.connection = function(coord, k, j) {
			const result = [];
			for (let i = 0; i < parameters.length; ++i) {
				result[i] = conn2nd[i][j][k].apply(undefined, coord);
			}
			return result;
		};

		//process TeX
		const symbols = ['theta', 'phi', 'rho', 'gamma', 'Gamma'];
		symbols.forEach((sym,i) => {
			propertiesHTML = propertiesHTML.replace(new RegExp(sym, 'g'), '\\'+sym);
		});
		ids.properties.innerHTML = propertiesHTML;
		typeset(() => [ids.properties]);

		reset();	//regenerate mesh
	};
	//wait for done lua loading
	let loadingInterval = setInterval(() => {
		if (luaDoneLoading) {
			clearInterval(loadingInterval);
			loadingInterval = undefined;
			doUpdateEquations();
		}
	}, 100);

};

let delayToUpdateEquations = undefined;
{
	let updateInterval = undefined;
	delayToUpdateEquations = () => {
		if (updateInterval !== undefined) {
			clearInterval(updateInterval);
		}
		updateInterval = setTimeout(() => {
			updateEquations();
		}, 1000);
	};
}

allInputs.forEach((input,i) => {
	input.addEventListener('change', delayToUpdateEquations);
	input.addEventListener('keypress', delayToUpdateEquations);
	input.addEventListener('paste', delayToUpdateEquations);
	input.addEventListener('input', delayToUpdateEquations);
});

Object.entries(coordCharts).forEach(entry => {
	const [name,coordChart] = entry;
	const option = DOM('option', {
		text : name,
		value : name,
		appendTo : ids.coordinateSystem,
	});
	if (coordChart[name] == currentCoordChart) option.setAttribute('selected', 'selected');
});
ids.coordinateSystem.addEventListener('change', e => {
	selectCoordChart(ids.coordinateSystem.value);
});
function selectCoordChart(name) {
	currentCoordChart = coordCharts[name];
	coordLabels.forEach((x,i) => { ids['equation'+x.toUpperCase()].value = '0'; });
	if (currentCoordChart.equations !== undefined) {
		currentCoordChart.equations.forEach((equation,i) => {
			ids['equation'+coordLabels[i].toUpperCase()].value = equation;
		});
	}
	ids.parameters.value = currentCoordChart.parameters.join(', ');
	ids.constants.value = '';
	if (currentCoordChart.constants !== undefined) {
		let constantText = [];
		Object.entries(currentCoordChart.constants).forEach(entry => {
			const [k,v] = entry;
			constantText.push(k+'='+v);
		});
		if (constantText.length > 0) {
			ids.constants.value = constantText.concat();
		}
	}
	updateEquations();
}

const tmpQ = quat.create();
mouse = new Mouse3D({
	pressObj : canvas,
	move : function(dx,dy) {
		if (inputState == 'select') {
			let bestCoord = findClickedCoord();

			if (bestCoord !== undefined) {
				selectCoord(bestCoord);
			}
		}

		if (inputState == 'rotate') {
			let rotAngle = Math.PI / 180 * .03 * Math.sqrt(dx*dx + dy*dy);
			quat.setAxisAngle(tmpQ, [dy, dx, 0], rotAngle);

			quat.mul(sceneObj.angle, tmpQ, sceneObj.angle);
			quat.normalize(sceneObj.angle, sceneObj.angle);
		}

		if (inputState == 'direction') {
			//re-integrate geodesic and uncorrected coordinates along surface
			let clickedCoord = findClickedCoord();
			if (clickedCoord !== undefined) {
				let mappedCurrentCoord = currentCoordChart.mapping(currentCoord);
				//TODO
				//rotate view ray and origin by inverse angle
				//intersect view ray with plane of surface normal and mapped coordinate


				let mappedClickedCoord = currentCoordChart.mapping(clickedCoord);
				let delta = [
					mappedClickedCoord[0] - mappedCurrentCoord[0],
					mappedClickedCoord[1] - mappedCurrentCoord[1],
					mappedClickedCoord[2] - mappedCurrentCoord[2]
				];
				let invAngle = quat.create();
				quat.conjugate(invAngle, glutil.view.angle);
				vec3.transformQuat(delta, delta, invAngle);
				let basis0 = currentCoordChart.unitDiff(currentCoord, 0);
				let basis1 = currentCoordChart.unitDiff(currentCoord, 1);
				let d0 = vec3.dot(basis0, delta);
				let d1 = vec3.dot(basis1, delta);
				let dl = vec2.length([d0, d1]);
				chooseDirection([d0/dl, d1/dl]);
			}
		}
	},
	zoom : function(dz) {
		glutil.view.pos[2] += .0001 * dz;
	}
});

const offset = [ [0,0], [1,0], [1,1], [1,1], [0,1], [0,0] ];
let indexes = [];
for (let itheta = 0; itheta < intDivs[0]-1; ++itheta) {
	for (let iphi = 0; iphi < intDivs[1]-1; ++iphi) {
		for (let k = 0; k < offset.length; ++k) {
			indexes.push(
				iphi + offset[k][0] + (itheta + offset[k][1]) * intDivs[1]
			);
		}
	}
}

glutil.view.pos[2] = 2;
glutil.view.zFar = 100;
glutil.view.zNear = .1;
const plainShader = new glutil.Program({
	vertexCode : `
in vec3 vertex;
uniform mat4 projMat;
uniform mat4 mvMat;
void main() {
	gl_Position = projMat * mvMat * vec4(vertex, 1.);
}
`,
	fragmentCode : `
uniform vec4 color;
out vec4 fragColor;
void main() {
	fragColor = color;
}
`
});
const meshShader = new glutil.Program({
	vertexCode : `
in vec3 vertex;
in vec2 intCoord;
in vec3 normal;
uniform mat4 projMat;
uniform mat4 mvMat;
out vec2 intCoordV;
out vec3 normalV;
out vec3 vertexV;
void main() {
	normalV = (mvMat * vec4(normal, 0.)).xyz;
	vec4 mvtx = mvMat * vec4(vertex, 1.);
	vertexV = mvtx.xyz;
	gl_Position = projMat * mvtx;
	intCoordV = intCoord;
}
`,
	fragmentCode : `
uniform vec4 color;
in vec2 intCoordV;
in vec3 normalV;
in vec3 vertexV;
out vec4 fragColor;
void main() {
	vec3 n = normalize(normalV);
	if (n.z < 0.) n = -n;	//backface lighting

	vec2 fc = mod(intCoordV.xy / 10., 1.); //grid
	float i = 1. - 8. * fc.x * fc.y * (1. - fc.x) * (1. - fc.y);
	i = pow(i, 50.);

	fragColor = vec4(.25, .5, .5, 1.);
	fragColor.rgb *= 1. - i;
	vec3 u = normalize(vertexV);
	float l = dot(n, u);
	fragColor.rgb *= max(abs(l), .3);
}
`
});
sceneObj = new glutil.SceneObject({
	static : false
});

meshObj = new glutil.SceneObject({
	parent : sceneObj,
	mode : gl.TRIANGLES,
	attrs : {
		vertex : new glutil.ArrayBuffer({
			count : intDivs[0] * intDivs[1],
			usgae : gl.DYNAMIC_DRAW
		}),
		coord : new glutil.ArrayBuffer({
			dim : 2,
			count : intDivs[0] * intDivs[1]
		}),
		intCoord : new glutil.ArrayBuffer({
			dim : 2,
			count : intDivs[0] * intDivs[1],
			usage : gl.DYNAMIC_DRAW
		}),
		normal : new glutil.ArrayBuffer({
			count : intDivs[0] * intDivs[1],
			usage : gl.DYNAMIC_DRAW
		})
	},
	indexes : new glutil.ElementArrayBuffer({
		data : new Uint16Array(indexes)
	}),
	uniforms : {
		color : [1,1,1,1]
	},
	shader : meshShader,
	static : false
});

selectionObj = new glutil.SceneObject({
	parent : sceneObj,
	mode : gl.LINE_LOOP,
	attrs : {
		vertex : new glutil.ArrayBuffer({
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

let basisColors = [
	[1,0,0,1],
	[0,1,0,1]
];
let connColors = [
	[
		[1,1,0,1],
		[0,1,1,1]
	],
	[
		[1,0,1,1],
		[1,.5,0,1]
	]
];
for (let i = 0; i < 2; ++i) {
	basisObjs[i] = new glutil.SceneObject({
		parent : sceneObj,
		mode : gl.LINES,
		attrs : {
			vertex : new glutil.ArrayBuffer({
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
	for (let j = 0; j < 2; ++j) {
		connObjs[i][j] = new glutil.SceneObject({
			parent : sceneObj,
			mode : gl.LINES,
			attrs : {
				vertex : new glutil.ArrayBuffer({
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

partialPathObj = new glutil.SceneObject({
	parent : sceneObj,
	mode : gl.LINE_STRIP,
	attrs : {
		vertex : new glutil.ArrayBuffer({
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
geodesicPathObj = new glutil.SceneObject({
	parent : sceneObj,
	mode : gl.LINE_STRIP,
	attrs : {
		vertex : new glutil.ArrayBuffer({
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

window.addEventListener('resize', resize);
resize();

update();
