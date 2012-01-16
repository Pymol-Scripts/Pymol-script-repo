/*
    GLmol - Molecular Viewer on WebGL/Javascript  BETA version (0.37)

   (C) Copyright 2011, biochem_fan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    This program uses
      Three.js 
         https://github.com/mrdoob/three.js
         Copyright (c) 2010-2011 three.js Authors. All rights reserved.
      jQuery
         http://jquery.org/
         Copyright (c) 2011 John Resig
 */

// Patch MeshLambertMaterial's shader to fix lighting in doubleSided mesh.
// This trick is not applicable to MeshPhongMaterial
THREE.ShaderLib.lambert.vertexShader = THREE.ShaderLib.lambert.vertexShader.replace("vec3 transformedNormal = normalize( normalMatrix * normal )", "vec3 transformedNormal = normalize( normalMatrix * normal );\nif (transformedNormal.z < 0.0) transformedNormal *= -1.0;");

THREE.Geometry.prototype.colorAll = function (color) {
   for (var i = 0; i < this.faces.length; i++) {
      this.faces[i].color = color;
   }
};

THREE.Matrix4.prototype.isIdentity = function() {
   for (var i = 1; i <= 4; i++) {
      for (var j = 1; j <= 4; j++) {
         var shouldBe = (i == j) ? 1 : 0;
         if (this['n' + i + j] != shouldBe) return false;
      }
   }
   return true;
};

var GLmol = (function() {
function GLmol(id, suppressAutoload) {
   this.Nucleotides = ['  G', '  A', '  T', '  C', '  U', ' DG', ' DA', ' DT', ' DC', ' DU'];
   this.ElementColors = {"H": 0xCCCCCC, "C": 0xAAAAAA, "O": 0xCC0000, "N": 0x0000CC, "S": 0xCCCC00, "P": 0x6622CC,
                         "F": 0x00CC00, "CL": 0x00CC00, "BR": 0x882200, "I": 0x6600AA,
                         "FE": 0xCC6600, "CA": 0x8888AA};
// Reference: A. Bondi, J. Phys. Chem., 1964, 68, 441.
   this.vdwRadii = {"H": 1.2, "Li": 1.82, "Na": 2.27, "K": 2.75, "C": 1.7, "N": 1.55, "O": 1.52,
                   "F": 1.47, "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "SE": 1.90,
                   "ZN": 1.39, "CU": 1.4, "NI": 1.63};

   this.id = id;
   this.aaScale = 1; // or 2

   this.container = $('#' + this.id);
   this.WIDTH = this.container.width() * this.aaScale, this.HEIGHT = this.container.height() * this.aaScale;
   this.ASPECT = this.WIDTH / this.HEIGHT;
   this.NEAR = 1, FAR = 800;
   this.CAMERA_Z = -150;
   this.renderer = new THREE.WebGLRenderer({antialias: true});
   // 'antialias: true' works only in Chrome. Firefox 10 will support this.
   // setting this.aaScale = 2 will enable antialias in Firefox but GPU load increases.
   this.renderer.domElement.style.width = "100%";
   this.renderer.domElement.style.height = "100%";
   this.container.append(this.renderer.domElement);
   this.renderer.setSize(this.WIDTH, this.HEIGHT);

   this.camera = new THREE.PerspectiveCamera(20, this.ASPECT, 1, 800); // will be updated anyway
   this.camera.position = new THREE.Vector3(0, 0, this.CAMERA_Z);
   this.camera.lookAt(new THREE.Vector3(0, 0, 0));
   this.perspectiveCamera = this.camera;
   this.orthoscopicCamera = new THREE.OrthographicCamera();
   this.orthoscopicCamera.position.z = this.CAMERA_Z;
   this.orthoscopicCamera.lookAt(new THREE.Vector3(0, 0, 0));

   var self = this;
   $(window).resize(function() { // only window can capture resize event
      self.WIDTH = self.container.width() * self.aaScale;
      self.HEIGHT = self.container.height() * self.aaScale;
      self.ASPECT = self.WIDTH / self.HEIGHT;
      self.renderer.setSize(self.WIDTH, self.HEIGHT);
      self.camera.aspect = self.ASPECT;
      self.camera.updateProjectionMatrix();
      self.show();
   });

   this.scene = null;
   this.rotationGroup = null; // which contains modelGroup
   this.modelGroup = null;

   this.bgColor = 0x000000;
   this.fov = 20;
   this.fogStart = 0.4;
   this.slabNear = -50; // relative to the center of rotationGroup
   this.slabFar = +50;

   // Default values
   this.sphereRadius = 1.5; 
   this.cylinderRadius = 0.2;
   this.lineWidth = 1.5 * this.aaScale;
   this.curveWidth = 3 * this.aaScale;
   this.defaultColor = 0xCCCCCC;
   this.sphereQuality = 16;
   this.cylinderQuality = 8;
   this.axisDIV = 5; // 3 still gives acceptable quality
   this.strandDIV = 6;
   this.nucleicAcidStrandDIV = 4;
   this.tubeDIV = 8;
   this.coilWidth = 0.3;
   this.helixSheetWidth = 1.3;
   this.nucleicAcidWidth = 0.8;
 
   // UI variables
   this.cq = new THREE.Quaternion(1, 0, 0, 0);
   this.dq = new THREE.Quaternion(1, 0, 0, 0);
   this.isDragging = false;
   this.mouseStartX = 0;
   this.mouseStartY = 0;
   this.currentModelPos = 0;
   this.cz = 0;
   this.enableMouse();

   if (suppressAutoload) return;
   this.loadMolecule();
}

GLmol.prototype.setupLights = function(scene) {
   var directionalLight =  new THREE.DirectionalLight(0xFFFFFF);
   directionalLight.position = new THREE.Vector3(30, 30, -150);
   directionalLight.intensity = 0.6;
   scene.add(directionalLight);
   var ambientLight = new THREE.AmbientLight(0x808080);
   scene.add(ambientLight);

   scene.add(directionalLight);  
};

GLmol.prototype.parseSDF = function(str) {
   var atoms = this.atoms;
   var protein = this.protein;

   var lines = str.split("\n");
   if (lines.length < 4) return;
   var atomCount = parseInt(lines[3].substr(0, 3));
   if (isNaN(atomCount) || atomCount <= 0) return; // might be a PDB file.
   var bondCount = parseInt(lines[3].substr(3, 3));
   var offset = 4;
   if (lines.length < 4 + atomCount + bondCount) return;
   for (var i = 1; i <= atomCount; i++) {
      var line = lines[offset];
      offset++;
      var atom = {};
      atom.serial = i;
      atom.x = parseFloat(line.substr(0, 10));
      atom.y = parseFloat(line.substr(10, 10));
      atom.z = parseFloat(line.substr(20, 10));
      atom.hetflag = true;
      atom.atom = atom.elem = line.substr(30, 3).replace(/ /g, "");
      atom.bonds = [];
      atom.bondOrder = [];
      atoms[i] = atom;
   }
   for (i = 1; i <= bondCount; i++) {
      var line = lines[offset];
      offset++;
      var from = parseInt(line.substr(0, 3));
      var to = parseInt(line.substr(3, 3));
      var order = parseInt(line.substr(6, 3));
      atoms[from].bonds.push(to);
      atoms[from].bondOrder.push(order);
      atoms[to].bonds.push(from);
      atoms[to].bondOrder.push(order);
   }

   protein.sdf = true;
};

GLmol.prototype.parsePDB2 = function(str) {
   var atoms = this.atoms;
   var protein = this.protein;
   var molID;

   var atoms_cnt = 0;
   lines = str.split("\n");
   for (var i = 0; i < lines.length; i++) {
      line = lines[i].replace(/^\s*/, ''); // remove indent
      var recordName = line.substr(0, 6);
      if (recordName == 'ATOM  ' || recordName == 'HETATM') {
         var atom, resn, chain, resi, x, y, z, hetflag, elem, serial, altLoc, b
         serial = parseInt(line.substr(6, 5));
         atom = line.substr(12, 4).replace(/ /g, "");
         altLoc = line.substr(16, 1);
         if (altLoc != ' ' && altLoc != 'A') continue; // FIXME: ad hoc
         resn = line.substr(17, 3);
         chain = line.substr(21, 1);
         resi = parseInt(line.substr(22, 5)); 
         x = parseFloat(line.substr(30, 8));
         y = parseFloat(line.substr(38, 8));
         z = parseFloat(line.substr(46, 8));
         b = parseFloat(line.substr(60, 8));
         elem = line.substr(76, 2).replace(/ /g, "");
         if (elem == '') { // for some incorrect PDB files
            elem = line.substr(12, 4).replace(/ /g,"");
         }
         if (line[0] == 'H') hetflag = true;
         else hetflag = false;
         if (!protein[chain]) protein[chain] = {'residues':[]};
         atoms[serial] = {'resn': resn, 'x': x, 'y': y, 'z': z, 'elem': elem,
  'hetflag': hetflag, 'chain': chain, 'resi': resi, 'serial': serial, 'atom': atom,
  'bonds': [], 'ss': 'c', 'color': 0xFFFFFF, 'bonds': [], 'bondOrder': [], 'b': b /*', altLoc': altLoc*/};
      } else if (recordName == 'SHEET ') {
         var startChain = line.substr(21, 1);
         var startResi = parseInt(line.substr(22, 4));
         var endChain = line.substr(32, 1);
         var endResi = parseInt(line.substr(33, 4));
         protein.sheet.push([startChain, startResi, endChain, endResi]);
     } else if (recordName == 'CONECT') {
// MEMO: We don't have to parse SSBOND, LINK because both are also 
// described in CONECT. But what about 2JYT???
         var from = parseInt(line.substr(6, 5));
         for (var j = 0; j < 4; j++) {
            var to = parseInt(line.substr([11, 16, 21, 26][j], 5));
            if (isNaN(to)) continue;
            if (atoms[from] != undefined) {
               atoms[from].bonds.push(to);
               atoms[from].bondOrder.push(1);
            }
         }
     } else if (recordName == 'HELIX ') {
         var startChain = line.substr(19, 1);
         var startResi = parseInt(line.substr(21, 4));
         var endChain = line.substr(31, 1);
         var endResi = parseInt(line.substr(33, 4));
         protein.helix.push([startChain, startResi, endChain, endResi]);
     } else if (recordName == 'CRYST1') {
         protein.a = parseFloat(line.substr(6, 9));
         protein.b = parseFloat(line.substr(15, 9));
         protein.c = parseFloat(line.substr(24, 9));
         protein.alpha = parseFloat(line.substr(33, 7));
         protein.beta = parseFloat(line.substr(40, 7));
         protein.gamma = parseFloat(line.substr(47, 7));
         protein.spacegroup = line.substr(55, 11);
         this.defineCell();
      }  else if (recordName.substr(0, 5) == 'SCALE') {
         var n = parseInt(recordName.substr(5, 1));
         protein.scaleMatrix['n' + n + '1'] = parseFloat(line.substr(10, 10));
         protein.scaleMatrix['n' + n + '2'] = parseFloat(line.substr(20, 10));
         protein.scaleMatrix['n' + n + '3'] = parseFloat(line.substr(30, 10));
         protein.scaleMatrix['n' + n + '4'] = parseFloat(line.substr(45, 10));
      } else if (recordName == 'REMARK') {
         if (line.substr(13, 5) == 'BIOMT') {
            var n = parseInt(line[18]);
            var m = parseInt(line.substr(21, 2));
            if (protein.biomtMatrices[m] == undefined) protein.biomtMatrices[m] = new THREE.Matrix4().identity();
            protein.biomtMatrices[m]['n' + n + '1'] = parseFloat(line.substr(24, 9));
            protein.biomtMatrices[m]['n' + n + '2'] = parseFloat(line.substr(34, 9));
            protein.biomtMatrices[m]['n' + n + '3'] = parseFloat(line.substr(44, 9));
            protein.biomtMatrices[m]['n' + n + '4'] = parseFloat(line.substr(54, 10));

         } else if (line.substr(13, 5) == 'SMTRY') {
            var n = parseInt(line[18]);
            var m = parseInt(line.substr(21, 2));
            if (protein.symmetryMatrices[m] == undefined) protein.symmetryMatrices[m] = new THREE.Matrix4().identity();
            protein.symmetryMatrices[m]['n' + n + '1'] = parseFloat(line.substr(24, 9));
            protein.symmetryMatrices[m]['n' + n + '2'] = parseFloat(line.substr(34, 9));
            protein.symmetryMatrices[m]['n' + n + '3'] = parseFloat(line.substr(44, 9));
            protein.symmetryMatrices[m]['n' + n + '4'] = parseFloat(line.substr(54, 10));
         }
      } else if (recordName == 'HEADER') {
         protein.pdbID = line.substr(62, 4);
      } else if (recordName == 'TITLE ') {
         if (protein.title == undefined) protein.title = "";
            protein.title += line.substr(10, 70) + "\n"; // CHECK: why 60 is not enough???
      } else if (recordName == 'COMPND') {
              // TODO: Implement me!
      }
   }

   // Assign secondary structures 
   for (i = 0; i < atoms.length; i++) {
      atom = atoms[i]; if (atom == undefined) continue;

      var found = false;
      // MEMO: Can start chain and end chain differ?
      for (j = 0; j < protein.sheet.length; j++) {
         if (atom.chain != protein.sheet[j][0]) continue;
         if (atom.resi < protein.sheet[j][1]) continue;
         if (atom.resi > protein.sheet[j][3]) continue;
         atom.ss = 's';
         if (atom.resi == protein.sheet[j][1]) atom.ssbegin = true;
         if (atom.resi == protein.sheet[j][3]) atom.ssend = true;
      }
      for (j = 0; j < protein.helix.length; j++) {
         if (atom.chain != protein.helix[j][0]) continue;
         if (atom.resi < protein.helix[j][1]) continue;
         if (atom.resi > protein.helix[j][3]) continue;
         atom.ss = 'h';
         if (atom.resi == protein.helix[j][1]) atom.ssbegin = true;
         else if (atom.resi == protein.helix[j][3]) atom.ssend = true;
      }
   } 
};

// Catmull-Rom subdivision
GLmol.prototype.subdivide = function(_points, DIV) { // points as Vector3
   var ret = [];
   var points = _points;
   /*points = new Array(); // Smoothing test for beta sheet.
   points.push(_points[0]);
   for (var i = 0, lim = _points.length - 1; i < lim; i++) {
      var p1 = _points[i], p2 = _points[i + 1];
      points.push(new THREE.Vector3((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2));
   }
   points.push(_points[_points.length - 1]);*/

   for (var i = -1, size = points.length; i <= size - 3; i++) {
      var p0 = points[(i == -1) ? 0 : i];
      var p1 = points[i + 1], p2 = points[i + 2];
      var  p3 = points[(i == size - 3) ? size - 1 : i + 3];
      var v0 = new THREE.Vector3().sub(p2, p0).multiplyScalar(0.5);
      var v1 = new THREE.Vector3().sub(p3, p1).multiplyScalar(0.5);
      for (var j = 0; j < DIV; j++) {
         var t = 1.0 / DIV * j;
         var x = p1.x + t * v0.x 
                  + t * t * (-3 * p1.x + 3 * p2.x - 2 * v0.x - v1.x)
                  + t * t * t * (2 * p1.x - 2 * p2.x + v0.x + v1.x);
         var y = p1.y + t * v0.y 
                  + t * t * (-3 * p1.y + 3 * p2.y - 2 * v0.y - v1.y)
                  + t * t * t * (2 * p1.y - 2 * p2.y + v0.y + v1.y);
         var z = p1.z + t * v0.z 
                  + t * t * (-3 * p1.z + 3 * p2.z - 2 * v0.z - v1.z)
                  + t * t * t * (2 * p1.z - 2 * p2.z + v0.z + v1.z);
         ret.push(new THREE.Vector3(x, y, z));
      }
   }
   ret.push(points[points.length - 1]);
   return ret;
};

GLmol.prototype.drawAtomsAsSpheres = function(group, atomlist, atomR, defaultRadius, solventRadius) {
   var points = [];
   var colors = [];
   var radii = [];

   for (var i = 0; i < atomlist.length; i++) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      var r = (this.vdwRadii[atom.elem] != undefined) ? this.vdwRadii[atom.elem] : defaultRadius;
      r += solventRadius;
      points.push([atom.x, atom.y, atom.z]); colors.push(atom.color); radii.push(r);
   }
   this.drawSpheres(group, points, colors, radii);
};

// 496msec for 2580atoms, 2580msec for 12641atoms (draw 7msec); div = 1, mesh
GLmol.prototype.drawSpheres = function(group, points, colors, radii) {
   var time = new Date();
   var pgeo = new THREE.Geometry();
   var ret = pgeo.vertices;
   var ret2 = pgeo.colors;
   var ret3 = pgeo.faces;

   var geo = new THREE.IcosahedronGeometry(1);
   var basePoints = geo.vertices, baseFaces = geo.faces;
   var jlim = basePoints.length;
   var klim = geo.faces.length;
   var offset = 0;
   for (var i = 0, lim = points.length; i < lim; i++) {
      var x = points[i][0], y = points[i][1], z = points[i][2], r = radii[i], c = new THREE.Color(colors[i]);
      for (var j = 0; j < jlim; j++) {
              ret.push(new THREE.Vertex(new THREE.Vector3(x + basePoints[j].position.x * r,
                       y + basePoints[j].position.y * r, z + basePoints[j].position.z * r)));
//              ret2.push(c);
      }
      for (var k = 0; k < klim; k++) {
         var tmp = baseFaces[k];
         var f = new THREE.Face3(tmp.a + offset, tmp.b + offset, tmp.c + offset);
         f.color = c; f.vertexNormals = tmp.vertexNormals;
         ret3.push(f);
      }
      offset += jlim;
   }
//   var ps = new THREE.ParticleSystem(pgeo, new THREE.ParticleBasicMaterial({vertexColors: true, size: 0.2}));
   var ps = new THREE.Mesh(pgeo, new THREE.MeshLambertMaterial({vertexColors: true}));
   console.log("drawDots: " + (+new Date() - time) + "ms for " + ret.length + "dots");
   group.add(ps);
};


GLmol.prototype.drawAtomsAsSphere = function(group, atomlist, defaultRadius, forceDefault) {
   var sphereGeometry = new THREE.SphereGeometry(1, this.sphereQuality, this.sphereQuality); // r, seg, ring
   for (var i = 0; i < atomlist.length; i++) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      var sphereMaterial = new THREE.MeshLambertMaterial({color: atom.color});
      var sphere = new THREE.Mesh(sphereGeometry, sphereMaterial);
      group.add(sphere);
      var r = (!forceDefault && this.vdwRadii[atom.elem] != undefined) ? this.vdwRadii[atom.elem] : defaultRadius
      sphere.scale.x = sphere.scale.y = sphere.scale.z = r;
      sphere.position.x = atom.x;
      sphere.position.y = atom.y;
      sphere.position.z = atom.z;
   }
};

// about two times faster than sphere when div = 2
GLmol.prototype.drawAtomsAsIcosahedron = function(group, atomlist, defaultRadius, forceDefault) {
   var geo = this.IcosahedronGeometry();
   for (var i = 0; i < atomlist.length; i++) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      var mat = new THREE.MeshLambertMaterial({color: atom.color});
      var sphere = new THREE.Mesh(geo, mat);
      sphere.scale.x = sphere.scale.y = sphere.scale.z = (!forceDefault && this.vdwRadii[atom.elem] != undefined) ? this.vdwRadii[atom.elem] : defaultRadius;
      group.add(sphere);
      sphere.position.x = atom.x;
      sphere.position.y = atom.y;
      sphere.position.z = atom.z;
   }
};

GLmol.prototype.isConnected = function(atom1, atom2) {
   var s = atom1.bonds.indexOf(atom2.serial);
   if (s != -1) return atom1.bondOrder[s];

   if (this.protein.sdf && (atom1.hetflag || atom2.hetflag)) return 0; // CHECK: or should I ?

   var distSquared = (atom1.x - atom2.x) * (atom1.x - atom2.x) + 
                     (atom1.y - atom2.y) * (atom1.y - atom2.y) + 
                     (atom1.z - atom2.z) * (atom1.z - atom2.z);

//   if (atom1.altLoc != atom2.altLoc) return false;
   if (isNaN(distSquared)) return 0;
   if (distSquared < 0.5) return 0; // maybe duplicate position. FIXME: Is this treatment correct?

   if (distSquared > 1.3 && (atom1.elem == 'H' || atom2.elem == 'H' || atom1.elem == 'D' || atom2.elem == 'D')) return 0;
   if (distSquared < 3.42 && (atom1.elem == 'S' || atom2.elem == 'S')) return 1;
   if (distSquared > 2.78) return 0;
   return 1;
};

GLmol.prototype.drawBondsAsStick = function(group, atomlist, bondR, atomR, ignoreNonbonded) {
   var sphereGeometry = this.IcosahedronGeometry();

   for (var _i in atomlist) {
      var i = atomlist[_i];
      var atom1 = this.atoms[i];
      if (atom1 == undefined) continue;
      var connected = false;

      for (var _j in atomlist) {
         var j = atomlist[_j];
         if (i == j) continue;
         var atom2 = this.atoms[j];
         if (atom1 == undefined || atom2 == undefined) continue;
         if (this.isConnected(atom1, atom2) == 0) continue;
         connected = true;

         this.drawCylinder(group, new THREE.Vector3(atom1.x, atom1.y, atom1.z),
                           new THREE.Vector3((atom1.x + atom2.x) / 2, (atom1.y + atom2.y) / 2, (atom1.z + atom2.z) / 2), bondR, atom1.color);
       }
       if (ignoreNonbonded && !connected) continue;
       var sphereMaterial = new THREE.MeshLambertMaterial({color: atom1.color});
       var sphere = new THREE.Mesh(sphereGeometry, sphereMaterial);
       sphere.scale.x = sphere.scale.y = sphere.scale.z = atomR * 0.95;
       group.add(sphere);
       sphere.position.x = atom1.x;
       sphere.position.y = atom1.y;
       sphere.position.z = atom1.z;
    }
};

GLmol.prototype.defineCell = function() {
    var protein = this.protein;
    if (protein.a == undefined) return;

    protein.ax = protein.a;
    protein.ay = 0;
    protein.az = 0;
    protein.bx = protein.b * Math.cos(Math.PI / 180.0 * protein.gamma);
    protein.by = protein.b * Math.sin(Math.PI / 180.0 * protein.gamma);
    protein.bz = 0;
    protein.cx = protein.c * Math.cos(Math.PI / 180.0 * protein.beta);
    protein.cy = protein.c * (Math.cos(Math.PI / 180.0 * protein.alpha) - 
               Math.cos(Math.PI / 180.0 * protein.gamma) 
             * Math.cos(Math.PI / 180.0 * protein.beta)
             / Math.sin(Math.PI / 180.0 * protein.gamma));
    protein.cz = Math.sqrt(protein.c * protein.c * Math.sin(Math.PI / 180.0 * protein.beta)
               * Math.sin(Math.PI / 180.0 * protein.beta) - protein.cy * protein.cy);
};

GLmol.prototype.drawUnitcell = function(group) {
    var protein = this.protein;
    if (protein.a == undefined) return;

    var vertices = [[0, 0, 0], [protein.ax, protein.ay, protein.az], [protein.bx, protein.by, protein.bz], [protein.ax + protein.bx, protein.ay + protein.by, protein.az + protein.bz],
          [protein.cx, protein.cy, protein.cz], [protein.cx + protein.ax, protein.cy + protein.ay,  protein.cz + protein.az], [protein.cx + protein.bx, protein.cy + protein.by, protein.cz + protein.bz], [protein.cx + protein.ax + protein.bx, protein.cy + protein.ay + protein.by, protein.cz + protein.az + protein.bz]];
    var edges = [0, 1, 0, 2, 1, 3, 2, 3, 4, 5, 4, 6, 5, 7, 6, 7, 0, 4, 1, 5, 2, 6, 3, 7];    

    var geo = new THREE.Geometry();
    for (var i = 0; i < edges.length; i++) {
       geo.vertices.push(new THREE.Vertex(new THREE.Vector3(vertices[edges[i]][0], vertices[edges[i]][1], vertices[edges[i]][2])));
    }
   var lineMaterial = new THREE.LineBasicMaterial({linewidth: 1, color: 0xcccccc});
   var line = new THREE.Line(geo, lineMaterial);
   line.type = THREE.Lines;
   group.add(line);
};


// TODO: Refactor!
GLmol.prototype.drawBondsAsLineSub = function(geo, atom1, atom2, order) {
   var midpoint = new THREE.Vector3((atom1.x + atom2.x) / 2,
                   (atom1.y + atom2.y) / 2, (atom1.z + atom2.z) / 2);
   var delta, dot;
   if (order > 1) { // Find the bond plane. TODO: Find inner side of a ring
      var axis = new THREE.Vector3(atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z);
      var found = null;
      for (var i = 0; i < atom1.bonds.length && !found; i++) {
         var atom = this.atoms[atom1.bonds[i]]; if (!atom) continue;
         if (atom.serial != atom2.serial && atom.elem != 'H') found = atom;
      }
      for (var i = 0; i < atom2.bonds.length && !found; i++) {
         var atom = this.atoms[atom2.bonds[i]]; if (!atom) continue;
         if (atom.serial != atom1.serial && atom.elem != 'H') found = atom;
      }
      if (found) {
         var tmp = new THREE.Vector3(atom1.x - found.x, atom1.y - found.y, atom1.z - found.z);
         dot = tmp.dot(axis);
         delta = new THREE.Vector3(tmp.x - axis.x * dot, tmp.y - axis.y * dot, tmp.z - axis.z * dot);
      }
      if (!found || Math.abs(dot - 1) < 0.001) {
         if (axis.x < 0.01 && axis.y < 0.01) {
            delta = new THREE.Vector3(0, - axis.z, axis.y);
         } else {
            delta = new THREE.Vector3(- axis.y, axis.x, 0);
         }
      }
      delta.normalize().multiplyScalar(0.15);
   }

   var color = new THREE.Color(atom1.color);
   geo.vertices.push(new THREE.Vertex(new THREE.Vector3(atom1.x, atom1.y, atom1.z)));
   geo.colors.push(color);
   geo.vertices.push(new THREE.Vertex(midpoint));
   geo.colors.push(color);
   if (order > 1) {
      geo.vertices.push(new THREE.Vertex(new THREE.Vector3(atom1.x + delta.x, atom1.y + delta.y, atom1.z + delta.z)));
      geo.colors.push(color);
      geo.vertices.push(new THREE.Vertex(new THREE.Vector3(midpoint.x + delta.x, midpoint.y + delta.y, midpoint.z + delta.z)));
      geo.colors.push(color);
   }
   if (order > 2) {
      geo.vertices.push(new THREE.Vertex(new THREE.Vector3(atom1.x + delta.x * 2, atom1.y + delta.y * 2, atom1.z + delta.z * 2)));
      geo.colors.push(color);
      geo.vertices.push(new THREE.Vertex(new THREE.Vector3(midpoint.x + delta.x * 2, midpoint.y + delta.y * 2, midpoint.z + delta.z * 2)));
      geo.colors.push(color);
   }
       
   color = new THREE.Color(atom2.color);
   geo.vertices.push(new THREE.Vertex(new THREE.Vector3(atom2.x, atom2.y, atom2.z)));
   geo.colors.push(color);
   geo.vertices.push(new THREE.Vertex(midpoint));
   geo.colors.push(color);
   if (order > 1) {
      geo.vertices.push(new THREE.Vertex(new THREE.Vector3(atom2.x + delta.x, atom2.y + delta.y, atom2.z + delta.z)));
      geo.colors.push(color);
      geo.vertices.push(new THREE.Vertex(new THREE.Vector3(midpoint.x + delta.x, midpoint.y + delta.y, midpoint.z + delta.z)));
      geo.colors.push(color);
   }
   if (order > 2) {
      geo.vertices.push(new THREE.Vertex(new THREE.Vector3(atom2.x + delta.x * 2, atom2.y + delta.y * 2, atom2.z + delta.z * 3)));
      geo.colors.push(color);
      geo.vertices.push(new THREE.Vertex(new THREE.Vector3(midpoint.x + delta.x * 2, midpoint.y + delta.y * 2, midpoint.z + delta.z * 2)));
      geo.colors.push(color);
   }
};

GLmol.prototype.drawBondsAsLine = function(group, atomlist, lineWidth) {
   var geo = new THREE.Geometry();   
   var nAtoms = atomlist.length;

   for (var _i = 0; _i < nAtoms; _i++) {
      var i = atomlist[_i];
      var  atom1 = this.atoms[i];
      if (atom1 == undefined) continue;
      for (var _j = _i + 1; _j < _i + 30 && _j < nAtoms; _j++) {
         var j = atomlist[_j];
         var atom2 = this.atoms[j];
         if (atom2 == undefined) continue;
         var order = this.isConnected(atom1, atom2);
         if (order == 0) continue;

         this.drawBondsAsLineSub(geo, atom1, atom2, order);
      }
      for (var _j = 0; _j < atom1.bonds.length; _j++) {
          var j = atom1.bonds[_j];
          if (j < i + 30) continue; // be conservative!
          if (atomlist.indexOf(j) == -1) continue;
          var atom2 = this.atoms[j];
          if (atom2 == undefined) continue;
          this.drawBondsAsLineSub(geo, atom1, atom2, atom1.bondOrder[_j]);
      }
    }
   var lineMaterial = new THREE.LineBasicMaterial({linewidth: lineWidth});
   lineMaterial.vertexColors = true;

   var line = new THREE.Line(geo, lineMaterial);
   line.type = THREE.Lines;
   group.add(line);
};

GLmol.prototype.drawSmoothCurve = function(group, _points, width, colors, div) {
   if (_points.length == 0) return;

   div = (div == undefined) ? 5 : div;

   var geo = new THREE.Geometry();
   var points = this.subdivide(_points, div);

   for (var i = 0; i < points.length; i++) {
      geo.vertices.push(new THREE.Vertex(points[i]));
       geo.colors.push(new THREE.Color(colors[(i == 0) ? 0 : Math.round((i - 1) / div)]));
  }
  var lineMaterial = new THREE.LineBasicMaterial({linewidth: width});
  lineMaterial.vertexColors = true;
  var line = new THREE.Line(geo, lineMaterial);
  line.type = THREE.LineStrip;
  group.add(line);
};

GLmol.prototype.drawAsCross = function(group, atomlist, delta) {
   var geo = new THREE.Geometry();
   var points = [[delta, 0, 0], [-delta, 0, 0], [0, delta, 0], [0, -delta, 0], [0, 0, delta], [0, 0, -delta]];
 
   for (var i = 0, lim = atomlist.length; i < lim; i++) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      var c = new THREE.Color(atom.color);
      for (var j = 0; j < 6; j++) {
         geo.vertices.push(new THREE.Vertex(new THREE.Vector3(atom.x + points[j][0], atom.y + points[j][1], atom.z + points[j][2])));
         geo.colors.push(c);
      }
  }
  var lineMaterial = new THREE.LineBasicMaterial({linewidth: this.lineWidth});
  lineMaterial.vertexColors = true;
  var line = new THREE.Line(geo, lineMaterial);
  line.type = THREE.Lines;
  group.add(line);
};

// FIXME: Winkled...
GLmol.prototype.drawSmoothTube = function(group, _points, colors, radii) {
   if (_points.length < 2) return;

   var circleDiv = this.tubeDIV, axisDiv = this.axisDIV;
   var geo = new THREE.Geometry();
   var points = this.subdivide(_points, axisDiv);
   var prevAxis1 = new THREE.Vector3(), prevAxis2;

   for (var i = 0, lim = points.length; i < lim; i++) {
      var r, idx = (i - 1) / axisDiv;
      if (i == 0) r = radii[0];
      else { 
         if (idx % 1 == 0) r = radii[idx];
         else {
            var floored = Math.floor(idx);
            var tmp = idx - floored;
            r = radii[floored] * tmp + radii[floored + 1] * (1 - tmp);
         }
      }
      var delta, axis1, axis2;

      if (i < lim - 1) {
         delta = new THREE.Vector3().sub(points[i], points[i + 1]);
         axis1 = new THREE.Vector3(0, - delta.z, delta.y).normalize().multiplyScalar(r);
         axis2 = new THREE.Vector3().cross(delta, axis1).normalize().multiplyScalar(r);
//      var dir = 1, offset = 0;
         if (prevAxis1.dot(axis1) < 0) {
                 axis1.negate(); axis2.negate();  //dir = -1;//offset = 2 * Math.PI / axisDiv;
         }
         prevAxis1 = axis1; prevAxis2 = axis2;
      } else {
         axis1 = prevAxis1; axis2 = prevAxis2;
      }

      for (var j = 0; j < circleDiv; j++) {
         var angle = 2 * Math.PI / circleDiv * j; //* dir  + offset;
         var c = Math.cos(angle), s = Math.sin(angle);
         geo.vertices.push(new THREE.Vertex(new THREE.Vector3(
         points[i].x + c * axis1.x + s * axis2.x,
         points[i].y + c * axis1.y + s * axis2.y, 
         points[i].z + c * axis1.z + s * axis2.z)));
      }
   }

   var offset = 0;
   for (var i = 0, lim = points.length - 1; i < lim; i++) {
      var c =  new THREE.Color(colors[Math.round((i - 1)/ axisDiv)]);

      if (false) { // FIXME
         var r, idx = (i - 1) / axisDiv;
         if (idx % 1 != 0 && i != 0) {
            var floored = Math.floor(idx);
            var tmp = idx - floored;
            var nextColor = new THREE.Color(colors[floored + 1]);
            c.r = c.r * tmp + nextColor.r * (1 - tmp);
            c.g = c.g * tmp + nextColor.g * (1 - tmp);
            c.b = c.b * tmp + nextColor.b * (1 - tmp);
         }
      }

      var reg = 0;
      var r1 = new THREE.Vector3().sub(geo.vertices[offset].position, geo.vertices[offset + circleDiv].position).lengthSq();
      var r2 = new THREE.Vector3().sub(geo.vertices[offset].position, geo.vertices[offset + circleDiv + 1].position).lengthSq();
      if (r1 > r2) {r1 = r2; reg = 1;};
      for (var j = 0; j < circleDiv; j++) {
          geo.faces.push(new THREE.Face3(offset + j, offset + (j + reg) % circleDiv + circleDiv, offset + (j + 1) % circleDiv));
          geo.faces.push(new THREE.Face3(offset + (j + 1) % circleDiv, offset + (j + reg) % circleDiv + circleDiv, offset + (j + reg + 1) % circleDiv + circleDiv));
          geo.faces[geo.faces.length -2].color = c;
          geo.faces[geo.faces.length -1].color = c;
      }
      offset += circleDiv;
   }
   geo.computeFaceNormals();
   geo.computeVertexNormals(false);
   var mat = new THREE.MeshLambertMaterial();// mat.wireframe = true;
   mat.vertexColors = THREE.FaceColors;
   var mesh = new THREE.Mesh(geo, mat);
   mesh.doubleSided = true;
   group.add(mesh);
};


GLmol.prototype.drawMainchainCurve = function(group, atomlist, curveWidth, atomName, div) {
   var points = [], colors = [];
   var currentChain, currentResi;
   if (div == undefined) div = 5;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == atomName) && !atom.hetflag) {
         if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
            this.drawSmoothCurve(group, points, curveWidth, colors, div);
            points = [];
            colors = [];
         }
         points.push(new THREE.Vector3(atom.x, atom.y, atom.z));
         colors.push(atom.color);
         currentChain = atom.chain;
         currentResi = atom.resi;
      }
   }
    this.drawSmoothCurve(group, points, curveWidth, colors, div);
};

GLmol.prototype.drawMainchainTube = function(group, atomlist, atomName) {
   var points = [], colors = [], radii = [];
   var currentChain, currentResi;
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == atomName) && !atom.hetflag) {
         if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
            this.drawSmoothTube(group, points, colors, radii);
            points = []; colors = []; radii = [];
         }
         points.push(new THREE.Vector3(atom.x, atom.y, atom.z));
         radii.push((atom.b > 0) ? atom.b / 100 : 0.3); // CHECK: is linear scaling correct?
         colors.push(atom.color);
         currentChain = atom.chain;
         currentResi = atom.resi;
      }
   }
   this.drawSmoothTube(group, points, colors, radii);
};

GLmol.prototype.drawStrip = function(group, points1, points2, colors, div) {
   if ((points1.length) < 2) return;
   div = div || this.axisDIV;
   points1 = this.subdivide(points1, div);
   points2 = this.subdivide(points2, div);

   var geo = new THREE.Geometry();
   var i;
   for (i = 0; i < points1.length; i++) {
      geo.vertices.push(new THREE.Vertex(points1[i])); // 2i
      geo.vertices.push(new THREE.Vertex(points2[i])); // 2i + 1
   }
   for (i = 1; i < points1.length; i++) {
      var f;
      var diff = new THREE.Vector3().sub(points1[i], points2[i]);
         var f = new THREE.Face4(2 * i, 2 * i + 1, 2 * i - 1, 2 * i - 2);
         f.color = new THREE.Color(colors[Math.round((i - 1)/ div)]);
         geo.faces.push(f);
   }
   geo.computeFaceNormals();
   geo.computeVertexNormals(false);
   var material =  new THREE.MeshLambertMaterial();
   material.vertexColors = THREE.FaceColors;
   var mesh = new THREE.Mesh(geo, material);
   mesh.doubleSided = true;
   group.add(mesh);
};


GLmol.prototype.IcosahedronGeometry = function() {
   if (!this.icosahedron) this.icosahedron = new THREE.IcosahedronGeometry(2);
   return this.icosahedron;
};

GLmol.prototype.drawCylinder = function(group, from, to, radius, color, cap) {
   if (!from || !to) return;

   var midpoint = new THREE.Vector3().add(from, to).multiplyScalar(0.5);
   var color = new THREE.Color(color);

   if (!this.cylinderGeometry) {
      this.cylinderGeometry = new THREE.CylinderGeometry(1, 1, 1, this.cylinderQuality, 1, !cap);
      this.cylinderGeometry.faceUvs = [];
      this.faceVertexUvs = [];
   }
   var cylinderMaterial = new THREE.MeshLambertMaterial({color: color.getHex()});
   var cylinder = new THREE.Mesh(this.cylinderGeometry, cylinderMaterial);
   cylinder.doubleSided = true;
   cylinder.rotation.x = Math.PI / 2;
   cylinder.scale.x = cylinder.scale.z = radius;
   cylinder.scale.y = from.distanceTo(to);
   var cylinder2 = new THREE.Object3D();
   cylinder2.add(cylinder);
   cylinder2.position = midpoint;
   cylinder2.lookAt(from);
   group.add(cylinder2);
};

// FIXME: transition!
GLmol.prototype.drawHelixAsCylinder = function(group, atomlist, radius) {
   var start = null;
   var currentChain, currentResi;

   var others = [];

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;
      if (atom.hetflag) continue;
      if (atom.ss != 'h' || atom.ssend || atom.ssbegin) others.push(atom.serial);
      if (atom.atom != 'CA') continue;

      if (atom.ss == 'h' && atom.ssend /* || currentChain != atom.chain || currentResi + 1 != atom.resi*/) {
              this.drawCylinder(group, new THREE.Vector3(start.x, start.y, start.z), new THREE.Vector3(atom.x, atom.y, atom.z), radius, atom.color, true);
         start = null;
      }
      currentChain = atom.chain;
      currentResi = atom.resi;
      if (start == null && atom.ss == 'h' && atom.ssbegin) start = atom;
   }
   if (start != null) this.drawCylinder(group, new THREE.Vector3(start.x, start.y, start.z), new THREE.Vector3(atom.x, atom.y, atom.z), radius, atom.color);

   this.drawCartoon(group, others);
};

// MEMO: This routine assumes chains and residues are ordered and
//  atoms N, CA, C and O come in this order. Is this OK?
GLmol.prototype.drawCartoon = function(group, atomlist, div) {
   this.drawStrand(group, atomlist, 2, div, true);
};

GLmol.prototype.drawStrand = function(group, atomlist, num, div, fill,  coilWidth, helixSheetWidth) {
   num = num || this.strandDIV;
   div = div || this.axisDIV;
   coilWidth = coilWidth || this.coilWidth;
   helixSheetWidth = helixSheetWidth || this.helixSheetWidth;
   var points = []; for (var k = 0; k < num; k++) points[k] = [];
   var colors = [];
   var currentChain, currentResi, currentCA;
   var prevCO = null, ss=null;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == 'O' || atom.atom == 'CA') && !atom.hetflag) {
         if (atom.atom == 'CA') {
            if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
               for (var j = 0; j < num; j++)
                  this.drawSmoothCurve(group, points[j], 1 ,colors, div);
               if (fill) this.drawStrip(group, points[0], points[num - 1], colors, div);
               var points = []; for (var k = 0; k < num; k++) points[k] = [];
               colors = [];
               prevCO = null; ss = null;
            }
            currentCA = new THREE.Vector3(atom.x, atom.y, atom.z);
            currentChain = atom.chain;
            currentResi = atom.resi;
            ss = atom.ss;
            colors.push(atom.color);
         } else { // O
            var O = new THREE.Vector3(atom.x, atom.y, atom.z);
            O.subSelf(currentCA);
            O.normalize(); // can be omitted for performance
            O.multiplyScalar((ss == 'c') ? coilWidth : helixSheetWidth); 
            if (prevCO != undefined && O.dot(prevCO) < 0) {
               O.negate();
            }
            prevCO = O;
            for (var j = 0; j < num; j++) {
               var delta = -1 + 2 / (num - 1) * j;
               points[j].push(new THREE.Vector3(currentCA.x + prevCO.x * delta, 
                 currentCA.y + prevCO.y * delta, currentCA.z + prevCO.z * delta));
            }                         
         }
      }
   }
   for (var j = 0; j < num; j++)
      this.drawSmoothCurve(group, points[j], 1 ,colors, div);
   if (fill) this.drawStrip(group, points[0], points[num - 1], colors, div);
};

GLmol.prototype.drawNucleicAcidLadderSub = function(geo, lineGeo, atoms, color) {
//        color.r *= 0.9; color.g *= 0.9; color.b *= 0.9;
   if (atoms[0] != undefined && atoms[1] != undefined && atoms[2] != undefined &&
       atoms[3] != undefined && atoms[4] != undefined && atoms[5] != undefined) {
      var baseFaceId = geo.vertices.length;
      for (var i = 0; i <= 5; i++) geo.vertices.push(new THREE.Vertex(atoms[i]));
          geo.faces.push(new THREE.Face3(baseFaceId, baseFaceId + 1, baseFaceId + 2));
          geo.faces.push(new THREE.Face3(baseFaceId, baseFaceId + 2, baseFaceId + 3));
          geo.faces.push(new THREE.Face3(baseFaceId, baseFaceId + 3, baseFaceId + 4));
          geo.faces.push(new THREE.Face3(baseFaceId, baseFaceId + 4, baseFaceId + 5));
          for (var j = geo.faces.length - 4, lim = geo.faces.length; j < lim; j++) geo.faces[j].color = color;
    }
    if (atoms[4] != undefined && atoms[3] != undefined && atoms[6] != undefined &&
       atoms[7] != undefined && atoms[8] != undefined) {
       var baseFaceId = geo.vertices.length;
       geo.vertices.push(new THREE.Vertex(atoms[4]));
       geo.vertices.push(new THREE.Vertex(atoms[3]));
       geo.vertices.push(new THREE.Vertex(atoms[6]));
       geo.vertices.push(new THREE.Vertex(atoms[7]));
       geo.vertices.push(new THREE.Vertex(atoms[8]));
       for (var i = 0; i <= 4; i++) geo.colors.push(color);
       geo.faces.push(new THREE.Face3(baseFaceId, baseFaceId + 1, baseFaceId + 2));
       geo.faces.push(new THREE.Face3(baseFaceId, baseFaceId + 2, baseFaceId + 3));
       geo.faces.push(new THREE.Face3(baseFaceId, baseFaceId + 3, baseFaceId + 4));
       for (var j = geo.faces.length - 3, lim = geo.faces.length; j < lim; j++) geo.faces[j].color = color;
    }
};

GLmol.prototype.drawNucleicAcidLadder = function(group, atomlist) {
   var geo = new THREE.Geometry();
   var lineGeo = new THREE.Geometry();
   var baseAtoms = ["N1", "C2", "N3", "C4", "C5", "C6", "N9", "C8", "N7"];
   var currentChain, currentResi, currentComponent = new Array(baseAtoms.length);
   var color = new THREE.Color(0xcc0000);
   
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;

      if (atom.resi != currentResi || atom.chain != currentChain) {
         this.drawNucleicAcidLadderSub(geo, lineGeo, currentComponent, color);
         currentComponent = new Array(baseAtoms.length);
      }
      var pos = baseAtoms.indexOf(atom.atom);
      if (pos != -1) currentComponent[pos] = new THREE.Vector3(atom.x, atom.y, atom.z);
      if (atom.atom == 'O3\'') color = new THREE.Color(atom.color);
      currentResi = atom.resi; currentChain = atom.chain;
   }
   this.drawNucleicAcidLadderSub(geo, lineGeo, currentComponent, color);
   geo.computeFaceNormals();
   var mat = new THREE.MeshLambertMaterial();
   mat.vertexColors = THREE.VertexColors;
   var mesh = new THREE.Mesh(geo, mat);
   mesh.doubleSided = true;
   group.add(mesh);
};

// MEMO: Seems the number of objects is the limiting factor. Will mesh merging help?
GLmol.prototype.drawNucleicAcidStick = function(group, atomlist) {
   var currentChain, currentResi, start = null, end = null;
   
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;

      if (atom.resi != currentResi || atom.chain != currentChain) {
         if (start != null && end != null)
            this.drawCylinder(group, new THREE.Vector3(start.x, start.y, start.z), 
                              new THREE.Vector3(end.x, end.y, end.z), 0.3, start.color, true);
         start = null; end = null;
      }
      if (atom.atom == 'O3\'') start = atom;
      if (atom.resn == '  A' || atom.resn == '  G' || atom.resn == ' DA' || atom.resn == ' DG') {
         if (atom.atom == 'N1')  end = atom; //  N1(AG), N3(CTU)
      } else if (atom.atom == 'N3') {
         end = atom;
      }
      currentResi = atom.resi; currentChain = atom.chain;
   }
   if (start != null && end != null)
      this.drawCylinder(group, new THREE.Vector3(start.x, start.y, start.z), 
                        new THREE.Vector3(end.x, end.y, end.z), 0.3, start.color, true);
};

GLmol.prototype.drawNucleicAcidLine = function(group, atomlist) {
   var currentChain, currentResi, start = null, end = null;
   var geo = new THREE.Geometry();

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;

      if (atom.resi != currentResi || atom.chain != currentChain) {
         if (start != null && end != null) {
            geo.vertices.push(new THREE.Vertex(new THREE.Vector3(start.x, start.y, start.z)));
            geo.colors.push(new THREE.Color(start.color));
            geo.vertices.push(new THREE.Vertex(new THREE.Vector3(end.x, end.y, end.z)));
            geo.colors.push(new THREE.Color(start.color));
         }
         start = null; end = null;
      }
      if (atom.atom == 'O3\'') start = atom;
      if (atom.resn == '  A' || atom.resn == '  G' || atom.resn == ' DA' || atom.resn == ' DG') {
         if (atom.atom == 'N1')  end = atom; //  N1(AG), N3(CTU)
      } else if (atom.atom == 'N3') {
         end = atom;
      }
      currentResi = atom.resi; currentChain = atom.chain;
   }
   if (start != null && end != null) {
      geo.vertices.push(new THREE.Vertex(new THREE.Vector3(start.x, start.y, start.z)));
      geo.colors.push(new THREE.Color(start.color));
      geo.vertices.push(new THREE.Vertex(new THREE.Vector3(end.x, end.y, end.z)));
      geo.colors.push(new THREE.Color(start.color));
    }
   var mat =  new THREE.LineBasicMaterial({linewidth: 1, linejoin: false});
   mat.linewidth = 1.5; mat.vertexColors = true;
   var line = new THREE.Line(geo, mat);
   line.type = THREE.Lines;
   group.add(line);
};

GLmol.prototype.drawCartoonNucleicAcid = function(group, atomlist, div) {
   this.drawStrandNucleicAcid(group, atomlist, 2, div, true);
};

GLmol.prototype.drawStrandNucleicAcid = function(group, atomlist, num, div, fill, nucleicAcidWidth) {
   nucleicAcidWidth = nucleicAcidWidth || this.nucleicAcidWidth;
   div = div || this.axisDIV;
   num = num || this.nucleicAcidStrandDIV;
   var points = []; for (var k = 0; k < num; k++) points[k] = [];
   var colors = [];
   var currentChain, currentResi, currentO3;
   var prevOO = null;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == 'O3\'' || atom.atom == 'OP2') && !atom.hetflag) {
         if (atom.atom == 'O3\'') { // to connect 3' end. FIXME: better way to do?
            if (currentChain != atom.chain || currentResi + 1 != atom.resi) {               
               if (currentO3) {
                  for (var j = 0; j < num; j++) {
                     var delta = -1 + 2 / (num - 1) * j;
                     points[j].push(new THREE.Vector3(currentO3.x + prevOO.x * delta, 
                      currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
                  }
               }
               if (fill) this.drawStrip(group, points[0], points[1], colors);
               for (var j = 0; j < num; j++)
                  this.drawSmoothCurve(group, points[j], 1 ,colors, div);
               var points = []; for (var k = 0; k < num; k++) points[k] = [];
               colors = [];
               prevOO = null;
            }
            currentO3 = new THREE.Vector3(atom.x, atom.y, atom.z);
            currentChain = atom.chain;
            currentResi = atom.resi;
            colors.push(atom.color);
         } else { // OP2
            if (!currentO3) {prevOO == null; continue;} // for 5' phosphate (e.g. 3QX3)
            var O = new THREE.Vector3(atom.x, atom.y, atom.z);
            O.subSelf(currentO3);
            O.normalize().multiplyScalar(nucleicAcidWidth);  // TODO: refactor
            if (prevOO != undefined && O.dot(prevOO) < 0) {
               O.negate();
            }
            prevOO = O;
            for (var j = 0; j < num; j++) {
               var delta = -1 + 2 / (num - 1) * j;
               points[j].push(new THREE.Vector3(currentO3.x + prevOO.x * delta, 
                 currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
            }
            currentO3 = null;
         }
      }
   }
   if (currentO3) {
      for (var j = 0; j < num; j++) {
         var delta = -1 + 2 / (num - 1) * j;
         points[j].push(new THREE.Vector3(currentO3.x + prevOO.x * delta, 
           currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
      }
   }
   if (fill) this.drawStrip(group, points[0], points[1], colors); 
   for (var j = 0; j < num; j++)
      this.drawSmoothCurve(group, points[j], 1 ,colors, div);
};

GLmol.prototype.drawDottedLines = function(group, points, color) {
    var geo = new THREE.Geometry();
    var step = 0.3;

    for (var i = 0, lim = Math.floor(points.length / 2); i < lim; i++) {
        var p1 = points[2 * i], p2 = points[2 * i + 1];
        var delta = p2.clone().subSelf(p1);
        var dist = delta.length();
        delta.normalize().multiplyScalar(step);
        var jlim =  Math.floor(dist / step);
        for (var j = 0; j < jlim; j++) {
           var p = new THREE.Vector3(p1.x + delta.x * j, p1.y + delta.y * j, p1.z + delta.z * j);
           geo.vertices.push(new THREE.Vertex(p));
        }
        if (jlim % 2 == 1) geo.vertices.push(new THREE.Vertex(p2));
    }

    var mat = new THREE.LineBasicMaterial({'color': color.getHex()});
    mat.linewidth = 2;
    var line = new THREE.Line(geo, mat);
    line.type = THREE.Lines;
    group.add(line);
};

GLmol.prototype.getAllAtoms = function() {
   var ret = [];
   for (var i in this.atoms) {
      ret.push(this.atoms[i].serial);
   }
   return ret;
};

// Probably I can refactor using higher-order functions.
GLmol.prototype.getHetatms = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.removeSolvents = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.resn != 'HOH') ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getProteins = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (!atom.hetflag) ret.push(atom.serial);
   }
   return ret;
};

// TODO: Test
GLmol.prototype.excludeAtoms = function(atomlist, deleteList) {
   var ret = [];
   var blackList = new Object();
   for (var _i in deleteList) blackList[deleteList[_i]] = true;

   for (var _i in atomlist) {
      var i = atomlist[_i];

      if (!blackList[i]) ret.push(i);
   }
   return ret;
};

GLmol.prototype.getSidechains = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (atom.atom == 'C' || atom.atom == 'O' || atom.atom == 'N') continue;
      ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getAtomsWithin = function(atomlist, extent) {
   var ret = [];

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.x < extent[0][0] || atom.x > extent[1][0]) continue;
      if (atom.y < extent[0][1] || atom.y > extent[1][1]) continue;
      if (atom.z < extent[0][2] || atom.z > extent[1][2]) continue;
      ret.push(atom.serial);      
   }
   return ret;
};

GLmol.prototype.getExtent = function(atomlist) {
   var xmin = ymin = zmin = 9999;
   var xmax = ymax = zmax = -9999;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      xmin = (xmin < atom.x) ? xmin : atom.x;
      ymin = (ymin < atom.y) ? ymin : atom.y;
      zmin = (zmin < atom.z) ? zmin : atom.z;
      xmax = (xmax > atom.x) ? xmax : atom.x;
      ymax = (ymax > atom.y) ? ymax : atom.y;
      zmax = (zmax > atom.z) ? zmax : atom.z;
   }
   return [[xmin, ymin, zmin], [xmax, ymax, zmax]];
};

GLmol.prototype.getResiduesById = function(atomlist, resi) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (resi.indexOf(atom.resi) != -1) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getResidueBySS = function(atomlist, ss) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (ss.indexOf(atom.ss) != -1) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getChain = function(atomlist, chain) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (chain.indexOf(atom.chain) != -1) ret.push(atom.serial);
   }
   return ret;
};

// for HETATM only
GLmol.prototype.getNonbonded = function(atomlist, chain) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag && atom.bonds.length == 0) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.colorByAtom = function(atomlist, colors) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      var c = colors[atom.elem];
      if (c == undefined) c = this.ElementColors[atom.elem];
      if (c == undefined) c = this.defaultColor;
      atom.color = c;
   }
};


// MEMO: Color only CA. maybe I should add atom.cartoonColor.
GLmol.prototype.colorByStructure = function(atomlist, helixColor, sheetColor, colorSidechains) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (!colorSidechains && (atom.atom != 'CA' || atom.hetflag)) continue;
      if (atom.ss[0] == 's') atom.color = sheetColor;
      else if (atom.ss[0] == 'h') atom.color = helixColor;
   }
};

GLmol.prototype.colorByBFactor = function(atomlist, colorSidechains) {
   var minB = 1000, maxB = -1000;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (colorSidechains || atom.atom == 'CA' || atom.atom == 'O3\'') {
         if (minB > atom.b) minB = atom.b;
         if (maxB < atom.b) maxB = atom.b;
      }
   }

   var mid = (maxB + minB) / 2;

   var range = (maxB - minB) / 2;
   if (range < 0.01 && range > -0.01) return;
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (colorSidechains || atom.atom == 'CA' || atom.atom == 'O3\'') {
         var color = new THREE.Color(0);
         if (atom.b < mid)
            color.setHSV(0.667, (mid - atom.b) / range, 1);
         else
            color.setHSV(0, (atom.b - mid) / range, 1);
         atom.color = color.getHex();
      }
   }
};

GLmol.prototype.colorByChain = function(atomlist, colorSidechains) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (colorSidechains || atom.atom == 'CA' || atom.atom == 'O3\'') {
         var color = new THREE.Color(0);
         color.setHSV((atom.chain.charCodeAt(0)) % 15 / 15.0, 1, 0.9);
         atom.color = color.getHex();
      }
   }
};

GLmol.prototype.colorByResidue = function(atomlist, residueColors) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      c = residueColors[atom.resn]
      if (c != undefined) atom.color = c;
   }
};

GLmol.prototype.colorAtoms = function(atomlist, c) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      atom.color = c;
   }
};

GLmol.prototype.colorByPolarity = function(atomlist, polar, nonpolar) {
   var polarResidues = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS'];
   var nonPolarResidues = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TYR', 'TRP'];
   var colorMap = {};
   for (var i in polarResidues) colorMap[polarResidues[i]] = polar;
   for (i in nonPolarResidues) colorMap[nonPolarResidues[i]] = nonpolar;
   this.colorByResidue(atomlist, colorMap);   
};

// TODO: Add near(atomlist, neighbor, distanceCutoff)
// TODO: Add expandToResidue(atomlist)

GLmol.prototype.colorChainbow = function(atomlist, colorSidechains) {
   var cnt = 0;
   var atom, i;
   for (i in atomlist) {
      atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if ((colorSidechains || atom.atom != 'CA' || atom.atom != 'O3\'') && !atom.hetflag)
         cnt++;
   }

   var total = cnt;
   cnt = 0;
   for (i in atomlist) {
      atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if ((colorSidechains || atom.atom != 'CA' || atom.atom != 'O3\'') && !atom.hetflag) {
         var color = new THREE.Color(0);
         color.setHSV(240.0 / 360 * cnt / total, 1, 0.9);
         atom.color = color.getHex();
         cnt++;
      }
   }
};

GLmol.prototype.drawSymmetryMates = function(group, atomlist, matrices) {
   if (matrices == undefined) return;

   for (var i = 0; i < matrices.length; i++) {
      var mat = matrices[i];
      if (mat == undefined || mat.isIdentity()) continue;
      var symmetryMate = new THREE.Object3D();
      this.drawMainchainCurve(symmetryMate, atomlist, this.curveWidth, 'CA');
      this.drawMainchainCurve(symmetryMate, atomlist, this.curveWidth, '03\'');
      symmetryMate.matrixAutoUpdate = false;
      symmetryMate.matrix = mat;
      group.add(symmetryMate);
   }
};


GLmol.prototype.drawSymmetryMatesWithTranslation = function(group, atomlist, matrices) {
   if (matrices == undefined) return;

   for (var i = 0; i < matrices.length; i++) {
      var mat = matrices[i];
      if (mat == undefined) continue;

      for (var a = -1; a <=0; a++) {
         for (var b = -1; b <= 0; b++) {
             for (var c = -1; c <= 0; c++) {
                var symmetryMate = new THREE.Object3D();
                this.drawMainchainCurve(symmetryMate, atomlist, this.curveWidth, 'CA');
                this.drawMainchainCurve(symmetryMate, atomlist, this.curveWidth, 'O3\'');
                symmetryMate.matrixAutoUpdate = false;
                var translationMat = new THREE.Matrix4().setTranslation(
                   this.protein.ax * a + this.protein.bx * b + this.protein.cx * c,
                   this.protein.ay * a + this.protein.by * b + this.protein.cy * c,
                   this.protein.az * a + this.protein.bz * b + this.protein.cz * c);
                symmetryMate.matrix = mat.clone().multiplySelf(translationMat);
                if (symmetryMate.matrix.isIdentity()) continue;
                group.add(symmetryMate);
             }
         }
      }
   }
};

GLmol.prototype.defineRepresentation = function() {
   var all = this.getAllAtoms();
   var hetatm = this.removeSolvents(this.getHetatms(all));
   this.colorByAtom(all, {});
   this.colorByChain(all);

   this.drawAtomsAsSphere(this.modelGroup, hetatm, this.sphereRadius); 
   this.drawMainchainCurve(this.modelGroup, all, this.curveWidth, 'P');
   this.drawCartoon(this.modelGroup, all, this.curveWidth);
};

GLmol.prototype.getView = function() {
   if (!this.modelGroup) return [0, 0, 0, 0, 0, 0, 0, 1];
   var pos = this.modelGroup.position;
   var q = this.rotationGroup.quaternion;
   return [pos.x, pos.y, pos.z, this.rotationGroup.position.z, q.x, q.y, q.z, q.w];
};

GLmol.prototype.setView = function(arg) {
   if (!this.modelGroup || !this.rotationGroup) return;
   this.modelGroup.position.x = arg[0];
   this.modelGroup.position.y = arg[1];
   this.modelGroup.position.z = arg[2];
   this.rotationGroup.position.z = arg[3];
   this.rotationGroup.quaternion.x = arg[4];
   this.rotationGroup.quaternion.y = arg[5];
   this.rotationGroup.quaternion.z = arg[6];
   this.rotationGroup.quaternion.w = arg[7];
   this.show();
};

GLmol.prototype.setBackground = function(hex, a) {
   a = a | 1.0;
   this.bgColor = hex;
   this.renderer.setClearColorHex(hex, a);
   this.scene.fog.color = new THREE.Color(hex);
};

GLmol.prototype.initializeScene = function() {
   // CHECK: Should I explicitly call scene.deallocateObject?
   this.scene = new THREE.Scene();
   this.scene.fog = new THREE.Fog(this.bgColor, 100, 200);

   this.modelGroup = new THREE.Object3D();
   this.rotationGroup = new THREE.Object3D();
   this.rotationGroup.useQuaternion = true;
   this.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
   this.rotationGroup.add(this.modelGroup);

   this.scene.add(this.rotationGroup);
   this.setupLights(this.scene);
};

GLmol.prototype.zoomInto = function(atomlist, keepSlab) {
   // TODO: expand if symmetry mates are present
   var tmp = this.getExtent(atomlist);
   var center = new THREE.Vector3((tmp[0][0] + tmp[1][0]) / 2, (tmp[0][1] + tmp[1][1]) / 2, (tmp[0][2] + tmp[1][2]) / 2);
   this.modelGroup.position = center.multiplyScalar(-1);
   var x = tmp[1][0] - tmp[0][0], y = tmp[1][1] - tmp[0][1], z = tmp[1][2] - tmp[0][2];
   var maxD = Math.sqrt(x * x + y * y + z * z);
   if (maxD < 25) maxD = 25;

   if (!keepSlab) {
      this.slabNear = -maxD / 1.9;
      this.slabFar = maxD / 3;
   }

   this.rotationGroup.position.z = maxD * 0.35 / Math.tan(Math.PI / 180.0 * this.camera.fov / 2) - 150;
   this.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
};

GLmol.prototype.rebuildScene = function() {
   time = new Date();

   var view = this.getView();
   this.initializeScene();
   this.defineRepresentation();
   this.setView(view);

   console.log("builded scene in " + (+new Date() - time) + "ms");
};

GLmol.prototype.loadMolecule = function(repressZoom) {
   var time = new Date();

   this.protein = {sheet: [], helix: [], scaleMatrix: new THREE.Matrix4().identity(),
                    biomtMatrices: [], symmetryMatrices: [], pdbID: '', title: ''};
   this.atoms = [];

   var source = $('#' + this.id + '_src').val();
   this.parsePDB2(source);
   this.parseSDF(source);
   console.log("parsed in " + (+new Date() - time) + "ms");
   
   var title = $('#' + this.id + '_pdbTitle');
   var titleStr = '';
   if (this.protein.pdbID != '') titleStr += '<a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=' + this.protein.pdbID + '">' + this.protein.pdbID + '</a>';
   if (this.protein.title != '') titleStr += '<br>' + this.protein.title;
   title.html(titleStr);

   this.rebuildScene(true);
   if (repressZoom == undefined || !repressZoom) this.zoomInto(this.getAllAtoms());

   this.show();
 };

GLmol.prototype.setSlabAndFog = function() {
   var center = this.rotationGroup.position.z - this.camera.position.z;
   if (center < 1) center = 1;
   this.camera.near = center + this.slabNear;
   if (this.camera.near < 1) this.camera.near = 1;
   this.camera.far = center + this.slabFar;
   if (this.camera.near + 1 > this.camera.far) this.camera.far = this.camera.near + 1;
   if (this.camera instanceof THREE.PerspectiveCamera) {
      this.camera.fov = this.fov;
   } else {
      this.camera.right = center * Math.tan(Math.PI / 180 * this.fov);
      this.camera.left = - this.camera.right;
      this.camera.top = this.camera.right / this.ASPECT;
      this.camera.bottom = - this.camera.top;
   }
   this.camera.updateProjectionMatrix();
   this.scene.fog.near = this.camera.near + this.fogStart * (this.camera.far - this.camera.near);
//   if (this.scene.fog.near > center) this.scene.fog.near = center;
   this.scene.fog.far = this.camera.far;
};

GLmol.prototype.enableMouse = function() {
   var parentObj = this, glDOM = $(this.renderer.domElement); 

   // TODO: Touch panel support. 
   // Contribution is needed as I don't own any iOS or Android device with WebGL support.
   glDOM.bind('mousedown touchstart', function(ev) {
      ev.preventDefault();
      if (!parentObj.scene) return;
      var x = ev.pageX, y = ev.pageY;
      if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches[0]) {
         x = ev.originalEvent.targetTouches[0].pageX;
         y = ev.originalEvent.targetTouches[0].pageY;
      }
      if (x == undefined) return;
      parentObj.isDragging = true;
      parentObj.mouseButton = ev.which;
      parentObj.mouseStartX = x;
      parentObj.mouseStartY = y;
      parentObj.cq = parentObj.rotationGroup.quaternion;
      parentObj.cz = parentObj.rotationGroup.position.z;
      parentObj.currentModelPos = parentObj.modelGroup.position.clone();
      parentObj.cslabNear = parentObj.slabNear;
      parentObj.cslabFar = parentObj.slabFar;
    });

   glDOM.bind('DOMMouseScroll mousewheel', function(ev) { // Zoom
      ev.preventDefault();
      if (!parentObj.scene) return;
      var scaleFactor = (parentObj.rotationGroup.position.z - parentObj.CAMERA_Z) * 0.85;
      if (ev.originalEvent.detail) { // Webkit
         parentObj.rotationGroup.position.z += scaleFactor * ev.originalEvent.detail / 10;
      } else if (ev.originalEvent.wheelDelta) { // Firefox
         parentObj.rotationGroup.position.z -= scaleFactor * ev.originalEvent.wheelDelta / 400;
      }
      console.log(ev.originalEvent.wheelDelta, ev.originalEvent.detail, parentObj.rotationGroup.position.z);
      parentObj.show();
   });
   glDOM.bind("contextmenu", function(ev) {ev.preventDefault();});
   $('body').bind('mouseup touchend', function(ev) {
      parentObj.isDragging = false;
   });

   glDOM.bind('mousemove touchmove', function(ev) { // touchmove
      ev.preventDefault();
      if (!parentObj.scene) return;
      if (!parentObj.isDragging) return;
      var mode = 0;
      var modeRadio = $('input[name=' + parentObj.id + '_mouseMode]:checked');
      if (modeRadio.length > 0) mode = parseInt(modeRadio.val());

      var x = ev.pageX, y = ev.pageY;
      if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches[0]) {
         x = ev.originalEvent.targetTouches[0].pageX;
         y = ev.originalEvent.targetTouches[0].pageY;
      }
      if (x == undefined) return;
      var dx = (x - parentObj.mouseStartX) / parentObj.WIDTH;
      var dy = (y - parentObj.mouseStartY) / parentObj.HEIGHT;
      var r = Math.sqrt(dx * dx + dy * dy);
      if (mode == 3 || (parentObj.mouseButton == 3 && ev.ctrlKey)) { // Slab
          parentObj.slabNear = parentObj.cslabNear + dx * 100;
          parentObj.slabFar = parentObj.cslabFar + dy * 100;
      } else if (mode == 2 || parentObj.mouseButton == 3 || ev.shiftKey) { // Zoom
         var scaleFactor = (parentObj.rotationGroup.position.z - parentObj.CAMERA_Z) * 0.85; 
         if (scaleFactor < 80) scaleFactor = 80;
         parentObj.rotationGroup.position.z = parentObj.cz - dy * scaleFactor;
      } else if (mode == 1 || parentObj.mouseButton == 2 || ev.ctrlKey) { // Translate
         var scaleFactor = (parentObj.rotationGroup.position.z - parentObj.CAMERA_Z) * 0.85;
         if (scaleFactor < 20) scaleFactor = 20;
         var translationByScreen = new THREE.Vector3(- dx * scaleFactor, - dy * scaleFactor, 0);
         var q = parentObj.rotationGroup.quaternion;
         var qinv = new THREE.Quaternion(q.x, q.y, q.z, q.w).inverse().normalize(); 
         var translation = qinv.multiplyVector3(translationByScreen);
         parentObj.modelGroup.position.x = parentObj.currentModelPos.x + translation.x;
         parentObj.modelGroup.position.y = parentObj.currentModelPos.y + translation.y;
         parentObj.modelGroup.position.z = parentObj.currentModelPos.z + translation.z;
      } else if ((mode == 0 || parentObj.mouseButton == 1) && r != 0) { // Rotate
         var rs = Math.sin(r * Math.PI) / r;
         parentObj.dq.x = Math.cos(r * Math.PI); 
         parentObj.dq.y = 0;
         parentObj.dq.z =  rs * dx; 
         parentObj.dq.w =  rs * dy;
         parentObj.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0); 
         parentObj.rotationGroup.quaternion.multiplySelf(parentObj.dq);
         parentObj.rotationGroup.quaternion.multiplySelf(parentObj.cq);
      }
      parentObj.show();
   });
};


GLmol.prototype.show = function() {
   if (!this.scene) return;

   var time = new Date();
   this.setSlabAndFog();
   this.renderer.render(this.scene, this.camera);
   console.log("rendered in " + (+new Date() - time) + "ms");
};

// For scripting
GLmol.prototype.doFunc = function(func) {
    func(this);
};

return GLmol;
}());
