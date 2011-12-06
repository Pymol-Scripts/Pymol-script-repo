// JavaScript to create a visualization page for APBS web service
// Bob Hanson 1:17 PM 9/23/2010; last updated 9:26 AM 9/24/2010

// APPLET_PATH should be defined in the HTML, not here. 

if(typeof(APPLET_PATH) == "undefined") APPLET_PATH = "http://chemapps.stolaf.edu/jmol/docs/examples-12"
if(typeof(GZIP) == "undefined") GZIP = "" // better: ".gz"


// This empty GIF image is required for the scroller to avoid a thin line in MSIE:

var TRANSP_GIF = APPLET_PATH + "/apbstransp.gif" 

// I think this should take care of the temp file path:

var FILE_PATH = (document.location.href.indexOf("kryptonite.nbcr.net") >= 0 ? "/pdb2pqr/tmp" 
		: document.location.href.indexOf("examples-12") >= 0 ? "http://cupid.wustl.edu/~yonghuang/pdb2pqr/tmp" 
		: "tmp");


////// testing only //////

var downSampleTest = false  // rough surface for testing -- only 1/4 of the points in each dimension
var TEST_JOB = "12852665591" // Bob's private apbsjmol.htm uses this

// CAPITALIZED words here will be replaced by the fixScripts() function

var LOAD_SCRIPT = "load PATH/JOBID/JOBID.pqr;"
var ISOSURFACE_SCRIPT = "set echo top left;echo loading surface data...; refresh; set solventProbeRadius 1.4;\
		isosurface s1 " + (downSampleTest ? "downsample 4 " : "") 
		+ "colorscheme \"SCHEME\" color absolute MIN MAX SURFACE map 'PATH/JOBID/JOBID-pot.dxGZIP';echo \"\";";

////// option state globals ///////

var thisJobId, thisMin, thisMax, thisScheme, thisSurface, min0, max0, scheme0


function fixScript(s) {
	s = s.replace(/PATH/g,FILE_PATH).replace(/JOBID/g,thisJobId).replace(/GZIP/, GZIP)
		.replace(/SURFACE/, thisSurface).replace(/SCHEME/, thisScheme)
		.replace(/MIN/, thisMin).replace(/MAX/, thisMax)
	return s
}

////////// page setup ////////////

function init() {

  // not used; could be linked to via <body onload=init()>

}

function createVisualization(jobid, min, max) {

	// code to create the applet and buttons 

	jmolInitialize(APPLET_PATH, "JmolAppletSigned0.jar")
	jmolSetDocument(0)  // return HTML, don't do document.write

	if (jobid == "TEST") jobid = TEST_JOB

	document.title = "APBS Visualization " + jobid

	thisJobId = jobid
	thisMin = min0 = min
	thisMax = max0 = max
	thisScheme = scheme0 = "bwr"
	thisSurface = "sasurface"

	var script = LOAD_SCRIPT + ISOSURFACE_SCRIPT + "javascript \"createIsosurfaceOptions();createDisplayOptions()\""

	var s = "<h3><a target=_blank href=http://www.poissonboltzmann.org/apbs/>APBS</a> Visualization <a target=_blank href=" + FILE_PATH + "/" + jobid + ">" + jobid + "</a></h3>"
		+ "<table><tr><td width=500>"
		+ jmolApplet(500,fixScript(script))
		+ "<br /><small>powered by <a target='_blank' href='http://jmol.org'>Jmol</a> <span id=versiondiv></span></small>"
		+ "</td><td valign='top' width=350>"
		+ "<div id=div1 style='width:280px;text-align:left'></div>"
		+ "<div id=div2 style='width:350px;text-align:left'></div>"
		+ "</td></tr></table>"

	document.write(s)

}

// div1 is created by a callback from the applet when file loading is complete.

function createIsosurfaceOptions() {
	document.getElementById("versiondiv").innerHTML = jmolEvaluate("getproperty('appletinfo.version')") + " <a href=\"javascript:jmolScript('console')\">console</a>"

	var s = "<table>\
		<tr><td>surface:</td><td colspan='2'><select id=selSurface onkeypress='doSurface()' onchange='doSurface()' style='width:140pt'>\
		<option value='vdw'"
			+ (thisSurface == "vdw" ? "selected='true'" : "") + ">van der Waals (fastest)</option>\
		<option value='sasurface' " 
			+ (thisSurface == "sasurface" ? "selected='true'" : "") + ">solvent-accessible</option>\
		<option value='solvent'"
			+ (thisSurface == "solvent" ? "selected='true'" : "") + ">solvent-excluded (slowest)</option>\
		</select></td></tr>\
		<tr><td>scheme:</td><td colspan='2'><select id=selScheme onkeypress='doRemap(1)' onchange='doRemap(1)' style='width:140pt'>\
		<option value='bwr'"
			+ (thisScheme == "bwr" ? "selected='true'" : "") + ">blue-white-red</option>\
		<option value='rwb'"
			+ (thisScheme == "rwb" ? "selected='true'" : "") + ">red-white-blue</option>\
		<option value='rgb'"
			+ (thisScheme == "rgb" ? "selected='true'" : "") + ">red-green-blue</option>\
		</select></td></tr>\
		<tr><td valign=bottom>minimum:</td><td valign=bottom><input id=txtmin value='" + thisMin + "' style='width:20pt;text-align:right' onkeypress='doRemap(1)'> kT/e\
			<img src=" + TRANSP_GIF + " height=30 width=1></td><td>" 
			+ newScroller("min","","doScroll",100,-1,-1,false,-50,50,thisMin,1, "mouseup") + "</td></tr>\
		<tr><td valign=bottom>maximum:</td><td valign=bottom><input id=txtmax value='" + thisMax + "' style='width:20pt;text-align:right' onkeypress='doRemap(1)'> kT/e\
			<img src=" + TRANSP_GIF + " height=30 width=1></td><td>" 
			+ newScroller("max","","doScroll",100,-1,-1,false,-50,50,thisMax,1, "mouseup") + "</td></tr>\
		<tr><td></td><td><a href='javascript:resetDefaults()'>reset</a></td></tr>\
		</table>"
	document.getElementById("div1").innerHTML = s
	initScrollers()

}

// div2 is also created by a callback from the applet when file loading is complete.

function createDisplayOptions() {

	_jmol.buttonCssText="style='width:130'"
	var s = ""
	s += "<br />surface: "
	s += "<br />"
	s += jmolButton("isosurface on", "on")
	s += jmolButton("isosurface off", "off")
	s += jmolButton("isosurface translucent", "translucent")
	s += jmolButton("isosurface opaque", "opaque")
	s += "<br />model: (requires <b>off</b> or <b>translucent</b>)"
	s += "<br />"
	s += jmolButton("color structure")
	s += jmolButton("color group")
	s += "<br />"
	s += jmolButton("color amino")
	s += jmolButton("color cpk")
	s += "<br />"
	s += jmolButton("trace only")
	s += jmolButton("cartoon only")
	s += "<br />"
	s += jmolButton("backbone only")
	s += jmolButton("spacefill only;spacefill 23%;wireframe 0.15","ball&stick")
	s += "<br />save/load options: <a href=javascript:doSaveNote()>note</a>"
	s += "<br />"
	s += jmolButton("write IMAGE ?.png","save image")
	s += jmolButton("zap;load ?.png","load image")
	s += "<br />"
	s += jmolButton("write ?.jmol","save Jmol")
	s += jmolButton("zap;load ?.jmol","load Jmol")
	document.getElementById("div2").innerHTML = s
}

function doSaveNote() {
	var s = "A .jmol file is simply a ZIP file that can be dragged back into this applet or the Jmol app to recreate the exact state when it was created. "
	+ "It contains all the files necessary -- the .pqr file and the APBS .dx file. \n\nSaving the image is faster, "
	+ "and it also can be loaded back into Jmol, but that will require the remote files sill being present, because, while the image encodes the state, it does not include the data files themselves."
	alert(s)
}

////////// event processing ////////////


function resetDefaults() {
	thisMin = min0
	thisMax = max0
	thisScheme = scheme0
	createIsosurfaceOptions()
	doRemap()
}

function doSurface() {

	// called by selecting a new surface

	var d = document.getElementById("selSurface")
	var s = d[d.selectedIndex].value
	if (s == thisSurface) return
	thisSurface = s
	jmolScript("isosurface delete;refresh;" + fixScript(ISOSURFACE_SCRIPT))
	createIsosurfaceOptions() 
}


function doScroll(name, value) {

	// callback from the scroller

	var v = parseInt(value)
	if (name == "min") {
		if (v == thisMin)return
		thisMin = v
	} else {
		if (v == thisMax)return
		thisMax = v
	}
	document.getElementById("txtmin").value = thisMin
	document.getElementById("txtmax").value = thisMax
	doRemap()
}
	
function doRemap(doDelay) {

	// issue the Jmol COLOR ISOSURFACE command to remap the isosurface
	// when called from a key press, this needs to have a short timeout.

	if (doDelay) {
		setTimeout("doRemap()", 100)
		return
	}
	var d = document.getElementById("selScheme")
	thisScheme = d[d.selectedIndex].value
	var script = "color isosurface '"+thisScheme+"' range " + document.getElementById("txtmin").value + " " + document.getElementById("txtmax").value
	jmolScriptWait(script)
}


//////////

//CHECKJS  E:\js\gca\site\scroller.js 8/18/2008 6:51:31 PM
// scroller.js gives simple div-based scroll
// Bob Hanson 9:22 PM 7/21/2005
// VERY simple div-based scroller.
// body needs onload=initScrollers()
//   then, somewhere: document.write(newScroller())
//   with optional parameters: newScroller(name,caption,fCallback,width,x,y,isvertical,minvalue,maxvalue,initialvalue,factor)
// fCallback must be string form of function, using _n, _v, _p (name, value, position)
// For example, "newChart(Chart,'_n','_v','_p')
// to set a scroller value from code, use resetScroller()
Scrollers={}
isScrollerInitialized=false

function checkScroll(name){
	var pos=0
	var value=0
	var d=0
	var S=Scrollers[name]
	d=document.getElementById("scr_"+name)
	pos=(S.isvertical?d.scrollTop:d.scrollLeft)
	if(pos!=S.pos){
		S.pos=pos //0 to 20*width
		S.value=value=scrollValue(name)
		setScrollCaption(name)
		if(S.isinitialized)setTimeout(strScrollValues(S,S.callback),10)
	}
}

function initScrollers(){
	setTimeout('resetScroller("",-1)',100)
}

function newScroller(name,caption,fCallback,width,x,y,isvertical,minvalue,maxvalue,initialvalue,factor,fmouseup){
	if(!name)name="scroll-test"
	//if(!caption)caption="testing: pos=_p value=_v"
	if(!fCallback)fCallback="testScroll"
	if(!width)width=300
	if(isNaN(x))x=100
	if(isNaN(y))y=100
	if(!isvertical)isvertical=0
	if(!minvalue)minvalue=0
	if(!maxvalue)maxvalue=100
	if(isNaN(initialvalue))initialvalue=(maxvalue + minvalue)/2
	if(!factor)factor=4
	var s=""
	var sout=""
	var S=Scrollers[name]={}
	S.name=name
	S.caption=caption
	S.width=width
	S.x=x
	S.y=y
	S.isvertical=isvertical
	S.maxvalue=maxvalue
	S.minvalue=minvalue
	S.callback=fCallback
	S.value=initialvalue
	S.factor=factor
	S.initialvalue=initialvalue

	// slimmed down version for the APBS application

	var pos = (x < 0 ? "" : "position:absolute:left:" +  (isvertical ? (x-20) + ";top:" + y : x+";top:"+(y-40)))
	sout="\n<div id='scr_"+name+"'  onscroll=\"checkScroll('"+name+"')\" style='font-size:2pt;height:25;width:"+width+";overflow:auto'>"
		+"<img src=" + TRANSP_GIF +" height=1 width="+(width*(factor+1))+">"
		+"</div>"
	S.div=sout
	return sout
}

function resetScroller(which,value){
	for(var name in Scrollers){if(!which||name==which){
			Scrollers[name].isinitialized=false
			setScrollValue(name,0)
			setScrollValue(name,(value>=0?value:Scrollers[name].initialvalue),1)
			Scrollers[name].isinitialized=true
	}}
}

function scrollValue(name){
	var S=Scrollers[name]
	var value=S.pos/(S.width*S.factor)*(S.maxvalue - S.minvalue) + S.minvalue
	if(S.magnitude>1)value=Math.floor(value)
	if(S.magnitude==1)value=Math.floor(value*100)/100
	return value
}

function setScrollCaption(name){
	var caption=""
	var S=Scrollers[name]
	if(S.caption.indexOf("_")>=0){
		caption=strScrollValues(S,S.caption)
		//document.title=name+" "+caption
		document.getElementById("scr_"+name+"_caption").innerHTML=caption
	}
}

function setScrollPosition(name,pos){
	var d=0
	var S=Scrollers[name]
	d=document.getElementById("scr_"+name)
	if(S.isvertical){
		d.scrollTop=pos
	}else{
		d.scrollLeft=pos
	}
	S.pos=pos
	S.value=scrollValue(name)
	setScrollCaption(name)
}

function setScrollValue(name,value,dotrigger){
	var S=Scrollers[name]
	var pos=Math.floor((value - S.minvalue)/(S.maxvalue - S.minvalue)*S.width*S.factor)
	S.pos=pos
	S.value=scrollValue(name)
	setScrollPosition(name,pos)
	//	if(dotrigger)checkScroll(name)
}

function strScrollValues(S,what){
	return (what + "('_n',_v,_p)").replace(/\_n/,S.name).replace(/\_p/,S.pos).replace(/\_v/,S.value)
}

function testScroll(name,value,position){
	document.title="this is function testScroll('"+name+"',"+value+","+position+")"
}

