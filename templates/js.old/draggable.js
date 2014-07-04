var firefox = document.getElementById && !document.all;
var x,y;
var tx,ty;
var isDrag=false;
var uniqueID = 0;
var currwid = null;
var clickstate = currwid;

document.onmousemove=mouseMove;
document.onmousedown=selectMouse;
document.onmouseup=function(){
 isDrag=false;
}

function newState(widgetname) //When there are multiple windows, this will keep track of which window is currently selected
{
 clickstate = widgetname;
}

function newMessage(message, drag, leftpos, toppos){
 var barId;
 if(!drag)
 {
  barId = "nodrag";
 }
 else
 {
  barId = "titlebar";
 }
 currwid = "widget"+uniqueID;
 clickstate = currwid;
 uniqueID++;
 var ref=document.getElementById("closewid");		//find the place to add widget to
 var wstr="<div onmousedown=javascript:newState('"+currwid+"') id='"+currwid+"' class='transON' onmouseover=this.className='transOFF' onmouseout=this.className='transON' style='position:absolute; top: "+toppos+"; left: "+leftpos+" ; width: 150px;'>"+
		"<div id='wdrag' style='";
 if(drag)
 {
  wstr+=        "cursor: move;";
 }
 wstr+=         "width: 100%; height: 16; "+
		"background-color: #FFC859; border-bottom: 1px solid #CACACA; '>"+
		"<table id='d' style='font-family: arial; font-size: 8pt;' width=100% cellpadding=0 cellspacing=0 border=0>"+
		"<tr><td id='"+barId+"' ALIGN=center><span style='padding-left: 5px;'>PDB Status</span>"+
		"</td><td align=right> <span style='padding-right: 5px;'>" +
		"<a href=javascript:closeWidget('"+currwid+"') style='text-decoration: none;'>"+
		"<b>X</b></a></span></td></tr></table></div>"+
		"<div id='id1_content' style='padding: 5px; font-family: arial; font-size: 8pt; width: 100%'>"+
		message + " </div></div>";
	ref.innerHTML+=wstr;	//place the widget on screen
}


function selectMouse(e)
{
 if(firefox)
 {
  var p=e.target;
  if(p.attributes['id'] && p.attributes['id'].value=="titlebar")
  {
   isDrag=true;
   x=e.clientX;
   y=e.clientY;
   tx=parseInt(document.getElementById(clickstate).style.left);
   ty=parseInt(document.getElementById(clickstate).style.top);
  }
 }
  else
  {
   var p=event.srcElement;
   if(p.attributes['id'] && p.attributes['id'].value=="titlebar")
   {
    isDrag=true;
    x=event.clientX;
    y=event.clientY;
    tx=parseInt(document.getElementById(clickstate).style.left);
    ty=parseInt(document.getElementById(clickstate).style.top);
   }
  }
}

function mouseMove(e)
{
 if(isDrag)
 {
  var box=document.getElementById(clickstate);
  if(firefox)
  {
   box.style.left = e.clientX + (tx-x);
   box.style.top = e.clientY + (ty-y);
  }
  else
  {
   box.style.left = event.clientX + (tx-x);
   box.style.top = event.clientY + (ty-y);
  }
 }
}

function closeWidget(id){	
 var widget1=document.getElementById('closewid');
 var node=document.getElementById(id);
 widget1.removeChild(node);
}

