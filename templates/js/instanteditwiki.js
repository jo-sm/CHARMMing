<!--
//script by yvoschaap.com
//freely useable
//optional link back would be very web 2.0 :)

function datosServidor() {
};
datosServidor.prototype.iniciar = function() {
	try {
		// Mozilla / Safari
		this._xh = new XMLHttpRequest();
	} catch (e) {
		// Explorer
		var _ieModelos = new Array(
		'MSXML2.XMLHTTP.5.0',
		'MSXML2.XMLHTTP.4.0',
		'MSXML2.XMLHTTP.3.0',
		'MSXML2.XMLHTTP',
		'Microsoft.XMLHTTP'
		);
		var success = false;
		for (var i=0;i < _ieModelos.length && !success; i++) {
			try {
				this._xh = new ActiveXObject(_ieModelos[i]);
				success = true;
			} catch (e) {
				// Implementar manejo de excepciones
			}
		}
		if ( !success ) {
			// Implementar manejo de excepciones, mientras alerta.
			return false;
		}
		return true;
	}
}

datosServidor.prototype.ocupado = function() {
	estadoActual = this._xh.readyState;
	return (estadoActual && (estadoActual < 4));
}

datosServidor.prototype.procesa = function() {
	if (this._xh.readyState == 4 && this._xh.status == 200) {
		this.procesado = true;
	}
}

datosServidor.prototype.enviar = function(urlget,datos) {
	if (!this._xh) {
		this.iniciar();
	}
	if (!this.ocupado()) {
		this._xh.open("GET",urlget,false);
		this._xh.send(datos);
		if (this._xh.readyState == 4 && this._xh.status == 200) {
			return this._xh.responseText;
		}
		
	}
	return false;
}


var urlBase = "/charmming/wiki/editmessage/";
var formVars = "";
var changing = false;


function fieldEnter(campo,evt,idfld) {
	evt = (evt) ? evt : window.event;
	if (evt.keyCode == 13) {
	        //The below if condition is used to see if the user
		//enters a blank field for the title
                if(campo.value == "")
                   campo.value = "Untitled";
		elem = document.getElementById( idfld );
		remotos = new datosServidor;
		nt = remotos.enviar(urlBase + "?fieldname=" +escape(elem.id)+ "&content="+escape(campo.value)+"&"+formVars,"");
		//remove glow
		noLight(elem);
		elem.innerHTML = nt;
		changing = false;
		return false;
	} else {
		return true;
	}


}

function fieldBlur(campo,idfld) {
	if (true) {
                if(campo.value == "")
                   campo.value = "Untitled";
		elem = document.getElementById( idfld );
		remotos = new datosServidor;
		//THIS IS WIKI SPECIFIC
	        //Check to make sure the title is unique to the board
		if(elem.id.startsWith('edittitle_'))
		{
		    same_names = document.getElementById('edittitle_' + campo.value);
		  i = 1;
                 if(same_names && elem != same_names)
		 {
		  while(same_names != null)
		   {
		     i++;
		     same_names = document.getElementById('edittitle_' + campo.value + '(' + i + ')');
		   }
		   campo.value = campo.value + '(' + i + ')';
		 }
		}
		nt = remotos.enviar(urlBase + "?fieldname=" +escape(elem.id)+ "&content="+escape(campo.value)+"&"+formVars,"");
		//THIS IS WIKI SPECIFIC
		//If the title is being changed then the id for the message and title has to change also
		//The rankings id must change too
		//So must the table name and rating na,e
		if(elem.id.startsWith('edittitle_'))
		{
                    //the title and url are necessary for the delete function to work
                    //if the elem has edittitle in it then it means the title is whats being edited
                        //if this is the case then campo.value is the new title
                    title = elem.id.replace('edittitle_','');
		    titleurlhref = document.getElementById('url_' + elem.id.replace('edittitle_','')).href;
		    document.getElementById('url_' + elem.id.replace('edittitle_','')).href = titleurlhref.replace(title,campo.value);
		    document.getElementById('url_' + elem.id.replace('edittitle_','')).id = 'url_' + campo.value;
		    document.getElementById('message_' + elem.id.replace('edittitle_','')).id = ('message_' + campo.value);
		    document.getElementById(elem.id.replace('edittitle_','') + "_stars").id = campo.value + "_stars";
		    document.getElementById(elem.id.replace('edittitle_','') + "_currate").id = campo.value + "_currate";
		    for(i=1; i<6; i++)
		    {
		     document.getElementById(elem.id.replace('edittitle_','') + '_' + i).id = campo.value + '_' + i;
		    }
		    elem.id = "edittitle_" + campo.value;
		}
		elem.innerHTML = nt;
		changing = false;
		return false;
	}
}

//edit field created
function cambia(actual) {
	if(!changing){
                //WIKI SPECIFIC
	        //Remove <br>s from appearing if there is a newline
                actual.innerHTML = actual.innerHTML.replace(/<br>/g,"");
		//width = widthEl(actual.id) + 2;
		width = 150;
		height =heightEl(actual.id) + 40;
		//height = 250;
		if(height < 40)
			actual.innerHTML = "<input id=\""+ actual.id +"_field\" style=\"font-family:arial; width: "+width+"px; height: "+height+"px;\" maxlength=\"254\" type=\"text\" value=\"" + actual.innerHTML + "\" onkeypress=\"return fieldEnter(this,event,'" + actual.id + "')\" onfocus=\"highLight(this);\" onblur=\"noLight(this); return fieldBlur(this,'" + actual.id + "');\" />";
		else
			actual.innerHTML = "<textarea name=\"textarea\" id=\""+ actual.id +"_field\" style=\"font-size:1.21em;width: "+width+"px; height: "+height+"px;\" onfocus=\"highLight(this);\" onblur=\"noLight(this); return fieldBlur(this,'" + actual.id + "');\">" + actual.innerHTML + "</textarea>";
	
		changing = true;
	}

		actual.firstChild.focus();
}


//find all span tags with class editText and id as fieldname parsed to update script. add onclick function
function editbox_init(){
	if (!document.getElementsByTagName){ return; }
	var spans = document.getElementsByTagName("span");

	// loop through all span tags
	for (var i=0; i<spans.length; i++){
		var spn = spans[i];

        	if (((' '+spn.className+' ').indexOf("editText") != -1) && (spn.id)) {
			spn.onclick = function () { cambia(this); }
			spn.style.cursor = "pointer";
			spn.title = "Click to edit!";	
       		}

	}


}

//crossbrowser load function
function addEvent(elm, evType, fn, useCapture)
{
  if (elm.addEventListener){
    elm.addEventListener(evType, fn, useCapture);
    return true;
  } else if (elm.attachEvent){
    var r = elm.attachEvent("on"+evType, fn);
    return r;
  } else {
    alert("Please upgrade your browser to use full functionality on this page");
  }
}

//get width of text element
function widthEl(span){

	if (document.layers){
	  w=document.layers[span].clip.width;
	} else if (document.all && !document.getElementById){
	  w=document.all[span].offsetWidth;
	} else if(document.getElementById){
	  w=document.getElementById(span).offsetWidth;
	}
return w;
}

//get height of text element
function heightEl(span){

	if (document.layers){
	  h=document.layers[span].clip.height;
	} else if (document.all && !document.getElementById){
	  h=document.all[span].offsetHeight;
	} else if(document.getElementById){
	  h=document.getElementById(span).offsetHeight;
	}
return h;
}

function highLight(span){
            span.parentNode.style.border = "2px solid #D1FDCD";
            span.parentNode.style.padding = "0";
            span.style.border = "1px solid #54CE43";
}

function noLight(span){
        span.parentNode.style.border = "0px";
        span.parentNode.style.padding = "2px";
        span.style.border = "0px";     
	
}

//sets post/get vars for update

function setVarsForm(vars){
	formVars  = vars;
}

addEvent(window, "load", editbox_init);
-->
