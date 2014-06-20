if(typeof displayDict == "undefined"){
  var displayDict = new Object(); //"Associative array", i.e. dict, i.e. object with properties accessible with dict-like syntax
}
displayDict['displaySolvent'] = true; //By default we do display solvent, this gets flipped to false in toggleSolvent()
displayDict['cartoons'] = false; //We need to keep track of window state so as to not freak the user out.
displayDict['fancyCartoons'] = false;
//By default those three above are always the case. We change them when they're not.
displayDict['ribbons'] = false;
displayDict['ballandstick'] = true;
displayDict['zoomed_in'] = false; //Tracks all our zoom to region stuff.
$(".background_white, .background_black").click(function(){
    set_background($(this).attr('class').split("_")[1]);
    return false;
  });
$(".display").click(function(){
    changeDisplay($(this).attr('class').split(" ")[1]);
    return false;
  });
$(".solvent").click(function(){
    toggleSolvent(true,displayDict['zoomed_in']);
    return false;
  });
$(".region").click(function(){
    var woof = $(this).attr('class').split(" ")[1];
    if(woof.indexOf("select") > -1){
      selectRegion();
    }else{
      var zoom = woof.replace("zoom","");
      zoomRegion(zoom);
    }
    return false;
});


function toggleSolvent(toggle,zoomed){ //Utility function for toggling solvent display
  //"toggle" parameter determines whether this actually toggles the solvent, or whether we just flip
  var current_display = displayDict['displaySolvent'];
    
  if(!(toggle)){
    current_display = !(current_display);
  }else{
    if (current_display){
      if(zoomed){
        Jmol.script(jmolApplet, "selectionhalos on;"); //make sure to do this always
      }
      Jmol.script(jmolApplet, "display displayed and not [TIP];"); //this is the same both zoomed and not
      if(!(zoomed)){
        Jmol.script(jmolApplet, "selectionhalos off;");
      }
      displayDict['displaySolvent'] = false;
    }else{
      if(zoomed){
        Jmol.script(jmolApplet, "selectionhalos on;display displayed or (within(4.0,within(1,GROUP,within(4.0,selected))) and [TIP]);");
      }else{
        Jmol.script(jmolApplet, "display displayed or [TIP];");
      }
        displayDict['displaySolvent'] = true;
    }
    return;
  }
  Jmol.script(jmolApplet, "selectionhalos off;");
  if (current_display){     
    if(displayDict['cartoons']){ //Only one of these should be true at any time
      Jmol.script(jmolApplet, "display displayed and not [TIP];select displayed;cpk off; wireframe off;ribbons off;cartoons on;set cartoonsfancy false;color structure;");
    }else if(displayDict['fancyCartoons']){
      Jmol.script(jmolApplet, "display displayed and not [TIP];select displayed;cpk off; wireframe off;ribbons off;set cartoonsfancy true;cartoons on;color structure;");
    }else if(displayDict['ribbons']){
      Jmol.script(jmolApplet, "display displayed and not [TIP];select displayed;cpk off; wireframe off;cartoons off;ribbons on;color structure;");
    }else if(displayDict['wireframe']){
      Jmol.script(jmolApplet, "display displayed and not [TIP];select not [TIP] and displayed;cpk off;cartoons off;ribbons off;wireframe 55;color jmol;");
      }else if(displayDict['ballandstick']){
      Jmol.script(jmolApplet, "display displayed and not [TIP];select displayed;cpk 30%;wireframe 55;ribbons off;cartoons off;color jmol;");
    }else{ //What...
      alert("Something has gone very wrong. Please report this message.");
    }
  }else{     
    if(displayDict['cartoons']){ //Only one of these should be true at any time
      Jmol.script(jmolApplet, "select not [TIP] and displayed;cpk off; wireframe off;ribbons off;cartoons on;set cartoonsfancy false;color structure;");
    }else if(displayDict['fancyCartoons']){
      Jmol.script(jmolApplet, "select not [TIP] and displayed;cpk off; wireframe off;ribbons off;set cartoonsfancy true;cartoons on;color structure;");
    }else if(displayDict['ribbons']){
      Jmol.script(jmolApplet, "select not [TIP] and displayed;cpk off; wireframe off;cartoons off;ribbons on;color structure;");
    }else if(displayDict['wireframe']){
      Jmol.script(jmolApplet, "select displayed and not [TIP];cpk off;cartoons off;ribbons off;wireframe 55;color jmol;");
      }else if(displayDict['ballandstick']){
      Jmol.script(jmolApplet, "select displayed and not [TIP];cpk 30%;wireframe 55;ribbons off;cartoons off;color jmol;");
    }else{ //What...
      alert("Something has gone very wrong. Please report this message.");
    }
    Jmol.script(jmolApplet, "display [TIP] or displayed;select [TIP];cpk 30%;wireframe 55;color jmol;select none;"); //select none is IMPORTANT
    //if you don't select none, zoomRegion() will start doing very strange zooms.
  }
}

function changeDisplay(displayVar){ //Nice and elegant.
  for (var x in displayDict){ //Hopefully this doesn't cause issues like most other forloops or attempted forloops on this site
    if (x != 'displaySolvent'){
      displayDict[x] = false;
    }
  }
  displayDict[displayVar] = true;
  toggleSolvent(false); //More generic form...
}

function refresh_select(){ //Fixes selectionhalos
  Jmol.script(jmolApplet, 'select displayed and not [TIP];cartoons on;selectionhalos on');
}

function selectRegion(){ //Does the zoom to region stuff we want for MSCALE
  displayDict['zoomed_in']=true;
  Jmol.script(jmolApplet, 'select all and not [TIP];cpk off;wireframe off;cartoons on;\
      ribbons off;color structure;select [TIP];cpk 30%;wireframe 45;color jmol;select none;selectionhalos on;set picking group;');
}

function zoomRegion(zoomlevel){
  displayDict['zoomed_in']=true;
  displayDict['zoomlevel']=zoomlevel;
  //the following script does the following:
  //first select all atoms that are within zoomlevel angstroms of the residue selected,
  //as well as the residues connected to those atoms, if any. Then make that region into ballandstick mode
  //(cpk, wireframe does this), and turn off ribbons and cartoons for that region, and put up regular Jmol atom coloring.
  //then zoom in to that region (250 is an arbitrary zoom constant, can be anything from 5-200000, 250 seems to work best
  //but we can modify this value depending on zoomlevel), then everything but the solvent and the region in question
  //gets colored in a translucent fashion to improve visibility. Then all polar hydrogens are selected (i.e. those attached
  //to nitrogen/oxygen (element numbers 7 and 8)) and the non-polar hydrogens are removed from the visualization because they
  //clutter up the visual space and are not very relevant for almost any calculations you'd want to do with proteins.
  //Then we turn off solvent everywhere but in a cage (size determined by zoomlevel) around the region we defined before
  //and make the solvent in that area wireframe instead of ball-and-stick
  //TODO: Modify this script so we include Ions, and make sure that those come out as balls rather than wireframe!
  Jmol.script(jmolApplet, 'save selection oldselect;\
      selectionhalos off;select selected or within(1,GROUP,within(' + zoomlevel + ',selected));cpk 20%;wireframe 35;\
      ribbons off;cartoons off;color jmol;zoomto 1 (selected) 250;\
      select not selected and not [TIP];cpk off;wireframe off;ribbons on;color translucent structure;select connected(elemno=7 or elemno=8) and hydrogen;\
      display selected or (not hydrogen);display displayed and not [TIP];restore selection oldselect;\
      select within(' + zoomlevel + ',within(1,GROUP,within(' + zoomlevel + ',selected))) and [TIP];cpk 0;wireframe 35;display displayed or selected;\
      restore selection oldselect;selectionhalos on;');
}

function set_background(color){
  //sets background to whatever color the user defines. RIght now we care about white and black.
  Jmol.script(jmolApplet, 'background ' + color);
  }
