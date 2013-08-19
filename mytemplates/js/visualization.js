var displayDict = new Object(); //"Associative array", i.e. dict, i.e. object with properties accessible with dict-like syntax
displayDict['displaySolvent'] = true; //By default we do display solvent, this gets flipped to false in toggleSolvent()
displayDict['cartoons'] = false; //We need to keep track of window state so as to not freak the user out.
displayDict['fancyCartoons'] = false;
//By default those three above are always the case. We change them when they're not.
displayDict['ribbons'] = false;
displayDict['ballandstick'] = true;

function toggleSolvent(toggle){ //Utility function for toggling solvent display
  //"toggle" parameter determines whether this actually toggles the solvent, or whether we just flip
  Jmol.script(jmolApplet, "selectionhalos off"); //For selection page safety
  var current_display = displayDict['displaySolvent'];
  if(!(toggle)){
    current_display = !(current_display);
  }else{
    if (current_display){
      Jmol.script(jmolApplet, "display displayed and not [TIP];");
      displayDict['displaySolvent'] = false;
    }else{
      Jmol.script(jmolApplet, "display displayed or [TIP];");
      displayDict['displaySolvent'] = true;
    }
    return;
  }
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
    Jmol.script(jmolApplet, "display [TIP] or displayed;select [TIP];cpk 30%;wireframe 55;color jmol;")
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
  Jmol.script(jmolApplet, 'select none;selectionhalos on');
}
