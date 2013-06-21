//Jmol helpful functions  
function off(applet){
  Jmol.script(applet, "set atompicking off");
}

// Note: this does not go true/false...it only changes if you click something else.
function drag(applet){
  Jmol.script(applet, "set atompicking on;set picking dragminimize");
}

function del(applet){
  Jmol.script(applet, "set atompicking on;set picking assignatom_x");
}

function load(applet, molecule){
  Jmol.script(applet, "load /charmming/js/jsmol/" + molecule + ".xyz");
}

function action_save(applet){
  Jmol.script(applet, "save STATE temp");
  alert("Temporary state saved.");
}

function action_restore(applet){
  Jmol.script(applet, "restore STATE temp");
}

function pickAtom(applet, atombutton){
  Jmol.script(applet, "set atompicking on;set picking assignatom_" + atombutton.id);
}

function setbonds(applet, bondbutton){
  if (bondbutton.id == "n"){
    Jmol.script(applet, "set bondpicking false");
  }else{
    Jmol.script(applet, "set atompicking false;set bondpicking true");
    Jmol.script(applet, "set picking assignbond_" + bondbutton.id);
  }
}

function normalFormSubmit(input){
  if(!(input)){
    ligname = document.getElementById("LigName");
    if (ligname.value.length == 0){
      alert("Please input a name for your ligand.");
      return;
    }
    var spec_chars = "!@#$%^&*()+=-[]\\\';,./{}|\":<>?";
    for (var i=0; i < ligname.value.length; i++){
      if (spec_chars.indexOf(ligname.value.charAt(i)) != -1){
        alert("You cannot have any of the following characters in your ligand name:\n" + spec_chars);
        return;
      }
    }
    var MOLdata = Jmol.scriptWaitOutput(jmolApplet, "write PDB");
    if (MOLdata.length == 0){
      alert("Please build a ligand.");
      return;
    }else{
      var regex = new RegExp("\n", "gi");
      document.getElementById("molinfo").value = MOLdata.replace(regex,"\\n"); 
    }
  }
  document.getElementById("ligand_form").submit(); //This sends it and hands it off to django/Python for processing
}


