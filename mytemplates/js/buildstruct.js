var selection = false; //This keeps track of if you have made a selection, so we don't have to keep going to getPropertyAsArray
var show_hydrogens = true;
var show_all_protein = false; //This is when it's selected...important for showing hydrogens

function change_proto_state(newresbox, oldresid, segid){
  if(show_hydrogens){ //If show_hydrogens is false, why should we go through all this work?
    var segment_number = segmentlist.indexOf(segid);
    var proto_state_identifier = $(newresbox).val(); //i.e. LSN, ASP, etc.
    //We can rely on formalCharge given that everything has hydrogens on it to balance out the valences
    //and formalCharge is updated by JSmol whenever we change the charge of something.
    //However, it should NOT be trusted for ANYTHING else!
    var selection_string = " and atomno >= " + chain_terminators[segment_number];
    if (segment_number < chain_terminators.length){
      selection_string = selection_string + " and atomno < " + chain_terminators[segment_number + 1];
    }
    //Now goes a series of if statements for every residue
    if (proto_state_identifier.toUpperCase() == "LSN" || proto_state_identifier.toUpperCase() == "LYS"){ //Lysine
      var NZ = Jmol.getPropertyAsArray(jmolApplet, "atominfo", ("("+oldresid + ".NZ" + selection_string+")"));
     if(NZ.length < 1){
      //Make a fancy error box
    }else{
      if (NZ.formalCharge > 0){
        if(proto_state_identifier.toUpperCase() == "LSN"){
          //Deprotonate
          Jmol.script(jmolApplet, 'assign atom ({' + oldresid + '.NZ' + selection_string + '}) "Mi"');
        }
        //Otherwise do nothing
      }else{
        if(proto_state_identifier.toUpperCase() == "LYS"){
          //Protonate
          Jmol.script(jmolApplet, 'assign atom ({' + oldresid + '.NZ' + selection_string + '}) "Pl"');
        }
        //If it's LSN and deprotonated, do nothing
      }
    }
  }else if(proto_state_identifier.toUpperCase() == "GLU" || proto_state_identifier.toUpperCase() == "GLUP"){ //Glutamic acid 
    var OE2 = Jmol.getPropertyAsArray(jmolApplet, "atominfo", ("("+oldresid+".OE2"+selection_string+")"));
    if(OE2.length < 1){
      //Fancy error box
    }else{
      if(OE2.formalCharge < 0){
        if(proto_state_identifier.toUpperCase() == "GLUP"){ 
          //Protonate
          Jmol.script(jmolApplet, 'assign atom ({'+oldresid+'.OE2'+selection_string+'}) "Pl"');
        }
        //Do nothing otherwise, GLU is -1 charge
      }else{
        if(proto_state_identifier.toUpperCase() == "GLU"){
          //Deprotonate
          Jmol.script(jmolApplet, 'assign atom ({'+oldresid+'.OE2'+selection_string+'}) "Mi"');
        }
      }
    }
//  }else if(proto_state_identifier.toUpperCase() == "HSD" || proto_state_identifier.toUpperCase() == "HSE" || proto_state_identifier.toUpperCase() == "HSP"){
//  var 
  }else{
    return false;
  }
  return false;
  }
  return false;
}


function select_residue(resid,segid){
  var segment_number = segmentlist.indexOf(segid); //This is a bit trickier than a dict can handle so we need to do this trick
  var selection_string = resid + " and atomno >= " + chain_terminators[segment_number];
  selection = true;
  if (segment_number < chain_terminators.length){
    selection_string = selection_string + " and atomno < " + chain_terminators[segment_number + 1];
  }
  Jmol.script(jmolApplet, "select resno=" + selection_string);
  Jmol.script(jmolApplet, "display within(5.0,selected); center selected; zoom 0;");
  if(!(show_hydrogens)){
    Jmol.script(jmolApplet, "display displayed and not hydrogen;");
  }
  $(".hidden_options").show();
}

function click_proto_box(){
  if($('#proto_box').is(':checked')){
    $('#proton_divs').css('display','inline');
    jmol_height = $("#protonation")[0].scrollHeight - 150;
    if (jmol_height < 400){
      jmol_height = 400;
    }
    //Since you only get scrollHeight once the element renders and if the right side is too small, adapt to JSmol's normal size of 400x400
    $("#protonation_jsmol").css("height",jmol_height);
  }else{
  $('#proton_divs').hide();}
}


function show_protein(){
if(show_all_protein){ //If you're already showing the full thing...
  if(show_hydrogens){
    Jmol.script(jmolApplet, "display within(5.0,selected)");
  }else{
  Jmol.script(jmolApplet, "display within(5.0,selected) and not hydrogen");
  Jmol.script(jmolApplet, "display displayed or within(5.0,selected) and solvent"); //So the solvent hydrogens aren't killed
  }
  show_all_protein = false;
}else{
  show_all_protein = true;
  if(show_hydrogens){
    Jmol.script(jmolApplet, "display all");
  }else{
    Jmol.script(jmolApplet, "display all and not hydrogen;display displayed or solvent;");
  }
}
  
}
function show_hide_hydrogens(){
if(!(show_hydrogens)){ //Here you already have something selected and want to kill off hydrogens
  show_hydrogens = true;
  if(selection){
    if(show_all_protein){
      Jmol.script(jmolApplet, "display all");
    }else{
      Jmol.script(jmolApplet, "display within(5.0,selected)");
    }
  }else{
    Jmol.script(jmolApplet, "display displayed or hydrogen;");
  }
}else{
  if(selection){
    if(show_all_protein){
      Jmol.script(jmolApplet, "display all and not hydrogen; display displayed or solvent;");
    }else{
      Jmol.script(jmolApplet, "display within(5.0,selected) and not hydrogen;display displayed or solvent and within(5.0,selected)");
    }
  }else{
    Jmol.script(jmolApplet, "display displayed and not hydrogen;display displayed or solvent");
  }
  show_hydrogens = false;
}
}
