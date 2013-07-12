//Color constants for JSmol
var cyan = '[x00ffff]';
var magenta = '[xc800c8]';
var green = '[x008000]';


function create_bynum(add, atominfo){ //Universal bynum creation to make things easier
  var atomsele = "bynum " + atominfo[0].atomno;
  var i=1;//Hack because forloops crash firefox via allocation overflow somehow
  var start_range = atominfo[0].atomno;
  var end_range = start_range; // Add 1 every time
  while(i <= atominfo.length){ //Atominfo is sorted! This is important.
    if (add && atominfo[i - 1].color == "[x00ffff]"){
      $("#dialog_bad_add").dialog("open");
      return false;
    }
    if((i < atominfo.length) && (atominfo[i].atomno == (atominfo[i - 1].atomno + 1))){ //IF there's a run, keep looping.
      end_range = end_range + 1;
    }else{
      if(end_range > start_range){ //i.e. if there's a run at all
        if(start_range == atominfo[0].atomno){ //If it's the first...
          atomsele = atomsele + ":" + end_range;
        }else{
          atomsele = atomsele + " .or. bynum " + start_range + ":" + end_range;
        }
      }else{
        atomsele = atomsele + " .or. bynum " + start_range;
      }
      if(i < atominfo.length){
        start_range = atominfo[i].atomno;
        end_range = start_range;
      }
  }
  i = i + 1;
  }
  if (!(add)){ //No selection
    $("#atomselection").val(atomsele);
  }else{ //Add selection, you shouldn't select things that are already in there
    $("#atomselection").val($("#atomselection").val() + " .or. " + atomsele);
  }
}

function submit_atomselect(){
  if ($("#add").is(":visible")){ //Stop playing with our code!
      reset_select();
      return false;
    }
  var atominfo = Jmol.getPropertyAsArray(jmolApplet, "atominfo", "selected");
  if (atominfo.length == 0){
    $("#dialog_noatoms_alert").dialog("open");
    return false;
  }
  create_bynum(false, atominfo); //Create bynum without add
//The logic is a little screwy but it should work.
//  document.getElementById("atomselection").value = atomsele;
  Jmol.script(jmolApplet, "color cyan;select none;");
  Jmol.script(jmolApplet, "set picking atom");
  $(".atomselectbutton").hide(); //Hide all atom selection buttons...
  $("#add").show();
//        document.getElementById("atomselect_form").submit();
}

function add_atomselect(){
  if ($(".atomselectbutton").is(":visible")){ //Hey! Stop it.
      reset_select();
      return false;
  }
  var atominfo = Jmol.getPropertyAsArray(jmolApplet, "atominfo", "selected");
  if (atominfo.length == 0){
    $("#dialog_noatoms_alert").dialog("open");
    return false;
  }
  create_bynum(true, atominfo);
  Jmol.script(jmolApplet, "color cyan; select none; set picking atom;");
}  

function isNumeric(num) {
  return !isNaN(parseInt(num)) && isFinite(num);
  //Returns true if num is numeric, false if anything else. Can even be a string
  //with value "100", for example. Returns false if the field contains special characters
  //or letters, or anything like that
}

function select_qmsele(divtype, divnum){ //divtype holds what type of input it came from, divnum which one
  var atomarray = Jmol.getPropertyAsArray(jmolApplet, "atominfo", "selected");
  if (atomarray.length != 1){ //Oh dear. That's no good.
    $("#dialog_too_many_links").dialog("open");
    return false;
  }else{
    var atom = atomarray[0]; //There's only supposed to be one
    if (atom.name.toUpperCase() != "CA" && atom.name.toUpperCase() != "C" && atom.name.toUpperCase() != "CT2" && atom.name.toUpperCase() != "CB")
/*        && atom.name.toUpperCase() != "CG" && atom.name.toUpperCase() != "CD" && atom.name.toUpperCase() != "CD1" && atom.name.toUpperCase() != "CD2"
        && atom.name.toUpperCase() != "CE")*/{
      $("#dialog_bad_bond").dialog("open"); //At this point, you have an atom selected, so either you unselect it, or you reset the whole thing
      return false;
    }else{
     /* if (atom.resname.toUpperCase() == "PRO" && (atom.name.toUpperCase() != "CA" && atom.name.toUpperCase() != "CB") || 
          atom.resname.toUpperCase() == "TRP" && (atom.name.toUpperCase() != "CA" && atom.name.toUpperCase() != "CB" && atom.name.toUpperCase() != "CG"){
            $("#dialog_bad_bond").dialog("open");
            return false;
        } */
    }
    var atom_segid = "";
    i=0; //Search chain terminators for the atom number you're at
    while ( i < chain_terminators.length){
      if (atom.atomno < chain_terminators[i]){
        atom_segid = segmentlist[i - 1];
        break;
      }else if (i == (chain_terminators.length - 1)){
        atom_segid = segmentlist[i];
        break;
      }
      i = i + 1;
    }
    var divid = "#" + divtype + "_" + divnum;
    var identifier = $(divid).val().split("\t"); //This string "identifies" the atom, thus identifier
    var selection_string = ""; //Holds the Jmol script
    //Need to add in code that checks chain_terminators for stuff
    if (!($(divid).val() == "" || $(divid).val() == null)) {
        var seg_loc = chain_terminators.indexOf(identifier[2]);
        if (!(seg_loc >= 0)){
          $("#dialog_bad_linkatom_input").dialog("open");
          return false;
        }
        selection_string = "select " + selestring[0] + "." + selestring[1] +
          " and atomno >= " + seg_loc  + ((seg_loc < chain_terminators.length) ? (" and atomno < " + (seg_loc + 1)):"");
        //The above mass of a string does a weird thing to get segment identification.
        //Since all we have is the atom number at which a residue starts, we select
        //by residue number, name, and atom number greater than or equal to the starting point of the segment
        //in the field, and less than the starting point of the residue following it (if any).
        //If you already have a selection, reset the color of what you had before so things don't get screwy.
        Jmol.script(jmolApplet, selection_string + ";color cyan");
        Jmol.script(jmolApplet, "select atomno == " + atom.atomno); //Select your atom again
    }
    Jmol.script(jmolApplet, "selectionhalos off");
    if (divtype == "qmhost"){
      if (!(atom.color == cyan || atom.color == green)){ //color == cyan or green (since you can have one atom be three QM "hosts". May be an issue. Try special vals?)
        $("#dialog_bad_qm").dialog("open");
        return false;
      }else{
        if (selection_string != ""){ //i.e. if you already have a selection
          Jmol.script(jmolApplet, selection_string + ";color cyan");
          Jmol.script(jmolApplet, "select atomno == " + atom.atomno);
        }
      }
      Jmol.script(jmolApplet, "color green;select none;");
      $("#mmbutton"+divnum).attr("disabled",false);
    }else if (divtype == "mmhost"){
      if (atom.color == cyan || atom.color == green){ //cyan, "magenta", or green
        $("#dialog_bad_mm").dialog("open");
        return false;
      }
      if (selection_string != ""){
        Jmol.script(jmolApplet, selection_string + ";color jmol;");
        Jmol.script(jmolApplet, "select atomno == " + atom.atomno);
      }
      var mmcoords = atom.coord; //So we store the coordinates of our MM host
      var qminfo = $("#qmhost_"+divnum).val().split("\t"); //Get the number which forms the first half of the qmhost input box
        //If your input is wrong, that's your problem, since it'll generate something bogus somewhere.
      //atom doesn't change - it's still atomarray[0]. 
      Jmol.script(jmolApplet, "select none; select " + qminfo[0] + "." + qminfo[1] );
      var qmhost = Jmol.getPropertyAsArray(jmolApplet, "atominfo", "selected");
      var qmcoords = qmhost[0].coord;
      var dist = distance(qmcoords, mmcoords);
      if (dist > 2.0){
        $("#dialog_long_link").dialog("open");
        return false;
      }
      Jmol.script(jmolApplet, "select atomno == " + atom.atomno + ";color " + magenta + ";select none;");
    }else{
      Jmol.script(jmolApplet, "select none");
    }
    $(divid).val(atom.resno + "\t" + atom.name + "\t" + atom_segid);
    Jmol.script(jmolApplet, "selectionhalos on");
  }
}

//Calculates n-dimensional distance for two arrays of matching coordinates
//e.g. (x, y, z), (x, y, z). If the indices don't match
//e.g (x, y, z), (z, x, y), then the results will be erroneous.
function distance(point1, point2){
  if (point1.length != point2.length){
    console.error("Error in distance function - Dimension vectors not the same length.");
    return false; //this will probably break the page, which is totally fine
  }
  var i = 0;
  var dist = 0;
  try{
    while (i < point1.length){
      dist = dist + Math.pow((point2[i] - point1[i]), 2);
      i = i + 1;
    }
    return Math.sqrt(dist);
  }catch (e){
    console.error("Error in distance function - Cannot calculate distance. Check for nulls.");
    console.error(e);
  }
}
  

function reset_select(){
  Jmol.script(jmolApplet, "selectionhalos off;select all;color jmol;select none;selectionhalos on;");
  $(".atomselectbutton").show(); //This way the user can't touch it till JSmol renders
  $("#atomselection").val("");
  $("#linkatom_num").val("");
  $("#linkatom_inputs").html("");
  $("#add").hide();
}


function submit_qm_mm_atoms(){
  var qm_mm_selectors = $(".qmmminput");
  var i = 1; //I refuse to use for loops anymore
  //i must start at 1 since 0 is the link atom number box
  var numre = new RegExp("[0-9]");
  var carbre = new RegExp("[A-G]");
  var segre = new RegExp("[d,o]");
  //Since we're testing separate parts of the string, let's not combine these regexps
  while (i < qm_mm_selectors.length){
    var selector = qm_mm_selectors[i];
    if (selector.value == "" || selector.value == null){
      $("#dialog_empty_fields").dialog("open");
      return false;
    }
    selector = selector.value.split("\t");
    if (qm_mm_selectors[i].value.length > 17 || numre.exec(selector[0]) == null || carbre.exec(selector[1]) == null || segre.exec(selector[2]) == null){
      $("#dialog_bad_linkatom_input").dialog("open");
      return false;
    }
    i = i + 1;
  }
  $("#atomselect_form").submit();
}


function submit_linkatoms(){
  if ($("#atomselection").val() == ""){ //empty atom selection
    $("#dialog_noatoms_alert").dialog("open");
    return false;
  }
  if (isNumeric($("#linkatom_num").val())){
    var boxenNumber = parseInt($("#linkatom_num").val());
    i=0;
    $("#linkatom_inputs").html("<hr><br>"); //blank it before something goes wrong
    //Let's use while loops again...I don't trust for anymore.
    while(i < boxenNumber){
      $("#linkatom_inputs").append('QM host:<input type="text" readonly class="qmmminput" name="qmhost_' + i + 
          '" id="qmhost_' + i + '" value="">&nbsp;<button id="qmbutton' + i + '" class="selectbutton" onclick="select_qmsele(\'qmhost\', ' + i + ');return false;">Select</button>&nbsp;');
      $("#linkatom_inputs").append('MM host:<input type="text" readonly class="qmmminput" name="mmhost_' + i + 
          '" id="mmhost_' + i + '" value="">&nbsp;<button class="selectbutton" id="mmbutton' + i + '"  disabled onclick="select_qmsele(\'mmhost\', ' + i + ');return false;">Select</button>&nbsp;<br>');
      i = i + 1;
    }
    $("#linkatom_inputs").append('<br><button class="selectbutton" onclick ="submit_qm_mm_atoms();return false;">Submit Selection with Link Atoms</button>');
    $("#linkatom_inputs").css("display", "inline");
  }else{
    $("#dialog_bad_linkatoms").dialog("open");
    return false;
  }
  }

//Incoming jQueryUI error messages

$(function($){
    $("#dialog_long_link").dialog({
      resizable:false,
      height:200,
      width:600,
      modal:true,
      autoOpen:false,
      buttons:{
      "OK":function(){
      Jmol.script(jmolApplet, "select none;selectionhalos on");
      $(this).dialog("close");
      }
    }
  });
});

$(function($){
    $("#dialog_bad_add").dialog({
      resizable:false,
      height:200,
      width:600,
      modal:true,
      autoOpen:false,
      buttons:{
      "Select different atoms":function(){
      Jmol.script(jmolApplet, "select none;selectionhalos on");
      $(this).dialog("close");
      },
      "Reset selection":function(){
      reset_select();
      $(this).dialog("close");
      }
    }
  });
});

$(function($){
    $("#dialog_bad_qm").dialog({
      resizable:false,
      height:200,
      width:600,
      modal:true,
      autoOpen:false,
      buttons:{
      "OK":function(){
        Jmol.script(jmolApplet, "select none;selectionhalos on;");
        $(this).dialog("close");
      }
    }
  });
});

$(function($){
    $("#dialog_bad_mm").dialog({
      resizable:false,
      height:200,
      width:600,
      modal:true,
      autoOpen:false,
      buttons:{
      "OK":function(){
      Jmol.script(jmolApplet, "select none;selectionhalos on");
      $(this).dialog("close");
      }
    }
  });
});

$(function($){
    $("#dialog_empty_fields").dialog({
      resizable:false,
      height:200,
      width:600,
      modal:true,
      autoOpen:false,
      buttons:{
      "OK":function(){
        Jmol.script(jmolApplet, "select none");
        $(this).dialog("close");
      }
    }
  });
});

$(function($){
    $("#dialog_bad_linkatom_input").dialog({
      resizable:false,
      height:200,
      width:600,
      modal:true,
      autoOpen:false,
      buttons:{
      "OK":function(){
        $("#linkatom_inputs").html("");
        $("#linkatom_num").val("");
        $(this).dialog("close");
      }
    }
  });
});

$(function($){
        $( "#dialog_noatoms_alert").dialog({
          resizable:false,
          height:200,
          width:600,
          modal:true,
          autoOpen:false,
          buttons:{
          "OK":function(){
          $(this).dialog("close");
          }
        }
      });
        });
$(function($){
    $("#dialog_too_many_links").dialog({
      resizable:false,
      height:200,
      width:600,
      modal:true,
      autoOpen:false,
      buttons:{
      "OK":function(){
      $(this).dialog("close");
      }
    }
  });
});
$(function($){
    $("#dialog_bad_linkatoms").dialog({
      resizable:false,
      height:200,
      width:600,
      modal:true,
      autoOpen:false,
      buttons:{
      "OK":function(){
      $("#linkatom_num").val("0");
      $(this).dialog("close");
      }
    }
  });
    });
$(function($){
  $("#dialog_bad_bond").dialog({
    resizable:false,
    height:200,
    width:600,
    modal:true,
    autoOpen:false,
    buttons:{
    "Select different link atom":function(){
    Jmol.script(jmolApplet, "select none;selectionhalos on");
    $(this).dialog("close");
    },
    "Reset selection":function(){
    $("#reset").click();
    $(this).dialog("close");
    }
  }
});
});     
