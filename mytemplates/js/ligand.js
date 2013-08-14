var molname = document.getElementById("LigName").value;
molname = molname.substr(0,4).toUpperCase();
var sdf_link = (document.getElementById("dialog_confirm_PDB") != null);
var same_residue = (document.getElementById("dialog_sameres_alert") != null);
var charge_alert = (document.getElementById("dialog_charge_alert") != null);
var struct_error = (document.getElementById("dialog_struct_error") != null);

if(struct_error){
  document.getElementById("dialog_struct_error").innerHTML =  '<p><span class="ui-icon ui-icon alert" style="float:left;margin 0 auto auto 0;"></span>There already exists a ligand named ' + molname + ' in your list of structures.<br />&nbsp;&nbsp;&nbsp;&nbsp;Please choose a different name for your ligand.';
  
  $(function($){
      $( "#dialog_struct_error").dialog({
        resizable:false,
        height:200,
        width:600,
        modal:true,
        buttons:{
        "OK":function(){
        $(this).dialog("close");
        }
      }
    });
    });
}
if(same_residue){
$(function($) {
    $( "#dialog_sameres_alert").dialog({
      resizable:false,
      height:200,
      width:600,
      modal:true,
      buttons:{
      "OK":function(){
        $(this).dialog("close");
        }
    }
  });
});
}


if(charge_alert){
$(function($) { 
    $( "#dialog_charge_alert" ).dialog({
      resizable:false,
      height:200,
      width:600,
      modal:true,
      buttons:{
      "Yes, I'm sure": function() {
        document.getElementById("force_charge").value = "true";
        normalFormSubmit(true);
      },
      "No, I'm not": function() {
        $(this).dialog("close");
        }
      }
    });
    });
}


$(function($){
    $("#dialog_no_name").dialog({
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
    $("#dialog_spec_char_alert").dialog({
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
    $("#dialog_no_ligand").dialog({
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


$(document).ready(function() {
    $("#LigName").keyup(function(event) {
      if (event.keyCode == 13) { 
        normalFormSubmit(false);
      }
    })
  });
    
//JSmol helpful functions  
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
  Jmol.loadFile(applet, "$"+molecule);
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
    var ligname = document.getElementById("LigName");
    if (ligname.value.length == 0){ 
      $("#dialog_no_name").dialog("open");
      return;
    }
    var spec_chars = "!@#$%^&*()+=-[]\\\';,./{}|\":<>?";
    for (var i=0; i < ligname.value.length; i++){
      if (spec_chars.indexOf(ligname.value.charAt(i)) != -1){
        $("#dialog_spec_char_alert p").html('You cannot have any of the following characters in your ligand name:\n' + spec_chars);
        $("#dialog_spec_char_alert").dialog("open");
        return;
      }
    }
    var MOLdata = Jmol.scriptWaitOutput(jmolApplet, "write PDB");
    if (MOLdata.indexOf("HETATM") == -1){
      $("#dialog_no_ligand").dialog("open");
      return;
    }else{
      var regex = new RegExp("\n", "gi");
      document.getElementById("molinfo").value = MOLdata.replace(regex,"\\n"); 
    }
  }
  document.getElementById("ligand_form").submit(); //This sends it and hands it off to django/Python for processing
}


