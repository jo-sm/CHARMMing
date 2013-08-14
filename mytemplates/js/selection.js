//Color constants for JSmol
var cyan = '[x00ffff]';
var magenta = '[xc800c8]';
var green = '[x008000]';

function writeLinkAtomLines(seg_ids,modelType,num_linkatoms,current_div_id)
{
 var seg_list = new Array();
 seg_list = seg_ids.split(' ');
 //the last index of the array will be a blank space, so splice it!
 seg_list.splice(seg_list.length-1,1);
 num_patches = num_linkatoms;
 var model_string = "";
 var layer_num = "";
 var real_div_id = "linkatom_lines"
 if (modelType == "oniom"){
  layer_num = current_div_id.split("_"); //e.g. layer_2 turns into [layer, 2]
  model_string = "_layer_" + layer_num[layer_num.length -1];
  real_div_id = real_div_id + model_string;
  }
 var text = '<table class="qmmm_table" style="margin-left:auto;margin-right:auto;">';
 optionvalues = "";
 for(var b =0; b < seg_list.length; b++)
 {
  optionvalues = optionvalues + '<option value="' + seg_list[b] + '">' + seg_list[b] + '</option>';
 }
 for(var i = 0; i < num_patches; i++)
 {
  tempi = i+1;
  text = text + '<tr><td>QMHost '+tempi+'.</td><td> First SEGID:';
  text = text + '<select size="1" name="linkqmsegid' + model_string +"_"+ i +'">';
text = text + optionvalues + '</select><td> QM RESID:<input type="text" id="linkqm'+model_string+"_"+i+'" name="linkqm'+model_string+"_"+i+'" size=4></td> <td> QM Atom Type: <input type="text" id="qmatomtype'+model_string+"_"+i+'" name="qmatomtype' +model_string+"_"+i + '" size=5> </td></tr><tr><td>MMHost '+tempi+'.</td><td> MM SEGID: <select size="1" name="linkmmsegid'+model_string+"_"+i +'">' + optionvalues + '</select> </td><td>MM RESID:<input type="text" id="linkmm'+model_string+"_"+i+'" name="linkmm'+model_string+"_"+i+'" size=4></td><td> MM Atom Type: <input type="text" id="mmatomtype'+model_string+"_"+i+'" name="mmatomtype'+model_string+"_"+i + '" size=5> </td></tr>';
 }
  text= text + '</table>';
  document.getElementById(real_div_id).innerHTML = text;
}

function showHideQM(){
    if($("#useqmmm").is(":checked")){
      //      $(".qmmm_params").show();
      $(".model_selection").show();
      if ($("#usepbc").length > 0){
        $("#usepbc").attr("disabled",true);
        $("#usepbc").attr("checked",false); //just in case
      }
    }else{
    //      $(".qmmm_params").hide();
      $(".qmmm_params").remove(); //If this is a 0-length, then it doesn't matter, it won't except
      $(".oniom_params").remove();
      $(".model_selection").hide();
      if($("#usepbc").length > 0){
        $("#usepbc").attr("disabled",false);
      }
    }
  }

//For atom selection (QM/MM)
function goto_atomselect(){
  var inputs = document.getElementsByName("ptask");
  var form = null;
  if (inputs.length > 0){
    for(i=0;i<inputs.length;i++){
      if(inputs[i].checked){
        document.getElementById("task_id").value = inputs[i].value;
        break;
      }
    }
    var action = document.URL.split("/");
    var source = action[action.length -2];
    document.getElementById("source").value = source;
    if(source == "energy"){
      form = document.getElementById("ener_form");
    }
    if(source == "minimize"){
      form = document.getElementById("min_form");
    }if(source == "normalmodes"){
      form = document.getElementById("nma_form");
    }
//Note: Update here to create more QM/MM boxen in other pages
    form.action="/charmming/selection/";
    form.onsubmit= function(event){return true;};
    form.submit();
  }else{
    $("#dialog_coords_alert").dialog("open");
//    alert("No coordinates present. Please run at least one calculation on the full atom set before performing any QM/MM operations.")
  }
}

function hideQMBoxes(qmbox_number){
  var number_to_check = 0;
  if (typeof qmbox_number == "number"){
    number_to_check = qmbox_number;
  }else{
    var foo = qmbox_number.split("_");
    number_to_check = parseInt(foo[foo.length - 1]);
  }
  //We can make this recursive, but it's not useful.
  //We use the "layers" global included in the qmmm_params page.
  var current_layer = number_to_check;
  while(current_layer <= layers){
    $("#qmmm_params_layer_"+current_layer.toString()).hide();
    $("#mm_box_layer_"+current_layer.toString() + " input").attr("checked",true);
    current_layer = current_layer + 1;
  }
  $("#highest_qm_layer").val(number_to_check - 1);
  return true;
}
    
  

function create_bynum(add, atominfo){ //Universal bynum creation to make things easier
  var atomsele = "bynum " + atominfo[0].atomno;
  var i=1;//Hack because forloops crash firefox via allocation overflow somehow
  var start_range = atominfo[0].atomno;
  var end_range = start_range; // Add 1 every time
  while(i <= atominfo.length){ //Atominfo is sorted! This is important.
    if (add && atominfo[i - 1].color == cyan){
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
        if (start_range != atominfo[0].atomno){ //Fixes double-bynum glitch
          atomsele = atomsele + " .or. bynum " + start_range;
        }
      }
      if(i < atominfo.length){
        start_range = atominfo[i].atomno;
        end_range = start_range;
      }
  }
  i = i + 1;
  }
  var model_string = "";
  if(modelType == "oniom"){ //Add code for more models here
    model_string = "_layer_"+current_layer;
  }
  if (!(add)){ //No selection
    $("#atomselection"+model_string).val(atomsele);
  }else{ //Add selection, you shouldn't select things that are already in there
    $("#atomselection"+model_string).val($("#atomselection"+model_string).val() + " .or. " + atomsele);
  }
}

function submit_atomselect(){
  if ($("#add").is(":visible")){ //Stop playing with our code!
      reset_select(true,false);
      return false;
    }
  var atominfo = Jmol.getPropertyAsArray(jmolApplet, "atominfo", "selected and displayed"); //the "displayed" bit is for layers.
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
  $(".addselectbutton").show();
//        document.getElementById("atomselect_form").submit();
}

function add_atomselect(){
  if ($(".atomselectbutton").is(":visible")){ //Hey! Stop it.
      reset_select(true,false);
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
  var carbregex = new RegExp("C[0-9]+","g");
  if (atomarray.length != 1){ //Oh dear. That's no good.
    $("#dialog_too_many_links").dialog("open");
    return false;
  }else{
    var atom = atomarray[0]; //There's only supposed to be one. Watch our for CALCIUM!
    var atomname = atomarray[0].info.split(".")[1].split(" ")[0]; //e.g. from [UNK]1.H6 #6 -> ([UNK]1),(H6 #6) -> (H6)
    var atomresno = atomarray[0].info.split("]")[1].split(".")[0]; //e.g. from [UNK]1.H6 #6 -> ([UNK]),(1.H6 #6) -> (1)
    if(atomname != null){
      if (atomname.toUpperCase() != "CA" && carbregex.exec(atomname.toUpperCase()) == null && atomname.toUpperCase() != "C" && atomname.toUpperCase() != "CT2" && atomname.toUpperCase() != "CB")
  /*        && atomname.toUpperCase() != "CG" && atomname.toUpperCase() != "CD" && atomname.toUpperCase() != "CD1" && atomname.toUpperCase() != "CD2"
          && atomname.toUpperCase() != "CE")*/{
        $("#dialog_bad_bond").dialog("open"); //At this point, you have an atom selected, so either you unselect it, or you reset the whole thing
        return false;
      }else{
       /* if (atom.resname.toUpperCase() == "PRO" && (atomname.toUpperCase() != "CA" && atomname.toUpperCase() != "CB") || 
            atom.resname.toUpperCase() == "TRP" && (atomname.toUpperCase() != "CA" && atomname.toUpperCase() != "CB" && atomname.toUpperCase() != "CG"){
              $("#dialog_bad_bond").dialog("open");
              return false;
          } */
      }

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
    var divid = "#" + divtype + divnum;
    var identifier = $(divid).val().split("\t"); //This string "identifies" the atom, thus identifier
    var selection_string = ""; //Holds the Jmol script
    //Need to add in code that checks chain_terminators for stuff
    if (!($(divid).val() == "" || $(divid).val() == null)) {
        var seg_loc = segmentlist.indexOf(identifier[2]);
        if (!(seg_loc >= 0)){
          $("#dialog_bad_linkatom_input").dialog("open");
          return false;
        }
        selection_string = "select " + identifier[0] + "." + identifier[1] +
          " and atomno >= " + chain_terminators[seg_loc]  + ((seg_loc + 1 < chain_terminators.length) ? (" and atomno < " + (chain_terminators[seg_loc + 1])):"");
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
      $("#mmbutton"+divnum).prop("disabled",false);
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
      var qminfo = $("#qmhost"+divnum).val().split("\t"); //Get the number which forms the first half of the qmhost input box
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
    $(divid).val(atomresno + "\t" + atomname + "\t" + atom_segid);
    Jmol.script(jmolApplet, "selectionhalos on");
  }
  return false;
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
  

function reset_select(default_color, wipe_all_below){
  //default_color means restore to original state. We DON'T want this for layered selections.
  if(default_color){
    Jmol.script(jmolApplet, "selectionhalos off;select all;color jmol;select none;selectionhalos on;");
  }
  $(".atomselectbutton").show(); //This way the user can't touch it till JSmol renders
  if(modelType == "oniom"){
    if (wipe_all_below != false){
      $("#atomselection_layer_"+wipe_all_below).val("");
      $("#linkatom_num_layer_"+wipe_all_below).val("");
      $("#linkatom_inputs_layer_"+wipe_all_below).html("");
    }else{
      var curr_modifier = total_layers - 1; //Wipe out all the fields associated to any selections
      while (curr_modifier >= current_layer){
        $("#atomselection_layer_"+curr_modifier).val("");
        $("#linkatom_num_layer_"+curr_modifier).val("");
        $("#linkatom_inputs_layer_"+curr_modifier).html("");
        curr_modifier = curr_modifer - 1;
      }
    }
  }else{
    $("#atomselection").val("");
    $("#linkatom_num").val("");
    $("#linkatom_inputs").html("");
  }
  $(".addselectbutton").hide();
  return "successful";
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
    console.log(qm_mm_selectors[i]);
    selector = selector.value.split("\t");
    if (qm_mm_selectors[i].value.length > 17 || numre.exec(selector[0]) == null || carbre.exec(selector[1]) == null || segre.exec(selector[2]) == null){
      $("#dialog_bad_linkatom_input").dialog("open");
      return false;
    }
    i = i + 1;
  }
  //So we just tested for whether the user is trying to destroy our site. We carry on, branching in half
  //We have input fields telling us what our model type is and how many layers we have, so when we do this we need to "respawn" - i.e., we need to make more boxen.
  //We also need a "layer counter" - i.e., how many layers we have gone through so far.
  if(modelType == "oniom" && total_layers > 2){
    //Verify that the current layer is not the maximum one, i.e. 1
    if(current_layer == 1){
      $("#atomselect_form").prop("onsubmit","");
      $("#atomselect_form").submit();
      return true;
    }else{
      var current_layer_input_html = $("#atom_selection_inputs_layer_"+current_layer).html();
      var modified_div = "atom_selection_inputs_layer_"+(current_layer - 1);
      var next_layer_div_string = "<div id='"+modified_div+"'></div>";
      $("#linkatom_inputs_layer_"+current_layer).after(next_layer_div_string);
      current_layer = current_layer - 1; //Update the currently active layer
      $("#"+modified_div).hide(); //SO that the removals don't get seen by the user
      $("#"+modified_div).html(current_layer_input_html);
      $("#"+modified_div+" .removable_copies").remove(); //Remove all excess fields...
      //THis field we have issues templating since it already has a class, therefore we update it by hardcoding.
      $("#" + modified_div+" #linkatom_num_layer_"+(current_layer+1)).prop("id","linkatom_num_layer_"+current_layer);
      $("#" + modified_div+" #linkatom_num_layer_"+(current_layer+1)).prop("name","linkatom_num_layer_"+current_layer);
      var divs_to_update = $("#" + modified_div + " .layer_update"); //Make all the layers change...
      $("#" + modified_div+" #linkatom_inputs_layer_"+(current_layer+1)).html(""); //Blank out the linkatom_inputs display
      var layer_regex = new RegExp("layer_[0-9]","g");
      var curr_div = 0;
      //There's probably a more jQuery way of doing this but I don't trust replacing en masse when each field is different
      while (curr_div < divs_to_update.length){
        if (divs_to_update[curr_div].tagName == "INPUT" || divs_to_update[curr_div].tagName == "BUTTON"){
          divs_to_update[curr_div].name = divs_to_update[curr_div].name.replace(layer_regex,"layer_"+current_layer);
        }
        divs_to_update[curr_div].id = divs_to_update[curr_div].id.replace(layer_regex,"layer_"+current_layer);
        if(divs_to_update[curr_div].tagName == "H2"){
          divs_to_update[curr_div].innerHTML = divs_to_update[curr_div].innerHTML.replace("Layer "+(current_layer + 1),"Layer "+current_layer);
        }
        curr_div = curr_div +1;
      }
      if(highest_qm_layer < current_layer){ //This is basically only here for 4-layer since it's otherwise quite impossible to get to this point.
        $("#header_layer_"+current_layer).remove();
        $("#linkatom_num_layer_"+current_layer).remove();
      }
      var current_layer_selectbuttons_html = $("#atomselect_buttons").html(); //Hack for Lee's button doubling.
      $("#atomselection_layer_"+current_layer).after(current_layer_selectbuttons_html);
      $("#"+modified_div+" .addselectbutton").remove(); //Extraneous buttons...
      $("#selectbutton_layer_"+(current_layer + 1)).remove()
      $("#"+modified_div).show();
      Jmol.script(jmolApplet, "select color='"+cyan+"' or color='"+green+"';display selected;color jmol;select none;selectionhalos on;center displayed;zoom 0;");
      //The above selects our QM/MM region, the link atoms, turns off selectionhalos
      //so that the user doesn't see JSmol changing things, displays only the selected
      //atoms, centers the display, zooms such that the displayed atoms fill the entire screen
      //selects nothing, then turns selectionhalos back on. Phew!
      reset_select(false,current_layer);
      return false;
    } 
  }
  $("#atomselect_form").prop("onsubmit","");
  $("#atomselect_form").submit();
}



function submit_linkatoms(){
  var model_string = ""
  if(total_layers > 1){
   model_string = "_layer_" + current_layer; //This is different from model_string below, note th 
  } //Add code for more models here. Using model_string allows us to use less if statements.
  if ($("#atomselection"+model_string).val() == ""){ //empty atom selection
    $("#dialog_noatoms_alert").dialog("open");
    return false;
  }
  if (isNumeric($("#linkatom_num"+model_string).val())){
    var boxenNumber = parseInt($("#linkatom_num"+model_string).val());
    i=0;
    $("#linkatom_inputs"+model_string).html("<hr><br>"); //blank it before something goes wrong
    model_html_string = model_string  + "_" //This adds a _ for clarity later
    //Let's use while loops again...I don't trust for anymore.
    while(i < boxenNumber){
      $("#linkatom_inputs"+model_string).append('QM host:<input type="text" readonly class="qmmminput" name="qmhost'+model_html_string + i + 
          '" id="qmhost'+model_html_string + i + '" value="">&nbsp;<button id="qmbutton' +model_html_string+ i + '" class="selectbutton" onclick="select_qmsele(\'qmhost\', \'' +model_html_string+ i + '\');return false;">Select</button>&nbsp;');
      $("#linkatom_inputs"+model_string).append('MM host:<input type="text" readonly class="qmmminput" name="mmhost' +model_html_string+ i + 
          '" id="mmhost'+model_html_string + i + '" value="">&nbsp;<button class="selectbutton" id="mmbutton'+model_html_string + i + '"  disabled onclick="select_qmsele(\'mmhost\',\'' +model_html_string+ i + '\');return false;">Select</button>&nbsp;<br>');
      i = i + 1;
    }
    $("#linkatom_inputs"+model_string).append('<br><button id="selectbutton_layer_' + current_layer + '" class="selectbutton layer_update" onclick ="submit_qm_mm_atoms();return false;">Submit Selection with Link Atoms</button>');
    $("#linkatom_inputs"+model_string).css("display", "inline");
  }else{
    $("#dialog_bad_linkatoms").dialog("open");
    return false;
  }
  }

//Incoming jQueryUI error messages


  $(function($){
    $("#dialog_bad_params").dialog({
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
        $( "#dialog_coords_alert").dialog({
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
      reset_select(true,false);
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
