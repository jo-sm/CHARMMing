var jmolCounter = -1;
var langevinCounter = -1;
var mdCounter = -1;
var selected_div = "";
// Tab variables...


//This is for badhet selection on the build structure page
function showUploads(segdiv){
    var segname = segdiv.split("_");
    segname = segname[segname.length - 1]; //Get the last part...
    if ($("#"+segdiv).val() == "upload"){
      $("#toppar_upload_"+segname).css("display","inline");
    }else{
      $("#toppar_upload_"+segname).css("display","none");
    }
}

// The following is taken from the select/edit structures page
function send_select_form(form) {
   $.ajax({
      url:"/charmming/select/" + form.id,
      type:"post",
    });
//   new Ajax.Request("/charmming/select/" + form.id, {method:'post', asynchronous:true});
}

function send_delete_form(filename) {
   $.ajax({
      url:"/charmming/deletefile/",
      type:"post",
      async:false,
      data:{'filename':filename}
    });
      
//   new Ajax.Request("/charmming/deletefile/", {method:'post', asynchronous:false, parameters: {'filename':filename}});
}

function send_ligand(){
  window.location = "/charmming/ligand_design/";
}

// The following is taken from fileuploadform...
var last_option_id ="";

function setOptionId(div)
{
 last_option_id = div.id;
}


function addResidue(residue)
{
 sequ_field = document.getElementById('id_sequ');
 sequ_text = document.getElementById('id_sequ').value;
 //If the line is not blank, add a space
 if(sequ_text.length > 1)
 {
  sequ_text += " ";
 }
 sequ_text += residue;
 sequ_field.value = sequ_text;
}

//This swaps two words separated by an underscore
//example: bob_pdb is returned as pdb_bob
//this is done because the form ID is the inverse
//of the div_id it controls
//upload_pdb controls pdb_upload
function getChoiceDiv(div_id)
{
 var re =  /[a-z]*/;
 var first_part = re.exec(div_id);
 var re2 = /_[a-z]*/;
 var second_part = re2.exec(div_id);
 second_part[0] = second_part[0].replace(/[_]/,'');
 new_div = second_part[0] + "_" + first_part[0];
 return new_div;
}

//This gets the lesson value then passes it to the parent
function getLessonValue()
{
  return document.getElementById("lesson").value;
}

// end of stuff copied over from fileuploadform

//displays the div_id sent to it, specially designed for
//fileuploadform page
function displayFields(div_id)
{
 if(div_id !='rtf_upload' && div_id != 'prm_upload')
 {
  var choices = document.getElementsByName('submit_structure_option');
  for(var i=0; i < choices.length; i++)
  {
   if(getChoiceDiv(choices[i].id) != div_id)
    choices[i].checked = false;
    setVisible(getChoiceDiv(choices[i].id),"none");
  }
 } 
 if(div_id == 'rtf_upload' || div_id == 'prm_upload')
 {
  if(isVisible(div_id))
  {
   setVisible(div_id,"none");
  }
  else
  {
   setVisible(div_id,"block");
  }
  checkId(last_option_id);
  setVisible(getChoiceDiv(last_option_id),"block");
  return;
 }
 if(isVisible(div_id))
  setVisible(div_id,"none");
 else if(div_id == 'crd_upload')
 {
  setVisible(div_id,"block");
  setVisible('psf_upload',"block");
 }
 else
  setVisible(div_id,"block");

}


function setCGDisplay(id)
{
  var i=0;
  var check_array = $("input[name=buildtype]");
  while (i < check_array.length){
    if(!(check_array[i].checked)){
      $("#"+check_array[i].value+"_display").hide();
    }else{
      $("#"+check_array[i].value+"_display").show();
    }
    i = i + 1;
  }
  //Go through all of them, makes the code nice and generic. Can't use enhanced for loop again...why?
}

function addResidue(residue)
{
 sequ_field = document.getElementById('id_sequ');
 sequ_text = document.getElementById('id_sequ').value;
 //If the line is not blank, add a space
 if(sequ_text.length > 1)
 {
  sequ_text += " ";
 }
 sequ_text += residue;
 sequ_field.value = sequ_text;
}

//This swaps two words separated by an underscore
//example: bob_pdb is returned as pdb_bob
//this is done because the form ID is the inverse
//of the div_id it controls
//upload_pdb controls pdb_upload
function getChoiceDiv(div_id)
{
 var re =  /[a-z]*/;
 var first_part = re.exec(div_id);
 var re2 = /_[a-z]*/;
 var second_part = re2.exec(div_id);
 second_part[0] = second_part[0].replace(/[_]/,'');
 new_div = second_part[0] + "_" + first_part[0];
 return new_div;
}

//Sends AJAX request to delete a job
function killJob(jobid)
{
  $.ajax({
      url:"/charmming/killjob/" + jobid + "/",
      type:"post",
      async:false
    })
// new Ajax.Request('/charmming/killjob/' + jobid + '/',{method:'post',asynchronous:false});
}

//Sends AJAX request to delete student
function deleteStudent(username)
{
  $.ajax({
      url:"/charmming/normalmodes/deletestudent" + username,
      type:"post",
      async:false
    });
// Ajax.Request('/charmming/teachers/deletestudent/' + username,{method:'post',asynchronous:false});
}

//This is used to display more options when a user wants to select
//solvation structure size
function setSolvSize()
{
 if($("#set_pref").is(":checked")){
  if($("#solv_struc").val() == 'cubic'){
    $("#set_x").show();
    $("#set_y").show();
    $("#set_z").show();
    return "cubic";
  }
  else if($("#solv_struc").val() == 'sphere' || $("#solv_struc").val() == 'rhdo'){
    $("#set_x").show();
    $("#set_y").show();
    $("#set_z").show();
    return "rhdo";
  }
  else if($("#solv_struc").val() == 'hexa' || $("#solv_struc").val() == 'tetra'){
    $("#set_x").show();
    $("#set_y").show();
    $("#set_z").hide();
    return "hexa";
  }
  $("#no_pref_input").hide();
  return "set_pref";
 }
}

//to make replica exchange visible
function checkReplicaExchange(div_id,change_div_id)
{
 if(document.getElementById(div_id).checked)
 {
  setVisible(change_div_id,"block");
 }
 else
 {
  setVisible(change_div_id,"none");
 }
}

//for use with viewing the output files
function checkDisplayFile(filename) 
{
  return 1;
}

function hideFileList(header)
{
  var filediv = document.getElementById(header);
  var fileheader = document.getElementById("collapse-" + header); 
  if (filediv.style.display == "block")
    {
      filediv.style.display = "none";
      fileheader.innerHTML = "<h2><button onclick=\"hideFileList('" + header + "');\">+</button>" + header + " files</h2>";
  }else{
      filediv.style.display = "block";
      fileheader.innerHTML = "<h2><button onclick=\"hideFileList('" + header + "');\">-</button>" + header + " files</h2>";
    }
}
      
//for use with Normal modes
//TODO: Replace this qwith showing QM/MM template!
function showHideQMMM()
{
 if(document.getElementById('usevibran').checked) {
    setVisible("qmmmform","none");
 } else {
    setVisible("qmmmform","block");
 }
}

function checkUncheckNMA(div_id,other_div,change)
{
 if(document.getElementById('usevibran').checked)
 {
  changeTabLabel(AjaxTabs.GetFocusedTabId(),'Normal Modes Analysis: Vibran');
  setVisible(div_id,change);
  setVisible(other_div,"none");
 }
 else if(document.getElementById('useenm').checked)
 {
  changeTabLabel(AjaxTabs.GetFocusedTabId(),'Normal modes Analysis: ENM');
  setVisible(div_id,change);
  setVisible(other_div,"none");
 }
}

//for use with md
function checkUncheck(div_id,other1_div,other2_div,other3_div,change)
{
  setVisible(div_id,change);
  setVisible(other1_div,"none");
  setVisible(other2_div,"none");
  setVisible(other3_div,"none");
}



function send_form_gener(form_name,url,div_to_change) 
{
  form = document.getElementById(form_name);
  $.ajax({
      url:url,
      type:"post",
      data:$(form).serialize(),
      success: function(responseData){
      div_to_change.innerHTML = responseData;}
    });
//  new Ajax.Updater(div_to_change,url,{method:'post', asynchronous:true, parameters:Form.serialize(form)});
}

function send_form_mdanal(form,link,divupdate)
{
  $.ajax({
      url:link,
      type:"post",
      data:$(form).serialize()
    });
//   new Ajax.Request(link, {method:'post', asynchronous:true, parameters:Form.serialize(form)});
   changeStatus("mdform","mdform","mdanal");
   return false;
}

function send_form_rms(form,divupdate,filenames)
{
  var nsegs = 0;
  for(i=0;i<filenames.length;i++) {
     try {
        if(document.getElementById('id_' + filenames[i]).checked) {
           nsegs++;
        }
     } catch(e) {
        continue;
     }
  }
  if(nsegs < 2) {
     Dialog.alert("<center>Please select at least two structures to compare RMSDs.</center>", {width:300, height:100, okLabel: "close"});
     return false;
  }
  divid = document.getElementById(divupdate);
  divid.innerHTML = 'Calculating RMSD Matrix...';
  $.ajax({
      url: "/charmming/analysis/rmsd/",
      type: "post",
      data:$(form).serialize(),
      success: function(returnData){
      divid.innerHTML = returnData;}
      });

}

function send_form_energy(form,divupdate,filenames)
{
    divid = document.getElementById(divupdate);
    divid.innerHTML = 'Calculating Energy...';
    $.ajax({
        url: "/charmming/energy/",
        type:"post",
        data:$(form).serialize(),
        success: function(returnData){
        divid.innerHTML = returnData;}
        });
//    new Ajax.Updater(divupdate,'/charmming/energy/', {method:'post', asynchronous:true, parameters:Form.serialize(form)});
}

function send_form_oxired(form,link,divchange,divupdate,filenames)
{

    var rdxchecked = 0;
    var rbut;  
 
    rbut = document.getElementsByName("picksite");
    for(i=0; i < rbut.length; i++) {
      if(rbut[i].checked) {
        rdxchecked = 1;
      }
    }

    if(rdxchecked == 0)
    {
     Dialog.alert("<center>Please select an oxidation/reduction site.</center>", {width:300, height:100, okLabel: "close"});
     return false;
    }
    $.ajax({
        url:link,
        type:"post",
        data:$(form).serialize()
      });
//    new Ajax.Request(link, {method:'post', asynchronous:true, parameters:Form.serialize(form)});
    changeStatus(divchange,divupdate,"oxired");
    return false
}

function send_form_nma(form,link,divchange,divupdate,filenames)
{

  var ifsegchecked = 0;

    if(document.getElementById('usevibran').checked || document.getElementById('useenm').checked) {
      $.ajax({
          url:link,
          type:"post",
          data:$(form).serialize()
        });
//      new Ajax.Request(link, {method:'post', asynchronous:true, parameters:Form.serialize(form)});
      changeStatus(divchange,divupdate,"nma");
      return false
    } else {
      Dialog.alert("<center>Please Select an Option</center>", {width:300, height:100, okLabel: "close"});
      return false;
    } 
}

//like writeDisuLines but for restraints
function writeRestraintLines(div_id)
{
 num_restraints = parseInt(document.getElementById('num_restraints').value);
 text = '<table class="qmmm_table" style="margin-left:auto;margin-right:auto">';
 for(var i=0;i<num_restraints;i++)
 {
  number = i+1;
  text = text + '<tr><td><p style="text-align:center">' + number + '.  cons harm bestfit mass force 100.0 select <input type="text" id="cons_selection' + i + '" name="cons_selection' + i +'"> end</p></td></tr>';
 }
 document.getElementById(div_id).innerHTML = text;
}

//pre: num_replica_div_id is the textbox where the user enters the number of replicas they want, edit_text_div_id is where the 
//     lines will get written out to
function writeRexTempLines(num_replica_div_id,edit_text_div_id,premadelines)
{
 num_replicas = parseInt(document.getElementById(num_replica_div_id).value);
 text = '<table cellspacing=0 cellpadding=0 border=0><tr><td>';
 for(i = 1; i < num_replicas + 1; i++)
 {
  text = text + '<tr><td>Temperature: ' + i + ' </td><td><input type="text" id="rextemp"' + i + '" name="rextemp' + i + '" size=4></td></tr>';
 }
 text = text + '</table>';
 document.getElementById(edit_text_div_id).innerHTML = text;
}

//when a user types in a number of disulfide bond patches on minimizeform.html
//it will print out the html for the number of patches
function writeDisuLines(seg_ids,num_disu_div_id,div_id,premadelines)
{
 var seg_list = new Array();
 seg_list = seg_ids.split(' ');

 //the last index of the array will be a blank space, so splice it!
 seg_list.splice(seg_list.length-1,1);
 num_patches = parseInt(document.getElementById(num_disu_div_id).value);
 text = '<table cellspacing=0 cellpadding=0 border=0><tr><td>';
 optionvalues = "";
 for(var b = 0; b < seg_list.length; b++)
 {
  optionvalues = optionvalues + '<option value="' + seg_list[b] + '">' + seg_list[b] + '</option>';
 }
 for(var i = 0; i < num_patches; i++)
 {
  tempi = i+1;
  text = text + '<tr><td>'+tempi+'.</td><td> Starting SEGID:';
  text = text + '<select size="1" name="disustartsegid' + i +'">';
text = text + optionvalues + '</select><td> Starting RESID:<input type="text" name="disustart'+i+'" size=4></td><td> Ending SEGID: <select size="1" name="disuendsegid' + i +'">' + optionvalues + '</select> </td><td>Ending RESID:<input type="text" name="disuend'+i+'" size=4></td></tr>';
 }
  text= text + '</td></tr></table>';
  document.getElementById(div_id).innerHTML = text;
}

// This code triggers the filling in of the disulfide list with values from the PDB
function geneDisulfideLines(dsbonds,seg_ids)
{
  var num_dsbonds = dsbonds.length / 4;
  var target_div = document.getElementById("disulfide_solv");

  // set up the segment list
  var seg_list = new Array();
  seg_list = seg_ids.split(' ');
  seg_list.splice(seg_list.length-1,1);

  text = '<table cellspacing=0 cellpadding=0 border=0>\n';
  for(var i = 0; i < num_dsbonds; i++) {
    tempi = i + 1;
    text = text + '<tr><td>' + tempi + '.</td><td> Starting SEGID: <select  size="1" name="disustartsegid' + i +'">';

    var ds_segid1 = dsbonds[i*4];
    var ds_resno1 = dsbonds[i*4+1];
    var ds_segid2 = dsbonds[i*4+2];
    var ds_resno2 = dsbonds[i*4+3];
    for(var j = 0; j < seg_list.length; j++) {
      if(seg_list[j] == ds_segid1) {
        text = text + '<option value="' + seg_list[j] + '" selected>' + seg_list[j] + '</option>\n';
      } else {
        text = text + '<option value="' + seg_list[j] + '">' + seg_list[j] + '</option>\n';
      }
    }
    text = text + '</select></td>\n';
    text = text + '<td> Starting RESID: <input type="text" name="disustart' + i + '" value="' + ds_resno1 + '" size=4></td>\n';
    text = text + '<td> Ending SEGID: <select size="1" name="disuendsegid' + i + '">';
    for(var j = 0; j < seg_list.length; j++) {
      if(seg_list[j] == ds_segid2) {
        text = text + '<option value="' + seg_list[j] + '" selected>' + seg_list[j] + '</option>\n';
      } else {
        text = text + '<option value="' + seg_list[j] + '">' + seg_list[j] + '</option>\n';
      }
    }
    text = text + '</select></td>\n';
    text = text + '<td> Ending RESID:  <input type="text" name="disuend' + i + '" value="' + ds_resno2 + '" size=4></td></tr>\n';
  } // end for i

  target_div.innerHTML = text;

  // set the number of patches correctly
  var mytxtname = makeUniqueDiv('num_patches');
  document.getElementById(mytxtname).value = "" + num_dsbonds;

  // make it show up for the user
  var mydivname = makeUniqueDiv('usepatch');
  document.getElementById(mydivname).checked = true;
  checkPatch(makeUniqueDiv('usepatch'),makeUniqueDiv('patch_area'));
}

//makes a unique div by appending the tabId
//This is used in the solvation/minimization forms
//because in order for the AJAX tabs to work, each div must have a unique ID
//For example: Two minimization tabs would not work properly without this function
//because the hideShowPatch() would see two divs with the same id and only change
//one of them
//This whole thing is not worth is now given that we don't and never will again use AJAX tabs.
function makeUniqueDiv(div_id)
{
 try
 {
  var doc =  document.getElementById(div_id);
  doc.id = doc.id + AjaxTabs.GetFocusedTabId();
 }
 catch (error)
 {
  try
  {
   var doc = document.getElementById(div_id + AjaxTabs.GetFocusedTabId());
  }
  catch (e2)
  {
   var doc = document.getElementById(div_id);
  }
 }
 return doc.id;
}

//hides/shows the hetatm_determination
function hetatmDeter(div_id,het_list)
{
 ifsegchecked = 0;
 for(i=0;i<het_list.length;i++)
 {
   if(document.getElementById('id_' + het_list[i])) {
     if(document.getElementById('id_' + het_list[i]).checked) {
       ifsegchecked = 1;
     }
   }
}
 text = document.getElementById(div_id) 
 if(ifsegchecked == 1)
  text.style.display = "block";
 else if(ifsegchecked == 0)
  text.style.display = "none";
}
//hides the patching feature if user selects certain pdb
function hideShowPatch(div_id,should_show)
{
 patch_area = document.getElementById(div_id);
 if(!should_show)
 {
  patch_area.style.display = "none";
 }
 else
  patch_area.style.display = "block";
}

//unchecks everything with the given name
//used if there are radios and checkboxes
function uncheckAllInList(filenames)
{
 filename_list = filenames.split(' ');
 for(var i = 0; i < filename_list.length; i++)
 {
   var temp_array = document.getElementsByName(filename_list[i]);
   for(var b = 0; b < temp_array.length; b++)
   {
    document.getElementsByName(filename_list[i])[b].checked = false;
   }
 }
}

//checks an element based on its ID
function check(div_id)
{
 document.getElementById(div_id).checked = true;
}
//unchecks an element based on its ID
function uncheck(div_id)
{
 document.getElementById(div_id).checked = false;
}

//used for when the Apply Restrains feature is selected
function checkRestraints(check_id,display_div)
{
 display = document.getElementById(display_div);
 if(document.getElementById(check_id).checked == true)
 {
  display.style.display = "block";
 }
 else
 {
  display.style.display = "none";
 }
}

//used for when the Apply Restrains feature is selected
function checkShake(check_id,display_div)
{
 display2 = document.getElementById(display_div);
 if(document.getElementById(check_id).checked == true)
 {
  display2.style.display = "block";
 }
 else
 {
  display2.style.display = "none";
 }
}

//used for the use_custom_shake area
function showCustomShake(shake_div,bool)
{
 if(bool)
  document.getElementById(shake_div).style.display = "block";
 else
  document.getElementById(shake_div).style.display = "none";
 
}

//used to display patching in minimizeform,solvation,md,ld
//can also be used in general cases
function checkPatch(usepatch_div_id,div_id)
{ 
 if(document.getElementById(usepatch_div_id).checked)
 {
  setVisible(div_id,"block");
 }
 else 
 {
  setVisible(div_id,"none");
 }
}

//returns whether the div element is visible or not
//can be replaced by jQuery's $("#elem_id").is(":visible")
function isVisible(div_id)
{
 if(document.getElementById(div_id).style.display == "none")
  return false;
 return true;
}


//Can change an element's visibility, for example (block or none)
//can be replaced by jQuery's $("#elem_id").show() or hide()
function setVisible(div_id,change)
{
 var style = document.getElementById(div_id).style; 
 style.display = change;
}

function report_error(form,divtoupdate)
{
 var ajax_div = document.getElementById(divtopdate);
 $.ajax({
    url:"/charmming/reporterror/",
    type:"post",
    data:$(form).serialize(),
    success: function(returnData){
    ajax_div.innerHTML = returnData;}
    });
//  new Ajax.Updater(divtoupdate,'/charmming/reporterror/',{method:'post', asynchronous:true, parameters:Form.serialize(form)});
 return false;
}
 

function send_form_changepassword(form,link,divupdate) 
{
  var ajax_div = document.getELementById(divtoupdate);
  $.ajax({
      url:link,
      type:"post",
      async:false,
      data:$(form).serialize(),
      success: function(returnData){
      ajax_div.innerHTML = returnData;}
      });
// new Ajax.Updater(divupdate,link, {method:'post', asynchronous:false, parameters:Form.serialize(form)});
   return false;
}

function send_form(form,link,divchange,divupdate,filenames) 
{
  $.ajax({
      url:link,
      type:"post",
      data:$(form).serialize()
      });
//  new Ajax.Request(link, {method:'post', asynchronous:true, parameters:Form.serialize(form)});
  changeStatus(divchange,divupdate,"prep");
  return false;
}

function changeStatus(divchange,divupdate,caller)
{
 divupdate = divupdate + AjaxTabs.GetFocusedTabId();
 if(caller != "mdanal") {
   var text = "<center><h2>Job Submitted! You can check the status in the sidebar.</h2></center><br /> <hr noshade width=80% size=3%><div id=" + divupdate + "></div>";
   document.getElementById(divchange).innerHTML=text;
 } else {
   divupdate = "mdform";
 }

 var wordsOfWisdom = "";
 if(caller == "md") {
   wordsOfWisdom = "<br /><br /><center><h2> Once Molecular Dynamics is finished, the output PDBs will be available from the view PDBs page.</h2></center> <br /> <hr noshade width=80% size=3%>";
 } else if(caller == "ld") {
   wordsOfWisdom = "<br /><br /><center><h2> Once Langevin Dynamics is finished, the output PDBs will be available from the view PDBs page.</h2></center> <br /> <hr noshade width=80% size=3%>";
 } else if(caller == "nma") {
   wordsOfWisdom = "<br /><br /><center><h2> Once the Normal Mode Analysis calculation is finished, the output PDBs will be available from the view PDBs page.</h2></center><br /> <hr noshade width=80% size=3%>";
 } else if(caller == "prep") {
   wordsOfWisdom = "<br /><br /><center><h2> Once your structure is prepared with minimization and solvation, you can run dynamics.</h2></center><br /> <hr noshade width=80% size=3%>";
 } else if(caller == "oxired") {
   wordsOfWisdom = "<br /><br /><center><h2> Oxidation/Reduction job submitted. Please check back on the redox page for results. </h2></center><br /> <hr noshade width=80% size=3%>";
 } else if(caller == "mdanal") {
   wordsOfWisdom = "<br /><br /><center><h2> The data file with the requested properties is available from the &quot;download files&quot; page</h2></center><br /> <hr noshade width=80% size=3%>";
 }
 document.getElementById(divupdate).innerHTML = wordsOfWisdom;
}

function checkId(div_id)
{
 document.getElementById(div_id).checked = true;
}

//Drug Design stuff
function send_delete_ligand_form(id) {
    var r=confirm("You are about to delete a ligand. This ligand will be removed from the all jobs, sets and the analysis!");
    if (r==true)
    {
	//alert ('ddeleting ligand'+id);
        $.ajax({
            url:"/charmming/dd_substrate/deleteligand",
            type:"post",
            data:$(form).serialize()
          });
//        new Ajax.Request("/charmming/dd_substrate/deleteligand/", {method:'post', asynchronous:false, parameters: {'id':id}});
    }
    else
    {
	//alert ('cancelled');
    }
    //alert ('ddeleting ligand'+id);
    //new Ajax.Request("/charmming/dd_substrate/deleteligand/", {method:'post', asynchronous:false, parameters: {'id':id}});
    
}

function doDDJobDetailsOnLoad(id){
    //alert ('jobdetailonload'+id);           
    //var code = "jobid={{job.id}}";
    currperiodicalupdater = $("jobinfo").PeriodicalUpdater(("/charmming/dd_infrastructure/viewjobinfo" + id), {method:'post', multiplier:5, minTimeout:100, maxTimeout:5000},
        function(remoteData){document.getElementById("jobinfo").innerHTML = remoteData;});
//    currperiodicalupdater = new Ajax.Updater('jobinfo', '/charmming/dd_infrastructure/viewjobinfo/' + id, {method: 'post'});
    //code="jobid={{job.owner_id}}";

    currperiodicalupdater2 = $("resultsdiv").PeriodicalUpdater("/charmming/dd_infrastructure/viewjobresults" + id, {method:'post', multiplier:5, minTimeout:100, maxTimeout:5000},function(remoteData){document.getElementById("resultsdiv").innerHTML = remoteData;});
//    currperiodicalupdater2 = new Ajax.PeriodicalUpdater('resultsdiv', '/charmming/dd_infrastructure/viewjobresults/'+id, {method: 'post', frequency: 5});
    $("#resultsdiv").hide();
        //alert ("blah");

}


function viewDDJobResults(job_id) {

    var params = 'jobid='+job_id;
    var link=$("#viewhideddjobresults"); //this is never used again?
    //alert (document.getElementById('ddjobresultsdiv').style.display);
    if ((document.getElementById('ddjobresultsdiv').style.display=='none')){
  $.ajax({
      url:"/charmming/dd_infrastructure/viewjobresults/"+job_id,
      type:"post",
      success:function(responseData){
      document.getElementById('ddjobresultsdiv').innerHTML=responseData;}
    })
      $("#ddjobresultsdiv").show();
        //var style = document.getElementById('ddjobresultsdiv').style;
        //style.display = change;
//	new Ajax.Updater('ddjobresultsdiv',"/charmming/dd_infrastructure/viewjobresults/"+job_id, {method:'post', asynchronous:true});
        //link.innerHTML="Hide Results"
    }
    else{
      $("#ddjobresultsdiv").hide();
	//link.innerHTML="View Results"
        //var style = document.getElementById('ddjobresultsdiv').style;
	//style.display = change;
	
    }
       
}

function dd_send_project_select_form(form) {
  $.ajax({
      url:"/charmming/dd_infrastructure/select_project/"+form.id,
      type:"post"
    });
//    new Ajax.Request("/charmming/dd_infrastructure/select_project/"+form.id, {method:'post', asynchronous:true});
}

function getCheckedAvailableConformationsIds()
{
    var nodes=document.getElementsByName("availableconformationsdiv")[0].childNodes;
    var checkedvalues="";
    for(i = 0;i < nodes.length;++i)
	if(nodes[i].id=="availableconformationsinsidediv"){
	    var checkboxes = nodes[i].getElementsByTagName('input');
	    //alert(checkboxes.length);
	    for(j = 0;j < checkboxes.length;++j){
		    //alert(checkboxes[j].value + ' ' + checkboxes[j].checked);
		    //alert("j " + j);
		if (checkboxes[j].checked){
		    if (checkedvalues=='') {
			checkedvalues=checkboxes[j].value;
			           
		    }
		    else {
			checkedvalues=checkedvalues + ',' + checkboxes[j].value;
		    }
		    //alert(checkedvalues);
		}
	    }
	}
      //alert(checkedvalues);
    return checkedvalues;
}

function getCheckedProjectConformationsIds()
{
    var nodes=document.getElementsByName("projectconformationsdiv")[0].childNodes;
    var checkedvalues="";
    for(i = 0;i < nodes.length;++i)
	if(nodes[i].id=="projectconformationsinsidediv"){
	    var checkboxes = nodes[i].getElementsByTagName('input');
	    //alert(checkboxes.length);
	    for(j = 0;j < checkboxes.length;++j){
		    //alert(checkboxes[j].value + ' ' + checkboxes[j].checked);
		    //alert("j " + j);
		if (checkboxes[j].checked){
		    if (checkedvalues=='') {
			checkedvalues=checkboxes[j].value;
			           
		    }
		    else {
			checkedvalues=checkedvalues + ',' + checkboxes[j].value;
		    }
		    //alert(checkedvalues);
		}
	    }
	}
    return checkedvalues;
}

function AddConformations(project_id)
{
    RefreshProjectConformations(project_id,getCheckedAvailableConformationsIds(),'');
    //var code = 'ptid=' + $("protein_types").getValue();
    //var myAjax = new Ajax.Updater('availableconformationsdiv', '/charmming/updatemanagesetsgrids/', {method: 'post', parameters: code})
}

function RemoveConformations(project_id)
{
  //alert (getCheckedProjectConformationsIds());
    RefreshProjectConformations(project_id,'',getCheckedProjectConformationsIds());
    //var code = 'ptid=' + $("protein_types").getValue();
    //var myAjax = new Ajax.Updater('availableconformationsdiv', '/charmming/updatemanagesetsgrids/', {method: 'post', parameters: code})
}
  
function RefreshProjectInfo(project_id)
{
    var code = '&name=' + document.getElementById("project_name").value + '&description=' + document.getElementById("project_description").value;
    var myAjax = $.ajax({
        url:"/charmming/dd_infrastructure/updateprojectinfo/"+project_id+"/update/",
        type:"post",
        data:code,
        success: function(responseData){
        document.getElementById("projectinfodiv").innerHTML=responseData;}
      });
//    var myAjax = new Ajax.Updater('projectinfodiv', '/charmming/dd_infrastructure/updateprojectinfo/'+project_id+'/update/', {method: 'post', parameters: code})
//    I don't know why Yuri makes it a var. In jQuery it shouldn't even be relevant...
}
    
function DisplayProjectInfo(project_id)
{
    //alert (document.getElementById("gridset_name").value);
    //var code = codestring;
    //alert (code);
    var myAjax = $.ajax({
        url:"/charmming/dd_infrastructure/updateprojectinfo/"+project_id+"/refresh/",
        type:"post",
        success: function(responseData){
        document.getElementById("projectinfodiv").innerHTML=responseData;}
      });
//    var myAjax = new Ajax.Updater('projectinfodiv', '/charmming/dd_infrastructure/updateprojectinfo/'+project_id+'/refresh/', {method: 'post'})
    //Shouldn't this be in RefreshProjectInfo?
}

function RefreshAvailableConformations()
{
    //alert ($("targets").getValue());
    //var code = "protein_id=" + $("targets").getValue() + "&setid=" + $("targets").getValue();
    var myAjax = $.ajax({
        url:"/charmming/dd_infrastructure/projectavailableconformations/"+document.getElementById("targets").value+"/",
        type:"post",
        success: function(requestData){
        document.getElementById("availableconformationsdiv").innerHTML = requestData;}
        });
//    var myAjax = new Ajax.Updater('availableconformationsdiv', '/charmming/dd_infrastructure/projectavailableconformations/'+document.getElementById("targets").value+'/', {method: 'post'});
    //var myAjax = new Ajax.Updater('availableconformationsdiv', 'www.google.com', {method: 'post'})
}


function RefreshProjectConformations(project_id,addedids,removedids)
{
    //alert (document.getElementById("Receptor").value);
    //var code = "project_id={{ project_id }}&addedids=" + addedids + "&removedids=" + removedids;
        //var code = 'ptid=' + $("protein_types").getValue();
    //alert('/charmming/dd_infrastructure/projectconformations/{{project_id}}/'+addedids+'/'+removedids);
    var myAjax = $.ajax({
        url:"/charmming/dd_infrastructure/projectconformations/"+project_id+"/"+addedids+"/"+removedids+"/",
        type:"post",
        success: function(requestData){
        document.getElementById("projectconformationsdiv").innerHTML = requestData;}
        });
//    var myAjax = new Ajax.Updater('projectconformationsdiv', '/charmming/dd_infrastructure/projectconformations/'+project_id+'/'+addedids+'/'+removedids+'/', {method: 'post'})

    //alert(document.getElementById('conformations').style)
}

function CheckUncheck()
{

    //alert(document.forms['projectdetail'].elements['conformation_id']);
    if ( document.forms['projectdetail'].elements['conformation_id'].length )
    {
	for (var x = 0; x < document.forms['projectdetail'].elements['gridfileid'].length; x++)
	{
	    if (document.forms['projectdetail'].elements['checkall'].checked)
	    {
		document.forms['projectdetail'].elements['gridfileid'][x].checked = true;   
	    }
	     else
	    {
		document.forms['projectdetail'].elements['gridfileid'][x].checked = false;
	    }
	     
	}
    }
       else
    {
	//alert('checkall'+document.forms['projectdetail'].elements['checkall']);
        if (document.forms['projectdetail'].elements['checkall'].checked)
	{
	    document.forms['projectdetail'].elements['conformation_id'].checked = true;            
	}
	      else
        // Now jmolApplet can be referenced from all the other stuff on this page.
	{
	    document.forms['projectdetail'].elements['conformation_id'].checked = false;
	}
    }
}

function SetCheckUncheck()
{

    if ( document.forms['projectdetail'].elements['project_conformation_id'].length )
    {
	for (var x = 0; x < document.forms['projectdetail'].elements['project_conformation_id'].length; x++)
	{
	    if (document.forms['projectdetail'].elements['setcheckall'].checked)
	    {
		document.forms['projectdetail'].elements['project_conformation_id'][x].checked = true;   
	    }
	     else
	    {
		document.forms['projectdetail'].elements['project_conformation_id'][x].checked = false;
	    }
	     
	}
    }
       else
    {
	if (document.forms['projectdetail'].elements['setcheckall'].checked)
	{
	    document.forms['projectdetail'].elements['project_conformation_id'].checked = true;            
	}
	      else
	{
	    document.forms['projectdetail'].elements['project_conformation_id'].checked = false;
	}
    }
}
function CheckRedirect()
{
    alert (document.getElementById("gridset_name").value);
}


function RefreshProjectInfo(projectid,action)
{
    //alert (document.getElementById("gridset_name").value);
    var code = 'name=' + document.getElementById("project_name").value + '&description=' + document.getElementById("project_description").value;
    //alert (code);
    var myAjax = $.ajax({
        url:"/charmming/dd_infrastructure/updateprojectinfo/"+projectid+"/"+action+"/",
        type:"post",
        data:code,
        success: function(requestData){
        document.getElementById("projectinfodiv").innerHTML = requestData;}
        });
//    var myAjax = new Ajax.Updater('projectinfodiv', '/charmming/dd_infrastructure/updateprojectinfo/'+projectid+'/'+action+'/', {method: 'post', parameters: code})
}
        
function DisplayBlankForm()
{
    //alert (document.getElementById("gridset_name").value);
    //var code = codestring;
    //alert (code);
    var myAjax = $.ajax({
        url:"/charmming/dd_infrastructure/updateprojectinfo/0/refresh/",
        type:"post",
        success: function(requestData){
        document.getElementById("availableconformationsdiv").innerHTML = requestData;}
        });
//    var myAjax = new Ajax.Updater('projectinfodiv', '/charmming/dd_infrastructure/updateprojectinfo/0/refresh/', {method: 'post'})
}



/////////////////////////////dd ligand set details
function send_delete_ligandset_form(setid) {
    //alert('deleting set'+setid);   
    if(setid == 'all_sets'){
      $.ajax({
          url:"/charmming/dd_substrate/deleteligandset/",
          type:"post", 
          async:false
        });
//	new Ajax.Request("/charmming/dd_substrate/deleteligandset/", {method:'post', asynchronous:false });
  } else
    {
	//alert ('deleting set: '+setid);
  $.ajax({
      url:"/charmming/dd_substrate/deleteligandset/"+setid,
      type:"post",
      async:false
    })
//	new Ajax.Request("/charmming/dd_substrate/deleteligandset/"+setid, {method:'post', asynchronous:false});
    }
    //var viewgridsetsframe = parent.document.getElementById('viewgridsetscontainer');
    //viewgridsetsframe.src = viewgridsetsframe.src;
       //parent.refreshBottomdock();
}

function getCheckedAvailableLigandIds()
{
    var checkedvalues= new Array();
    $("#availableligandsinsidediv input:checked").each(function() {
        checkedvalues.push($(this).val());
        });
    return checkedvalues;
}

function getCheckedSetLigandIds()
{
    var checkedvalues = new Array();
    $("#setligandsinsidediv input:checked").each(function() {
        checkedvalues.push($(this).val());
      }); //This goes through every checked value and puts it into the array.
      //No more heavy forlooping and confusing structure required.
    return checkedvalues;
}

function AddLigands(set_id)
{
    RefreshSetLigands(set_id,getCheckedAvailableLigandIds(),'');
    //var code = 'ptid=' + $("protein_types").getValue();
    //var myAjax = new Ajax.Updater('availableconformationsdiv', '/charmming/updatemanagesetsgrids/', {method: 'post', parameters: code})
}

function RemoveLigands(set_id)
{
  //alert (getCheckedProjectConformationsIds());
    RefreshSetLigands(set_id,'',getCheckedSetLigandIds());
    //var code = 'ptid=' + $("protein_types").getValue();
    //var myAjax = new Ajax.Updater('availableconformationsdiv', '/charmming/updatemanagesetsgrids/', {method: 'post', parameters: code})
}
  
function RefreshSetInfo(set_id)
{
    var code = '&name=' + document.getElementById("set_name").value + '&description=' + document.getElementById("set_description").value;
    var myAjax = $.ajax({
        url:"/charmming/dd_substrate/updateligandsetinfo/"+set_id+"/update/",
        type:"post",
        data:code,
        success: function(requestData){
        document.getElementById("setinfodiv").innerHTML = requestData;}
        });
//    var myAjax = new Ajax.Updater('setinfodiv', '/charmming/dd_substrate/updateligandsetinfo/'+set_id+'/update/', {method: 'post', parameters: code})
}
    
function DisplaySetInfo(set_id)
{
    //alert (document.getElementById("gridset_name").value);
    //var code = codestring;
    //alert (code);
    var myAjax = $.ajax({
        url:"/charmming/dd_substrate/updateligandsetinfo/"+set_id+"/refresh/",
        type:"post",
        success: function(requestData){
        document.getElementById("ligandsetinfodiv").innerHTML = requestData;}
        });
//    var myAjax = new Ajax.Updater('ligandsetinfodiv', '/charmming/dd_substrate/updateligandsetinfo/'+set_id+'/refresh/', {method: 'post'})
}

function RefreshAvailableLigands()
{
    //alert ($("targets").getValue());
    //var code = "protein_id=" + $("targets").getValue() + "&setid=" + $("targets").getValue();
    //alert ('/charmming/dd_substrate/setavailableligands/'+document.getElementById("avaiableligandsets").value+'');
    var myAjax = $.ajax({
        url:"/charmming/dd_substrate/setavailableligands/"+$("#availableligandsets").val()+'',
        type:"post",
        success: function(requestData){
        document.getElementById("availableligandsdiv").innerHTML = requestData;}
        });
//    var myAjax = new Ajax.Updater('availableligandsdiv', '/charmming/dd_substrate/setavailableligands/'+document.getElementById("avaiableligandsets").value+'', {method: 'post'})
    //var myAjax = new Ajax.Updater('availableconformationsdiv', 'www.google.com', {method: 'post'})
}


function RefreshSetLigands(set_id,addedids,removedids)
{
    //alert (document.getElementById("Receptor").value);
    //var code = "project_id={{ project_id }}&addedids=" + addedids + "&removedids=" + removedids;
        //var code = 'ptid=' + $("protein_types").getValue();
    //alert('/charmming/dd_infrastructure/projectconformations/{{project_id}}/'+addedids+'/'+removedids);
    var myAjax = $.ajax({
        url:"/charmming/dd_substrate/setligands/"+set_id+"/"+addedids+"/"+removedids+"/",
        type:"post",
        success: function(requestData){
        document.getElementById("setligandsdiv").innerHTML = requestData;}
        });
    //var myAjax = new Ajax.Updater('setligandsdiv', '/charmming/dd_substrate/setligands/'+set_id+'/'+addedids+'/'+removedids+'/', {method: 'post'})

    //alert(document.getElementById('conformations').style)
}

function LigandSet_CheckUncheck()
{

    if ( document.forms['setdetail'].elements['ligandid'].length )
    {
	for (var x = 0; x < document.forms['setdetail'].elements['ligandid'].length; x++)
	{
	    if (document.forms['setdetail'].elements['checkall'].checked)
	    {
		document.forms['setdetail'].elements['ligandid'][x].checked = true;   
	    }
	     else
	    {
		document.forms['setdetail'].elements['ligandid'][x].checked = false;
	    }
	     
	}
    }
       else
    {
	if (document.forms['setdetail'].elements['checkall'].checked)
	{
	    document.forms['setdetail'].elements['ligandid'].checked = true;            
	}
	      else
	{
	    document.forms['setdetail'].elements['ligandid'].checked = false;
	}
    }
}

function LigandSet_SetCheckUncheck()
{

    if ( document.forms['setdetail'].elements['set_ligand_id'].length )
    {
	for (var x = 0; x < document.forms['setdetail'].elements['set_ligand_id'].length; x++)
	{
	    if (document.forms['setdetail'].elements['setcheckall'].checked)
	    {
		document.forms['setdetail'].elements['set_ligand_id'][x].checked = true;   
	    }
	     else
	    {
		document.forms['setdetail'].elements['set_ligand_id'][x].checked = false;
	    }
	     
	}
    }
       else
    {
	if (document.forms['setdetail'].elements['setcheckall'].checked)
	{
	    document.forms['setdetail'].elements['set_ligand_id'].checked = true;            
	}
	      else
	{
	    document.forms['setdetail'].elements['set_ligand_id'].checked = false;
	}
    }
}

function RefreshLigandSetInfo(ligandsetid,action,form)
{
    //alert (document.getElementById("gridset_name").value);
    //var code = 'name=' + document.getElementById("ligandset_name").value + '&description=' + document.getElementById("ligandset_description").value;
//    var code = 'name=' + $("#ligandset_name").value + '&description=' + $("i#ligandset_description").value;
    //TODO: Will this have to change for jQuery...?
    //alert (code);
    $.ajax({
        url:'/charmming/dd_substrate/updateligandsetinfo/'+ligandsetid+'/'+action+'/',
        type:"post",
        data:form, //Pre-serialized. Use serializeArray() to force POST
        success: function(requestData){
        document.getElementById("ligandsetinfodiv").innerHTML = requestData;}
        });
//    var myAjax = new Ajax.Updater('ligandsetinfodiv', '/charmming/dd_substrate/updateligandsetinfo/'+ligandsetid+'/'+action+'/', {method: 'post', parameters: code})
}
        
function DisplayBlankLigandForm()
{ //Can't you just do this by calling RefreshLigandSetInfo with 0, refresh?
    //alert (document.getElementById("gridset_name").value);
    //var code = codestring;
    //alert (code);
    var myAjax = $.ajax({
        url:'/charmming/dd_substrate/updateligandsetinfo/0/refresh/',
        type:"post",
        success: function(requestData){
        document.getElementById("ligandsetinfodiv").innerHTML = requestData;}
        });
//    var myAjax = new Ajax.Updater('ligandsetinfodiv', '/charmming/dd_substrate/updateligandsetinfo/0/refresh/', {method: 'post'})
}




/////dsf form
function RefreshDSFLigands()
{
    //alert (document.getElementById("Receptor").value);
    //var code = "ptid=" + $("protein_types").getValue() + "&setid=" + $("gridsets").getValue();
    //alert (code);
//good    var myAjax = new Ajax.Updater('ligandsdiv', '/charmming/dd_infrastructure/updatedsfligands/'+$("ligandsets").getValue()+'', {method: 'post'})
    var myAjax = $.ajax({
        url:'/charmming/dd_infrastructure/updatedsfligands/'+$("#ligandsets").val()+'',
        type:"post",
        success: function(requestData){
        document.getElementById("ligandsdiv").innerHTML = requestData;}
        });
//    var myAjax = new Ajax.Updater('ligandsdiv', '/charmming/dd_infrastructure/updatedsfligands/'+$("ligandsets").getValue()+'', {method: 'post'});
    //alert(document.getElementById('conformations').style)
}

function doDDJobDetailsOnLoad(job_id){
               
    //var code = "jobid={{job.id}}";
    currperiodicalupdater = $.ajax({
        url:'/charmming/dd_infrastructure/viewjobinfo/'+job_id,
        type:"post",
        success: function(requestData){
        document.getElementById("jobinfo").innerHTML = requestData;}
        });
//    currperiodicalupdater = new Ajax.Updater('jobinfo', '/charmming/dd_infrastructure/viewjobinfo/'+job_id, {method: 'post'});
    //code="jobid={{job.owner_id}}";
    currperiodicalupdater2 = $("ddjobresultsdiv").PeriodicalUpdater('/charmming/dd_infrastructure/viewjobresults/'+job_id, {method: 'post', multiplier:5, minTimeout:100, maxTimeout:5000}, function(remoteData){document.getElementById("ddjobresultsdiv").innerHTML = remoteData;});
//    currperiodicalupdater2 = new Ajax.PeriodicalUpdater('ddjobresultsdiv', '/charmming/dd_infrastructure/viewjobresults/'+job_id, {method: 'post', frequency: 5});
    $("#ddjobresultsdiv").hide();

        //alert ("blah");

}
        
function getCheckedDSFLigands()
{
//    var nodes=document.getElementsByName("ligandsinsidediv")[0].childNodes; 
    var nodes = $("#ligandsinsidediv")[0].childNodes;
    var checkedvalues="";
    for(i = 0;i < nodes.length;++i)
	if(nodes[i].id=="availableconformationsinsidediv"){
	    var checkboxes = nodes[i].getElementsByTagName('input');
	        //alert(checkboxes.length);
	    for(j = 0;j < checkboxes.length;++j){
		    //alert(checkboxes[j].value + ' ' + checkboxes[j].checked);
		    //alert("j " + j);
		if (checkboxes[j].checked){
		    if (checkedvalues=='') {
			checkedvalues=checkboxes[j].value;
			           
		    }
		    else {
			checkedvalues=checkedvalues + ',' + checkboxes[j].value;
		    }
		        //alert(checkedvalues);
		}
	    }
	}
      //alert(checkedvalues);
    return checkedvalues;
}


function DSFLigands_CheckUncheck()
{

    if ( document.forms['dsfform'].elements['id_ligand_file'].length )
    {
	for (var x = 0; x < document.forms['dsfform'].elements['id_ligand_file'].length; x++)
	{
	    if (document.forms['dsfform'].elements['dsfligandscheckall'].checked)
	    {
		document.forms['dsfform'].elements['id_ligand_file'][x].checked = true;   
	    }
	         else
	    {
		document.forms['dsfform'].elements['id_ligand_file'][x].checked = false;
	    }
	         
	}
    }
       else
    {
	if (document.forms['dsfform'].elements['dsfligandscheckall'].checked)
	{
	    document.forms['dsfform'].elements['id_ligand_file'].checked = true;            
	}
	      else
	{
	    document.forms['dsfform'].elements['id_ligand_file'].checked = false;
	}
    }
}
