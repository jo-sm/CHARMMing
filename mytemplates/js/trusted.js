function sendScriptChanges(textarea)
{
 filename = textarea.id;
 text = textarea.value;
 new Ajax.Request('/charmming/sendscriptchanges/', {method:'post', asynchronous:false, parameters:{'filename':filename,'text':text}});
}

function getInputScript(radio,div_id)
{
 filename = radio.value;
 new Ajax.Updater(div_id,'/charmming/getinputdata/', {method:'post', asynchronous:false, parameters:{'filename':filename}});
}

function send_form_editscript(form,link,div_id)
{
 jobtype = document.getElementById('jobtype')
 jobtype = document.getElementById('jobtype').value;
 if(jobtype != "energy")
 {
   new Ajax.Request(link, {method:'post', asynchronous:true, parameters:Form.serialize(form)});
 }
 else
 {
   document.getElementById(div_id).innerHTML = "Calculating...";
   new Ajax.Updater(div_id,link, {method:'post', asynchronous:true, parameters:Form.serialize(form)});
 }
 if(jobtype == 'minimization' || jobtype == 'solvation')
   changeStatus(div_id,'status','prep');
 else if(jobtype == 'nma')
   changeStatus(div_id,'status','nma');
 else if(jobtype == 'md')
   changeStatus(div_id,'status','md');
 else if(jobtype == 'ld' || jobtype == 'sgld')
   changeStatus(div_id,'status','ld');
 return false;

}
