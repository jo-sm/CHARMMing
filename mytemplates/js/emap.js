function send_emap_form(form,step)
{
   new Ajax.Request("/emap/" + step + "/", {method:'post', asynchronous:false,  parameters:Form.serialize(form)});
   if(step==3)
   {
    document.getElementById('emap').innerHTML = "Done";    
    document.getElementById('download_emap').style.display = "block";    
   }
   else
   {
    emapScript(step+1);
   }
}
function genUniqueStrucList(form_id)
{
 num = form_id.value;
 text = "";
 for(var i=0; i <num; i++)
 {

   text = text + '<tr><td>' + (i+1) + '.<input type="text" name="p'+i+'"></td><td><input type="text" name="n'+i+'" size="4"></td></tr>';
 }
 document.getElementById('unique_struc').innerHTML = text;
}

function emapScript(step)
{
 if(step == 0)
 {
  text = '<form method="post"  action"." onsubmit="send_emap_form(this,0);">';
  text = text + 'Step one: Type in the MAP_FILENAME in UPPERCASE letters (.CCP4 or .MRC):<br>';
  text = text + '<input type="text" name="em" value="FILENAME.CCP4"><br><br>';
  text = text + 'Next: How many unique structure files will there be?  ';
  text = text + '<input type="text" name="numofuniquestruc" value="0" size="4" onKeyUp="genUniqueStrucList(this);">';
  text = text + '</table><br>';
  text = text + '<table cellspacing=0 cellpadding=0>';
 text =  text +"<tr><td><center>Filename</center></td><td><center>Number of Times Structure Appears </center></td></tr>";
  text = text + '<div id="unique_struc"> </div></table>';
  text = text + '<input type="submit" value="Next Step">';
  text = text + '</form>';
  document.getElementById('emap').innerHTML = text;
 }
 else if(step == 1)
 {
  text = '<form method="post"  action"." onsubmit="return send_emap_form(this,1);">';
  text = text + 'Step Two: Choose the correlation function type: <br>';
  text = text + '<table><tr><td><input type="radio" name="f" value="1"></td><td>Density Correlation (DC)</td></tr>';  text = text + '<tr><td><input type="radio" name="f" value="2" CHECKED></td><td>Core-Weighted Density Correlation (CWDC)</td></tr>';
  text = text + '<tr><td><input type="radio" name="f" value="3"></td><td>Laplacian Correlation (LC)</td></tr>';
  text = text + '<tr><td><input type="radio" name="f" value="4"></td><td>Core-Weighted Laplacian Correlation (CWLC)</td></tr>';
  text = text + '</table><br><input type="submit" value="Next Step">';
  document.getElementById('emap').innerHTML = text;
 }
 else if(step == 2)
 {
  text = '<form method="post"  action"." onsubmit="return send_emap_form(this,2);">';
  text = text + 'Step three: Fill out the information on the map and grid:<br>';
  text = text + '<table>';
  text = text + '<tr><td>Resolution (Å):</td><td><input type="text" name="r" value="15"></td></tr>';
  text = text + '<tr><td>Grid interval (Can be left blank):</td><td><input type="text" name="d"></td></tr>';
  text = text + '<tr><td>Density cutoff:</td><td><input type="text" name="cut" value="0.001"></td></tr>';
  text = text + '<tr><td>Correlation tolerance:</td><td><input type="text" name="c" value="0.01"></td></tr>';
  text = text + '<tr><td>Translation grid (Å):</td><td><input type="text" name="gt" value="2"></td></tr>';
  text = text + '<tr><td>Rotation grid (Å):</td><td><input type="text" name="gr" value="2"></td></tr>';

  text = text + '</table><br><input type="submit" value="Next Step">';
  document.getElementById('emap').innerHTML = text;

 }
 else if(step == 3)
 {
  text = '<form method="post"  action"." onsubmit="return send_emap_form(this,3);">';
  text = text + 'Step three: Fill out the Monte Carlo Calculation information:<br>';
  text = text + '<table>';
  text = text + '<tr><td>Monte Carlo cycles:</td><td><input type="text" name="nc" value="10"></td></tr>';
  text = text + '<tr><td>Monte Carlo steps:</td><td><input type="text" name="ns" value="100"></td></tr>';
  text = text + '<tr><td>Monte Carlo max translation (Å):</td><td><input type="text" name="u" value="20.0"></td></tr>';   text = text + '<tr><td>Monte Carlo max rotation (Å):</td><td><input type="text" name="v" value="10.0"></td></tr>';  text = text + '<tr><td>Monte Carlo max temperature factor:</td><td><input type="text" name="t" value="0.01"></td></tr>';

  text = text + '</table><br><input type="submit" value="Done">';
  document.getElementById('emap').innerHTML = text;
 }

}

