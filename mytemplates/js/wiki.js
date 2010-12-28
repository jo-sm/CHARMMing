boardname = ""
function setBoardName(board_name)
{
 boardname = board_name;
}

function highlightStar(tdobj)
{
 var starvals = tdobj.id.split('_');
 var starnum = parseInt(starvals[1]);
 stars_table = document.getElementById(starvals[0] + '_stars');
 tds = stars_table.getElementsByTagName('td');
 for(i = 0;i < tds.length ; i++)
 {
   tds[i].style.backgroundImage="url(/charmming/images/wiki/smallstargrey.gif)";
 } 
 for(i = 1;i < starnum+1; i++)
 {
  document.getElementById(starvals[0] + '_' + i).style.backgroundImage="url(/charmming/images/wiki/smallstarred.gif)";
 }
 var topic_title = starvals[0];
 if(starnum == 1)
  document.getElementById(topic_title + "_currate").innerHTML = "Unhelpful";
 else if(starnum == 2)
  document.getElementById(topic_title + "_currate").innerHTML = "Somehwat helpful";
 else if(starnum == 3)
  document.getElementById(topic_title + "_currate").innerHTML = "Helpful";
 else if(starnum == 4)
  document.getElementById(topic_title + "_currate").innerHTML = "Very Helpful";
 else if(starnum == 5)
  document.getElementById(topic_title + "_currate").innerHTML = "Add to CHARMMING";
}

function resetStars(tableobj)
{
 tds = tableobj.getElementsByTagName('td');
 for(i = 0;i < tds.length; i++)
 {
  if(tds[i].className != 'check')
   tds[i].style.backgroundImage="url(/charmming/images/wiki/smallstargrey.gif)";
  else
   tds[i].style.backgroundImage="url(/charmming/images/wiki/smallstarred.gif)";
 }
 var table_name = tableobj.id.split('_');
 var topic_title = table_name[0];
 document.getElementById(topic_title + "_currate").innerHTML = "Rating:";
}

function changeRating(tdobj)
{
 var starvals = tdobj.id.split('_');
 new Ajax.Request('/charmming/wiki/ratetopic/',{method:'post', asynchronous:true, parameters: {'topic_title':starvals[0],'user_rating':starvals[1],'board_name': boardname}});
 stars_table = document.getElementById(starvals[0] + '_stars');
 tds = stars_table.getElementsByTagName('td');
 //This essentially saves the state the user rated the topic as until they refresh it
 for(i = 0;i < tds.length ; i++)
 {
  var starvals2 = tds[i].id.split('_');
  if(parseInt(starvals2[1]) <=  parseInt(starvals[1]))
  {
   tds[i].style.backgroundImage="url(/charmming/images/wiki/smallstarred.gif)";
   tds[i].className="check";
  }
  else if(parseInt(starvals2[1]) >=  parseInt(starvals[1]))
  {
   tds[i].style.backgroundImage="url(/charmming/images/wiki/smallstargrey.gif)";
   tds[i].className="";
  }
 } 
}

