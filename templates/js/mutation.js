//For mutationselect.html
  function getSegID(segment){ //Returns the index of the specified segment in segnames. A bit hacky, but it makes sending the form easier
    for (var i=0;i<segnames.length;i++){
      if(segnames[i] == segment)
        return i; 
    }
  }

  changeSegBox(); //We call these before anything executes so the initial choice gets itself updated
  changeResBox();
  
  //Following two are the event handlers...
  
  function changeSegBox(){ //Changes the residue list based on which segment you choose
    var currentseg = getSegID($("#MutSegi").val());
    $("#MutResName").html("");
    var resnameString = ""
    if (!(jQuery.isEmptyObject(myDictList[currentseg]))){    //Don't get keys on an empty dict...
      var sortlist = []
        for(var key in myDictList[currentseg]){
          sortlist[sortlist.length] = key //Lovely little hack this is
        }
      sortlist.sort()
      for (i=0;i<sortlist.length;i++){
        resnameString = resnameString + '<option value="' + sortlist[i] + '">' + sortlist[i] + '</option>';
      }
      $("#MutResName").html(resnameString);
    }
    changeResBox(); //Maybe this helps...it seems change listeners don't stack
  };


  function changeResBox(){
      var currentseg = getSegID($("#MutSegi").val());
      var currentres = $("#MutResName").val();
      $("#MutResi").val("");
      var resiString = ""
      //Sanity check so we don't try to handle things with an empty dict
      if (!(jQuery.isEmptyObject(myDictList[currentseg]))){
        resids = myDictList[currentseg][currentres]; //Get entries under key
        for (var i=0;i<resids.length;i++){ //resids is a list, can't do for...in
          resiString = resiString + '<option value="' + resids[i] + '">' + resids[i] + '</option>';
        }
      $("#MutResi").html(resiString);
      }
      changeResNumBox();
  };

  function changeResNumBox(){
    var segletter = $("#MutSegi").val()[0]; //first letter, to get chains...
    var resname = $("#MutResName").val(); 
    var resnum = $("#MutResi").val();
    var currentsegid = getSegID($("#MutSegi").val());
    if (segletter != null){
      if (currentsegid < (segnames.length - 1)){
        var script = "select " + resnum + " and atomno < " + chain_terminators[currentsegid + 1] + " and atomno >= " + chain_terminators[currentsegid];
      }else{
        var script = "select " + resnum + " and atomno >= " + chain_terminators[currentsegid];
      }
      Jmol.script(jmolApplet, script);
    }
  }
    
  //And these are the eventlisteners.
  document.getElementById("MutResName").addEventListener('change', changeResBox);
  document.getElementById("MutSegi").addEventListener('change', changeSegBox);
  document.getElementById("MutResi").addEventListener('change', changeResNumBox);
  $('input:radio').on('change', function(){
      Jmol.script(jmolApplet, "selectionHalos off;load " + filepath + "-" + this.id + ".pdb;select none;selectionHalos on;");
      changeSegBox(); //We call these so the right things get selected
      changeResBox();
      });
