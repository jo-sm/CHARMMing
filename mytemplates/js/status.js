var scroll_to_bottom = true;

function dialog_open(task_id,task_status){
  if ($("#dialog_output_"+task_id).dialog().length > 0){
    $("#dialog_output_"+task_id).PeriodicalUpdater("stop");
    $("#dialog_output_"+task_id).dialog("close"); //If it exists, close it.
  }
  var opt = {
                      resizable:true,
                      height:600,
                      width:800,
                      modal:false,
                      autoOpen:false,
                      buttons:{
                        "Close":function(){
                          $(this).dialog("close");
                        }
                      }};
  if(task_status == "Running"){
    opt['buttons'] = {
                        "Disable autoscroll":function(){
                            scroll_to_bottom = false;
                        },
                        "Enable autoscroll":function(){
                            scroll_to_bottom = true;
                          },
                        "Close":function(){
                          $(this).dialog("close");
                        }
  };
  }

  $("#dialog_output_"+task_id).dialog(opt).dialog("open");
  //No idea if scope is messing me up here
  if (task_status == "Running"){
    var updater =  $("#dialog_output_"+task_id).PeriodicalUpdater('/charmming/view_task_output/'+task_id+'/', {method:'get', multiplier:2, minTimeout:100, maxTimeout:500}, function(remoteData){$("#dialog_output_"+task_id).html(remoteData);if(scroll_to_bottom){$("#dialog_output_"+task_id).scrollTop($("#dialog_output_"+task_id)[0].scrollHeight);}});
  }else{
    $.ajax({
                  url:"/charmming/view_task_output/"+task_id+"/",
                  type:"get",
                  async:true,
                  success:function(responseData){
                  $("#dialog_output_"+task_id).html(responseData);
                  }
                  });
  }
  }
  
