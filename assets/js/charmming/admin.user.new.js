jQuery('.admin-user-new').ready(function($) {
  $('form').on('submit', function(event) {
    $('.alert').addClass('hide');
    $.ajax({
      url: '/admin/user/new',
      data: $('form').serialize(),
      type: 'post',
      dataType: 'json'
    }).done(function() {
      window.href = '/admin/user';
    }).fail(function(jqXHR) {
      $('.alert').removeClass('hide').text(jqXHR.responseJSON.error);
    });
    event.preventDefault();
  });

  $('#random').on('click', function() {
      var text = "";
      var possible = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
      for (var i = 0; i < 8; i++)
        text += possible.charAt(Math.floor(Math.random() * possible.length));
      $('input[name="password"]').val(text);
  });
});
