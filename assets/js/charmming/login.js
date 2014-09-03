jQuery('html.login').ready(function() {
  $('form.login').on('submit', function(e) {
    $('.alert').addClass('hide');
    $.ajax({
      url: '/login',
      data: $('form').serialize(),
      dataType: 'json',
      type: 'post'
    }).done(function(data) {
      // Successfully logged in, redirect
      window.location.href = '/';
    }).fail(function(jqXHR) {
      $('.alert').removeClass('hide').text(jqXHR.responseJSON.error);
    });
    e.preventDefault();
    return false;
  });
});
