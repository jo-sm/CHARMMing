$(function() {
  $('form').on('submit', function(e) {
    $.ajax({
      url: '/login',
      data: $('form').serialize(),
      dataType: 'json',
      type: 'post'
    }).done(function(data) {
      // Successfully logged in, redirect
      window.location.href = '/';
    }).fail(function(jqXHR) {
      // Something went wrong!
    });
    e.preventDefault();
    return false;
  });
});
