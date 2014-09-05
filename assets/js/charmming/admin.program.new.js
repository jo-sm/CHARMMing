jQuery('.admin-program-new').ready(function() {
  var path = '';
  $('#verify').on('click', function() {
    path = $('input[name="program_path"]').val();
    $.ajax({
      url: '/admin/programs/verify',
      type: 'post',
      dataType: 'json',
      data: 'path=' + path
    }).done(function() {
      $('#verify').removeClass('secondary').addClass('success').html('Valid <i class="ion-checkmark-round"></i>');
    }).fail(function() {
      $('#verify').removeClass('secondary').addClass('alert').html('Invalid <i class="ion-close-round"></i>');
    });
  });

  $('input[name="program_path"]').on('keyup', function() {
    if ($(this).val() !== path) {
      $('#verify').removeClass('success alert').addClass('secondary').text('Verify');
    }
  });

  $('form').on('submit', function(e) {
    $('.alert').addClass('hide');
     $.ajax({
      url: '/admin/programs/new',
      data: $('form').serialize(),
      dataType: 'json',
      type: 'post'   
    }).done(function() {
      window.location.href = '/admin';
    }).fail(function(xhr) {
      $('.alert').removeClass('hide').text(xhr.responseJSON.error);
    });
    e.preventDefault();
  });
});
