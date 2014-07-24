jQuery('.structure-index').ready(function($) {
  $('.delete').on('click', function() {
    $this = $(this);
    // Ask if user is sure they want to delete!
    $.ajax({
      url: '/structure/' + $this.attr('id'),
      type: 'delete',
      dataType: 'json'
    }).done(function(data) {
      $this.closest('tr').remove();
    });
  });
});
