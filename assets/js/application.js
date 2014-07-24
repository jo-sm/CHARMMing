/* Dependencies:
 *= require vendor/modernizr
 *= require vendor/jquery
 *= require vendor/fastclick
 *= require vendor/jquery.textcomplete.min
 *= require vendor/jquery.ready
 *= require foundation.min
 *= require_tree charmming
 *= require_self
 */
(function($) {
  $(document).foundation();

  // Socket for notifications
  var sock = SockJS('http://131.247.213.199:8080/ns');
  sock.onopen = function() {
    console.log('Socket open');
  };

  sock.onmessage = function(event) {
    console.log('message', event.data);
  };

  sock.onclose = function() {
    console.log('close');
  };

  // Side menu expanded click
  $('.expand').on('click', function() {
    if ($(this).hasClass('ion-plus-round')) {
      // Is currently not expanded
      $(this).removeClass('ion-plus-round').addClass('ion-minus-round');
      $(this).parent().addClass('expanded');
    } else {
      // Is currently expanded
      $(this).removeClass('ion-minus-round').addClass('ion-plus-round');
      $(this).parent().removeClass('expanded');
    }
  });
})(jQuery);
