/* Dependencies:
 *= require vendor/modernizr
 *= require vendor/jquery
 *= require vendor/fastclick
 *= require vendor/jquery.textcomplete.min
 *= require vendor/jquery.ready
 *= require_tree charmming
 *= require_self
 */
(function($) {

  $('[data-dropdown]').each(function(el) {
    var dropdown = $(el).next('ul');
    dropdown.on('mouseover', function() {
      $(this).show();
    });
    dropdown.on('mouseout', function() {
      $(this).hide();
    });
  });

  // Socket for notifications
  var sock = SockJS('http://131.247.213.199:8080/ns');
  sock.onopen = function() {
  };

  sock.onmessage = function(event) {
  };

  sock.onclose = function() {
  };

  $('.tabs').on('click', '.tab', function(e) {
    $('.tab').removeClass('active');
    $(this).addClass('active');
    
    $('.tabs-content .tab-content').removeClass('active');
    $('.tabs-content ' + $(this).attr('href')).addClass('active');
    e.preventDefault();
  });

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
