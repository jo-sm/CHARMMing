(function ($) {
  var ready = $.fn.ready;
  $.fn.ready = function (fn) {
    if (this.context === undefined) {
      // The $().ready(fn) case.
      ready(fn($));
    } else if (this.selector) {
      if($(this.selector, this.context).length) {
        ready(fn($));
      }
    } else {
      ready($.proxy(function(){
        $(this).each(fn);
      }, this));
    }
  }
})(jQuery);
