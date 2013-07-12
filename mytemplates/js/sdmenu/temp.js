// $Id: common.js,v 1.25 2006/05/02 17:28:39 tvroom Exp $

function createPopup(xy, titlebar, name, contents, message) {
        var body = document.getElementsByTagName("body")[0];
        var div = document.createElement("div");
        div.id = name + "-popup";
        div.style.position = "absolute";

        var leftpos = xy[0] + "px";
        var toppos  = xy[1] + "px";

        div.style.left = leftpos;
        div.style.top = toppos;
        div.style.zIndex = "100";
        contents = contents || "";
        message  = message || "";

        div.innerHTML = '<div id="' + name + '-title" class="popup-title">' + titlebar + '</div>' +
                        '<div id="' + name + '-contents" class="popup-contents">' + contents + '</div>' +
                        '<div id="' + name + '-message" class="popup-message">' + message + '</div>';

        body.appendChild(div);
        div.className = "popup";
        return div;
}

function createPopupButtons() {
        var buttons = "";
        if (arguments.length > 0) {
                buttons = '<span class="buttons">';
        }
        for (var i=0; i<arguments.length; i++) {
                buttons =  buttons + "<span>" + arguments[i] + "</span>";
        }

        buttons = buttons + "</span>";
        return buttons;
}

function closePopup(id, refresh) {
        var el = $(id);
        if (el) {
                el.parentNode.removeChild(el);
        }
        if (refresh) {
                window.location.reload();
        }
}

function handleEnter(ev, func, arg) {
        if (!ev) {
                ev = window.event;
        }
        var code = ev.which || ev.keyCode;
        if (code == 13) { // return/enter
                func(arg);
                ev.returnValue = true;
                return true;
        }
        ev.returnValue = false;
        return false;
}


function moveByObjSize(div, addOffsetWidth, addOffsetHeight) {
        if (addOffsetWidth) {
                div.style.left = parseFloat(div.style.left || 0) + (addOffsetWidth * div.offsetWidth) + "px";
        }
        if (addOffsetHeight) {
                div.style.top = parseFloat(div.style.top || 0) + (addOffsetHeight * div.offsetHeight) + "px";
        }
}

function moveByXY(div, x, y) {
        if (x) {
                div.style.left = parseFloat(div.style.left || 0) + x + "px";
        }
        if (y) {
                div.style.top = parseFloat(div.style.top || 0) + y + "px";
        }
}

function getXYForId(id, addWidth, addHeight) {
        var div = $(id);
        var xy = Position.cumulativeOffset(div);
        if (addWidth) {
                xy[0] = xy[0] + div.offsetWidth;
        }
        if (addHeight) {
                xy[1] = xy[1] + div.offsetHeight;
        }
        return xy;
}

function toggleIntro(id, toggleid) {
        var obj = $(id);
        var toggle = $(toggleid);
        if (obj.className == 'introhide') {
                obj.className = "intro"
                toggle.innerHTML = "[-]";
        } else {
                obj.className = "introhide"
                toggle.innerHTML = "[+]";
        }
}

function tagsToggleStoryDiv(id, is_admin, type) {
        var bodyid = 'toggletags-body-' + id;
        var tagsbody = $(bodyid);
        if (tagsbody.className == 'tagshide') {
                tagsShowBody(id, is_admin, '', type);
        } else {
                tagsHideBody(id);
        }
}

function tagsHideBody(id) {
        // Make the body of the tagbox vanish
        var tagsbodyid = 'toggletags-body-' + id;
        var tagsbody = $(tagsbodyid);
        tagsbody.className = "tagshide"

        // Make the title of the tagbox change back to regular
        var titleid = 'tagbox-title-' + id;
        var title = $(titleid);
        title.className = "tagtitleclosed";

        // Make the tagbox change back to regular.
        var tagboxid = 'tagbox-' + id;
        var tagbox = $(tagboxid);
        tagbox.className = "tags";

        // Toggle the button back.
        var tagsbuttonid = 'toggletags-button-' + id;
        var tagsbutton = $(tagsbuttonid);
        tagsbutton.innerHTML = "[+]";
}

function tagsShowBody(id, is_admin, newtagspreloadtext, type) {

        type = type || "stories";

        //alert("Tags show body / Type: " + type );

        // Toggle the button to show the click was received
        var tagsbuttonid = 'toggletags-button-' + id;
        var tagsbutton = $(tagsbuttonid);
        tagsbutton.innerHTML = "[-]";

        // Make the tagbox change to the slashbox class
        var tagboxid = 'tagbox-' + id;
        var tagbox = $(tagboxid);
        tagbox.className = "tags";

        // Make the title of the tagbox change to white-on-green
        var titleid = 'tagbox-title-' + id;
        var title = $(titleid);
        title.className = "tagtitleopen";

        // Make the body of the tagbox visible
        var tagsbodyid = 'toggletags-body-' + id;
        var tagsbody = $(tagsbodyid);

        tagsbody.className = "tagbody";

        // If the tags-user div hasn't been filled, fill it.
        var tagsuserid = 'tags-user-' + id;
        var tagsuser = $(tagsuserid);
        if (tagsuser.innerHTML == "") {
                // The tags-user-123 div is empty, and needs to be
                // filled with the tags this user has already
                // specified for this story, and a reskey to allow
                // the user to enter more tags.
                tagsuser.innerHTML = "Retrieving...";
                var params = [];
                if (type == "stories") {
                        params['op'] = 'tags_get_user_story';
                        params['sidenc'] = id;
                } else if (type == "urls") {
                        //alert('getting user urls ' + id);
                        params['op'] = 'tags_get_user_urls';
                        params['id'] = id;
                }
                params['newtagspreloadtext'] = newtagspreloadtext
                ajax_update(params, tagsuserid);
                //alert('after ajax_update ' + tagsuserid);

                // Also fill the admin div.  Note that if the user
                // is not an admin, this call will not actually
                // return the necessary form (which couldn't be
                // submitted anyway).  The is_admin parameter just
                // saves us an ajax call to find that out, if the
                // user is not actually an admin.
                if (is_admin) {
                        var tagsadminid = 'tags-admin-' + id;
                        params = [];
                        if (type == "stories") {
                                params['op'] = 'tags_get_admin_story';
                                params['sidenc'] = id;
                        } else if (type == "urls") {
                                params['op'] = 'tags_get_admin_url';
                                params['id'] = id;
                        }
                        ajax_update(params, tagsadminid);
                }

        } else {
                if (newtagspreloadtext) {
                        // The box was already open but it was requested
                        // that we append some text to the user text.
                        // We can't do that by passing it in, so do it
                        // manually now.
                        var textinputid = 'newtags-' + id;
                        var textinput = $(textinputid);
                        textinput.value = textinput.value + ' ' + newtagspreloadtext;
                }
        }
}

function tagsOpenAndEnter(id, tagname, is_admin, type) {
        // This does nothing if the body is already shown.
        tagsShowBody(id, is_admin, tagname, type);
}

function reportError(request) {
        // replace with something else
        alert("error");
}

function tagsCreateForStory(id) {
        var toggletags_message_id = 'toggletags-message-' + id;
        var toggletags_message_el = $(toggletags_message_id);
        toggletags_message_el.innerHTML = 'Saving tags...';

        var params = [];
        params['op'] = 'tags_create_for_story';
        params['sidenc'] = id;
        var newtagsel = $('newtags-' + id);
        params['tags'] = newtagsel.value;
        var reskeyel = $('newtags-reskey-' + id);
        params['reskey'] = reskeyel.value;

        ajax_update(params, 'tags-user-' + id);

        // XXX How to determine failure here?
        toggletags_message_el.innerHTML = 'Tags saved.';
}

function tagsCreateForUrl(id) {
        var toggletags_message_id = 'toggletags-message-' + id;
        var toggletags_message_el = $(toggletags_message_id);
        toggletags_message_el.innerHTML = 'Saving tags...';

        var params = [];
        params['op'] = 'tags_create_for_url';
        params['id'] = id;
        var newtagsel = $('newtags-' + id);
        params['tags'] = newtagsel.value;
        var reskeyel = $('newtags-reskey-' + id);
        params['reskey'] = reskeyel.value;

        ajax_update(params, 'tags-user-' + id);

        // XXX How to determine failure here?
        toggletags_message_el.innerHTML = 'Tags saved.';
}

// helper functions
function ajax_update(params, onsucc, options, url) {
        var h = $H(params);

        if (!url)
                url = 'http://slashdot.org/ajax.pl';

        if (!options)
                options = {};

        options.method = 'post';
        options.parameters = h.toQueryString();

        var ajax = $.ajax({
                url:url,
                type:"post",
                data:h.toQueryString(),
                success: onsucc
                });
        //new Ajax.Updater(
        //        { success: onsucc },
        //        url,
        //        options
        //);
}

function ajax_periodic_update(secs, params, onsucc, options, url) {
        var h = $H(params);

        if (!url)
                url = 'http://slashdot.org/ajax.pl';

        if (!options)
                options = {};

        options.frequency = secs;
        options.method = 'post';
        options.parameters = h.toQueryString();

        var ajax = $.PeriodicalUpdater(url, {method: 'post', data=h.toQueryString(), multiplier:5, minTimeout:100, maxTimeout:(secs * 1000)}, onsucc);
//        var ajax = new Ajax.PeriodicalUpdater({ success: onsucc }, url, options);
}


/* Override Moo.fx cookie function, so that the cookie works across subdomains */
var Remember = new Object();
Remember = function(){};
Remember.prototype = {
        setCookie: function(value) {
                var date = new Date();
                date.setTime(date.getTime()+(this.days*24*60*60*1000));
                var expires = "; expires="+date.toGMTString();
                document.cookie = this.el+this.el.id+this.prefix+"="+value+expires+"; path=/; domain=.slashdot.org";
        }
}

function getElementsByClass(searchClass, node, tag) {
  var classElements = new Array();
  if (node == null) node = document;
  if (tag == null) tag = '*';
  var els = node.getElementsByTagName(tag);
  var elsLen = els.length;
  var pattern = new RegExp("(^|\s)" + searchClass + "(\s|$)");
  for (i = 0, j = 0; i < elsLen; i++) {
    if (pattern.test(els[i].className))
      classElements[j++] = els[i];
  }
  return classElements;
}

function toggleClass(element, firstClass, secondClass) {
  var removeClass = (Element.hasClassName(element, firstClass) ? firstClass : secondClass)
  var addClass = (Element.hasClassName(element, firstClass) ? secondClass : firstClass)

  Element.removeClassName(element, removeClass);
  Element.addClassName(element, addClass);
}

/* Available options: remember, expanded, collapsed, none */
var menuConfig = ['remember', 'expanded', 'remember', 'remember', 'remember', 'expanded'];

function init() {
  var blocks = $S('div#links div.block');
  blocks.each(function(block, index) {
    block.id = 'block' + index;
    var title = getElementsByClass('title', block, 'div')[0];
    var content = getElementsByClass('content', block, 'div')[0];
    if (menuConfig[index] == 'remember') {
      var slideEffect = new fx.RememberHeight(content, {duration: 300});
      Element.addClassName(block, parseInt(slideEffect.readCookie()) == 0 ? 'collapsed' : 'expanded');
    } else if (menuConfig[index] == 'expanded') {
      var slideEffect = new fx.Height(content, {duration: 300});
      Element.addClassName(block, 'expanded');
    } else if (menuConfig[index] == 'collapsed') {
      var slideEffect = new fx.Height(content, {duration: 300});
      slideEffect.hide();
      Element.addClassName(block, 'collapsed');
    } else {
      var slideEffect = null;
      Element.addClassName(block, 'static');
    }

    if (slideEffect) {
      title.onclick = function(event) {
        try { var element = Event.element(event); }
        catch(err) { var element = Event.element(window.event); }
        if (element.tagName != 'A') {
          slideEffect.toggle();
          toggleClass(this.parentNode, 'collapsed', 'expanded');
        }
      }
    }
  });

  var loginLink = $('loginlink');
  if (loginLink) {
    loginLink.onclick = function() {
      var slideEffect = new fx.Height('loginform', {duration: 250});
      slideEffect.toggle();
      toggleClass(this, 'collapsed', 'expanded');
      return false;
    }
  }
}
Event.observe(window, 'load', init);
