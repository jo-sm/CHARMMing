jQuery('.admin-task-new').ready(function() {
  var ractive = new Ractive({
    el: 'parameters',
    template: '#parameters-form',
    complete: function() {
      $.getJSON('/admin/program_sets').done(function(data) {
        ractive.set('program_sets', data);
      });
    },
    data: {
      groups: [],
      fieldValues: function() {
        if (!this.get('currentConditional')) return;
        var field = this.get(this.get('currentConditional')).field.split('-');
        return this.get('groups.' + (field[0]-1) + '.fields.' + (field[1]-1) + '.type.values');
      }
    },
    computed: {
      fieldType: function() {
        var field = this.get(this.get('currentConditional')).field.split('-');
        return this.get('groups.' + (field[0]-1) + '.fields.' + (field[1]-1) + '.type.name');
      }
    }
  });
  

  ractive.on('new-field', function(event, id) {
    var fields = this.get('groups.' + (id-1) + '.fields');
    fields.push({
      name: '',
      displayName: 'Display Name',
      id: fields.length+1,
      description: '',
      type: {
        name: 'text'
      },
      required: false,
      conditionals: [],
      order: fields.length+1,
      active: false
    });
    this.set('groups.' + (id-1) + '.fields', fields);
  });

  ractive.on('new-group', function() {
    var groups = this.get('groups');
    groups.push({
      id: groups.length + 1,
      order: groups.length + 1,
      name: 'Group Name',
      fields: []
    });
    this.set('groups', groups);
  });

  /* There is a bug where if you delete group #1 then #2, it stays as a ghost.
   * Need to do more thorough removal of ids and orders
   */
  ractive.on('delete', function(event, id) {
    var groups = this.get('groups');
    groups.splice((id-1), 1);
    this.set('groups', groups);
  });

  ractive.on('copy', function(event, id) {
    var group = this.get('groups.' + (id-1));
    var groups = this.get('groups');
    var copy = jQuery.extend({}, group);
    copy.active = false;
    copy.id = groups.length+1;
    copy.order = groups.length+1;
    groups.push(copy);
    this.set('groups', groups);
  });

  ractive.on('select-group', function(event, id) {
    // I'm not sure if I want to keep this syntax
    // It's a little difficult to follow for the non-JS and DOM aware developers, but it gets the job done without the need for jQuery
    // TODO: since event.node returns a DOM element, wrap that in a jQuery/Zepto wrapper
    this.set('groups.*.active', false);
    this.set('groups.' + (id-1) + '.active', true);
    this.set('selectedGroup', id);
    for(var i = 0; i < event.node.parentNode.children.length; i++) {
      var classes = event.node.parentNode.children[i].childNodes[0].className.split(' ');
      if (classes.indexOf('active') !== -1) {
        classes.splice(classes.indexOf('active'), 1);
      }
      event.node.parentNode.children[i].childNodes[0].className = classes.join(' ');
    }
    event.node.childNodes[0].className = event.node.childNodes[0].className + ' active';
  });

  ractive.on('select-field', function(event, id, groupId) {
    this.set('groups.' + (groupId-1) + '.fields.*.active', false);
    this.set('groups.' + (groupId-1) + '.fields.' + (id-1) + '.active', true);
    for(var i = 0; i < event.node.parentNode.parentNode.children.length; i++) {
      var classes = event.node.parentNode.parentNode.children[i].childNodes[0].className.split(' ');
      if (classes.indexOf('active') !== -1) {
        classes.splice(classes.indexOf('active'), 1);
      }
      event.node.parentNode.parentNode.children[i].childNodes[0].className = classes.join(' ');
    }
    event.node.className = event.node.className + ' active';
  });

  ractive.on('submit', function(event) {
    $.ajax({
      url: '/admin/tasks/new',
      type: 'post',
      dataType: 'json',
      data: "data=" + JSON.stringify(this.data)
    }).done(function(data) {
      window.location.href = '/admin/tasks';
    });
    event.original.preventDefault();
  });

  ractive.on('add-radio-value', function(el, event) {
    var data = this.get(el.keypath);
    if(!data.values) data.values = [];
    data.values.push({});
    this.set(el.keypath, data);
  });

  ractive.on('add-conditional', function(el, event) {
    var data = this.get(el.keypath);
    if(!data.conditionals) data.conditionals = [];
    data.conditionals.push({});
    this.set(el.keypath, data);
  });

  ractive.on('update-current-conditional', function(event) {
    this.set('currentConditional', event.keypath);
  });

  ractive.on('update-conditional', function(event, currentId, i, attribute) {
    currentId = currentId.split('-');
    currentId[0] = currentId[0]-1;
    currentId[1] = currentId[1]-1;
    var data = this.get('groups.' + currentId[0] + '.fields.' + currentId[1]);
    data.conditionals[i-1][attribute] = event.node.value;
    this.set('groups.' + currentId[0] + '.fields.' + currentId[1], data);
  });

  ractive.observe('groups.*.fields.*.order', function(newValue, oldValue, keypath) {
    var path = keypath.split('.');
    path.pop(); 
    path.pop();
    path = path.join('.');
    if(typeof oldValue !== 'undefined') {
      for(var i = 0; i < this.get(path).length; i++) { // Fields in this group
        var data = this.get(path)[i].order;
        if(data == newValue && keypath !== path + '.' + i + '.order') {
          this.set(path + '.' + i + '.order', oldValue);
        }
      }
    }
  });

  ractive.observe('groups.*.order', function(newValue, oldValue, keypath) {
    if(typeof oldValue !== 'undefined') {
      for(var i = 0; i < this.get('groups').length; i++) {
        var data = this.get('groups')[i].order;
        if(data == newValue && keypath !== 'groups.' + i + '.order') {
          this.set('groups.' + i + '.order', oldValue);
        }
      }
    }
  });

  /* It's necessary to observe this keypath due to the use of a pseudo-computed property
   * that happens when a Type of Input value is changed (i.e. adding a property to the 
   * radio button.
   */
  ractive.observe('groups.*.fields.*.type', function(newValue, oldValue, keypath) {
    // Call fieldValues (http://docs.ractivejs.org/latest/computed-properties)
    this.get('fieldValues').call(this);
  });
});
