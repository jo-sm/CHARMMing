jQuery('.admin-task-new').ready(function($) {
  var parameter = Ractive.extend({
    template: '#parameter',
    setTemplate: function(options) {
      options.partials.element = options.data.template;
    },
    beforeInit: function(options) {
      this.setTemplate(options);
    }
  });

  var ractive = new Ractive({
    el: 'parameters',
    template: '#parameters-form',
    components: {
      parameter: parameter
    },
    data: {
      parameters: []
    }
  });

  ractive.on('new-parameter', function(el, event) {
    this.data.parameters.push({
      name: '',
      displayName: 'Display Name',
      id: this.data.parameters.length + 1,
      description: '',
      type: 'text',
      required: false,
      conditional: 'none',
      order: this.data.parameters.length + 1
    });

    this.update('parameters');
  });

  ractive.observe('parameters.*.order', function(newValue, oldValue, keypath) {
    if(typeof oldValue !== 'undefined') {
      for(var i = 0; i < this.data.parameters.length; i++) {
        if(this.data.parameters[i].order == newValue && keypath !== 'parameters.' + i + '.order') {
          this.data.parameters[i].order = oldValue;
        }
      }
      this.update('parameters');
    }
  });
});
