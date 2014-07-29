jQuery('.admin-task-new').ready(function($) {
  var ractive = new Ractive({
    el: 'parameters',
    template: '#parameters-form',
    data: {
      parameters: []
    }
  });

  ractive.on('new-parameter', function(el, event) {
    var parameters = this.get('parameters');
    parameters.push({
      name: '',
      displayName: 'Display Name',
      id: parameters.length + 1,
      description: '',
      type: {
        name: 'text'
      },
      required: false,
      conditional: 'none',
      order: parameters.length + 1
    });
    this.set('parameters', parameters);
  });

  ractive.on('add-radio-value', function(el, event) {
    console.log(el.keypath);
    var data = this.get(el.keypath);
    if(!data.values) data.values = [];
    data.values.push({});
    this.set(el.keypath, data);
  });

  ractive.observe('parameters.*.order', function(newValue, oldValue, keyPath) {
    if(typeof oldValue !== 'undefined') {
      for(var i = 0; i < this.get('parameters').length; i++) {
        var dataPath = 'parameters.' + i + '.order';
        if(this.get(dataPath) == newValue && keyPath !== dataPath) {
          this.set('parameters.' + i + '.order', oldValue);
        }
      }
    }
  });

  ractive.observe('parameters.*.type.name', function(newValue, oldValue, keyPath) {
      var dataPath = keyPath.split('.');
      dataPath = dataPath.slice(0, dataPath.length-1).join('.');
      var data = this.get(dataPath);
      this.update('parameters');
  });
});
