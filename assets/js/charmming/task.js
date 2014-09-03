jQuery('.task').ready(function($) {
  function conditionalArray(field, template) {
    var conditional = {
        display: {
            and: [],
            or: []
        },
        disable: {
            and: [],
            or: []
        }
    }
    
    for(var m = 0; m < field.conditionals.length; m++) {
        var _c = field.conditionals[m];
        
        var _cS = ''; // Conditional string
        
        // Sometimes a condition may not be "fully baked"
        if(_c.field && _c.thenDo) {
            var f = _c.field.split('-');
            f[0] = parseFloat(f[0]);
            f[1] = parseFloat(f[1]);
            
            var _cC = template.groups[(f[0]-1)].fields[(f[1]-1)];
            
            if(_cC.conditionals && _cC.conditionals.length) {
                var c = conditionalArray(_cC, template);
                conditional.display.and = c.display.or.concat(c.display.and);
                conditional.disable.and = c.disable.or.concat(c.disable.and);
            } 
            var _cS = _cC.name;
            
            switch(_c.condition) {
                case 'is':
                    if(_c.thenDo === 'show' || _c.thenDo === 'disable') _cS += ' == ';
                    if(_c.thenDo === 'hide') _cS += ' != ';
                    break;
                case 'is-not':
                    if(_c.thenDo === 'show' || _c.thenDo === 'disable') _cS += ' != ';
                    if(_c.thenDo === 'hide') _cS += ' == ';
                    break;
                case 'gt':
                    if(_c.thenDo === 'show' || _c.thenDo === 'disable') _cS += ' > ';
                    if(_c.thenDo === 'hide') _cS += ' < ';
                    break;
                case 'lt':
                    if(_c.thenDo === 'show' || _c.thenDo === 'disable') _cS += ' < ';
                    if(_c.thenDo === 'hide') _cS += ' > ';
                    break;
                case 'checked':
                    if(_c.thenDo === 'show' || _c.thenDo === 'disable') _cS += ' == true';
                    if(_c.thenDo === 'hide') _cS += ' == false';
                    break;
                case 'unchecked':
                    if(_c.thenDo === 'show' || _c.thenDo === 'disable') _cS += ' == false';
                    if(_c.thenDo === 'hide') _cS += ' == true';
                    break;
                    
            }
            
            if(_c.condition !== 'checked' && _c.condition !== 'unchecked')
                _cS += '"' + _c.value + '"';
            
            if(_c.thenDo !== 'disable') conditional.display.or.push(_cS);
            if(_c.thenDo === 'disable') conditional.disable.or.push(_cS);
        }
    }
    
    return conditional;
  }

  window.loadTemplate = function loadTemplate(template) {
    var renderedTemplate = '';
    var data = {};
    
    renderedTemplate += ('<h2>' + template.name + '</h2>');
    renderedTemplate += ('<h3>' + template.description + '</h2>');
    renderedTemplate += ('<form on-submit="submit">');
    for(var i = 0; i < template.groups.length; i++) {
        var group = template.groups[i];
        var conditional = null;
        
        // Handle group's conditionals
        if(group.conditionals.length) {
            conditional = conditionalArray(group, template);
        }
        
        var condStr = [];
        
        if(conditional) {
            condStr.push(conditional.display.and.join(' && '));
            if(conditional.display.and.length) {
              condStr.push('&&');
            }
            condStr.push(conditional.display.or.join(' || '));
            renderedTemplate += '{{#if (' + condStr.join(' ') + ') }}';
        }
        
        // Beginning of group        
        renderedTemplate += "<div class='group'>" + "\n";

        renderedTemplate += '<h4>' + group.name + '</h4>' + "\n";
        for(var j = 0; j < group.fields.length; j++) {
            var field = group.fields[j];
            var fConditional = null;
            
            // Handle field's conditionals
            if(field.conditionals) {
                fConditional = conditionalArray(field, template);
            }
                                    
            if(fConditional && (fConditional.display.or.length || fConditional.display.and.length)) {
                condStr = fConditional.display.or.join(' || ');
                if(fConditional.display.and.length) {
                     condStr += ' && ' + fConditional.display.and.join(' && ');
                }
                renderedTemplate += '{{#if (' + condStr + ')}}';
            }
            
            renderedTemplate += '<div>';
            renderedTemplate += field.displayName;
            
            condStr = '';
            if(fConditional && (fConditional.disable.and.length || fConditional.disable.or.length)) {
                condStr = fConditional.disable.or.join(' || ');
                if(fConditional.disable.and.length) {
                    condStr += ' && ' + fConditional.disable.and.join(' && ');
                }
            }

            
            if(field.type.name === 'text') {                
                renderedTemplate += '<input type="text" value="{{' + field.name + '}}"' + (condStr ? ' disabled={{(' + condStr + ')}}' : '' ) + '>' + "\n";
            } else if(field.type.name === 'radio') {
                for(var k = 0; k < field.type.values.length; k++) {
                    var value = field.type.values[k];
                    
                    renderedTemplate += value.name + ': <input type="radio" name="{{' + field.name + '}}" value="' + value.value + '"' + (condStr ? ' disabled={{(' + condStr + ')}}' : '' ) + '>' + "\n";
                }
            } else if(field.type.name === 'checkbox') {
                renderedTemplate += '<input type="checkbox" checked="{{' + field.name + '}}"' + (condStr ? ' disabled={{#if (' + condStr + ')}}disabled{{/if}}' : '' ) + '>' + "\n";
            }
            
            renderedTemplate += '</div>' + "\n";
            if(fConditional && (fConditional.display.and.length || fConditional.display.or.length)) renderedTemplate += '{{/if}}';
        }
        
        // End of group        
        renderedTemplate += '</div>' + "\n";
        
        // End of conditional
        if(conditional) renderedTemplate += '{{/if}}';
    }
    
    renderedTemplate += '<input type="submit"></form>';
    
    //console.log(renderedTemplate);
    
    var ractive =  new Ractive({
        el: '#rendered',
        template: renderedTemplate,
        data: data
    });

    ractive.on('submit', function(event) {
      console.log(this.data);
      event.original.preventDefault();
    });
  }
});
