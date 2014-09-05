jQuery('.structure-new').ready(function() {
  $('.structure-upload input[type="radio"]').on('click', function() {
    $('.structure-upload .upload-options').addClass('active');
    $('.structure-upload .upload-options > div').removeClass('selected');
    $('.structure-upload .upload-options > div[id="' + $(this).attr('id') + '"]').addClass('selected');
  });

  var residue_codes = {
    'ASP': 'ASP - Aspartic Acid',
    'ASPP': 'ASPP - Aspartic Acid',
    'GLU': 'GLU - Glutamic Acid',
    'GLUP': 'GLUP - Glutamic Acid',
    'HIS': 'HIS - Histidine',
    'HSP': 'HSP - δ+ε Histidine',
    'ARG': 'ARG - Arginine',
    'LYS': 'LYS - Lysine',
    'ALA': 'ALA - Alanine',
    'ASN': 'ASN - Asparagine',
    'CYS': 'CYS - Cysteine',
    'GLN': 'GLN - Glutamine',
    'GLY': 'GLY - Glycine',
    'HSD': 'HSD - δ Histidine',
    'HSE': 'HSE - ε Histidine',
    'ILE': 'ILE - Isoleucine',
    'LEU': 'LEU - Leucine',
    'LSN': 'LSN - Lysine',
    'MET': 'MET - Methionine',
    'PHE': 'PHE - Phenylaniline',
    'PRO': 'PRO - Proline',
    'SER': 'SER - Serine',
    'THR': 'THR - Threonine',
    'TRP': 'TRP - Tryptophan',
    'TYR': 'TYR - Tyrosine',
    'VAL': 'VAL - Valine'
  };
  $('#sequence_textarea').textcomplete([{
    match: /[AaCcGgHhLlIiMmPpSsTtVv](\w*)$/,
    search: function(term, callback) {
      list = [];
      exp = new RegExp('^' + term.toUpperCase());
      for (var code in residue_codes)
        if (code.match(exp))
          list.push(code);
      callback(list);
    },
    replace: function(value) {
      return value + ' ';
    },
    template: function(value) {
      return residue_codes[value];
    },
    index: 0
  }]);
 /* 
  $('#sequence_textarea').on('keydown', function(event) {
    var lastInput = $(this).val().split(' ').pop();
    // We should not search unless it matches a known code
    if (lastInput.match(/^[AaCcGgHhLlIiMmPpSsTtVv]/)) {
      match = new RegExp('^' + lastInput.toUpperCase());
      list = [];
      for(var code in residue_codes) {
        if (code.match(match)) {
          list.push(residue_codes[code]);
        }
      }
      console.log(list);
    }
  });
  */
});
