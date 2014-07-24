jQuery('.working-structure-new').ready(function($) {
  // Load the selected segments and patches  
  loadSegments($('select#structure').val());

  $('select#structure').on('change', function() {
    loadSegments($(this).val());
  });

});

function loadSegments(id) {
  $('.loading').addClass('ion-loading-c');
  $('table > tbody').text('');
  $.getJSON('/structure/' + id + '/segments').done(function(data) {
    for (var model in data) {
      var _model = data[model];
      for (var segment in _model) {
        var _segment = _model[segment];
        $('table > tbody').append(patchTableRow(segment, _segment.segmentType, _segment.firstPatches, _segment.lastPatches, _segment.residue, _segment.additionalToppar));
      }
    }
    $('.loading').removeClass('ion-loading-c');
  });
}

function patchTableRow(segmentName, segmentType, firstPatches, lastPatches, residue, additionalToppar) {
  var str = '<tr>';
  str += '<td><input name="selected[' + segmentName + ']" type="checkbox"></td>'
  str += '<td>' + segmentName  + '</td>'; 
  str += '<td>' + segmentType + '</td>';
  str += '<td>' + (firstPatches.length == 0 ? 'None' : patchList(segmentName, 'first', firstPatches)) + '</td>';
  str += '<td>' + (lastPatches.length == 0 ? 'None' : patchList(segmentName, 'last', lastPatches)) + '</td>';
  str += '<td>' + (segmentType === 'bad' ? newTopology(segmentName, additionalToppar) : 'CHARMM36 CGenFF topology &amp; parameters' ) + '</td>';
  str += '<td class="center">' + ( residue ? '<strong>' + residue + '</strong><br><a href="http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/' + residue + '" target="_blank"><img src="http://www.rcsb.org/pdb/images/' + residue + '_120.gif"></a>' : '') + '</td>';
  str += '</tr>';
  return str;
}

function newTopology(id, additionalTopology) {
  var str = '<select id="toppar" name="toppar_' + id + '">';
  str += '<option value="auto">Attempt to automatically generate</option>';
  str += '<option value="upload">Upload my own</option>';
  for (var toppar in additionalTopology) {
    str += '<option value="' + toppar + '">' + additionalTopology[toppar] + '</option>';
  }
  str += '</select>';
  str += '<div class="upload-toppar">';
  str += '<div><label>Topology</label><input type="file" name="' + id + '_topology"></div>';
  str += '<div><label>Parameter</label><input type="file" name="' + id + '_parameter"></div>';
  str += '</div>';

  setTimeout(function() {
    $("select[name^='toppar_']").on('change', function() {
      if ($(this).val() == 'upload') {
        $(this).next('div').show();
      } else {
        $(this).next('div').hide();
      }
    });
  }, 5);

  return str;
}

function patchList(id, first_last, patches) {
  var str = '<select class="inline" name="' + id + '_' + first_last + '">';
  for (var i = 0; i < patches.length; i++) {
    str += '<option value="' + patches[i] + '">' + patches[i] + '</option>';
  }
  str += '</select>'
  return str;
}
