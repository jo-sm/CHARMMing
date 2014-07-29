# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting field 'TaskParameterOption.task_parameters'
        db.delete_column('charmming_taskparameteroption', 'task_parameters_id')

        # Adding field 'TaskParameterOption.task_parameter'
        db.add_column('charmming_taskparameteroption', 'task_parameter',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, to=orm['charmming.TaskParameter']),
                      keep_default=False)


    def backwards(self, orm):
        # Adding field 'TaskParameterOption.task_parameters'
        db.add_column('charmming_taskparameteroption', 'task_parameters',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, to=orm['charmming.TaskParameter']),
                      keep_default=False)

        # Deleting field 'TaskParameterOption.task_parameter'
        db.delete_column('charmming_taskparameteroption', 'task_parameter_id')


    models = {
        'charmming.group': {
            'Meta': {'object_name': 'Group'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '150'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '45'})
        },
        'charmming.inputscripttemplate': {
            'Meta': {'object_name': 'InputScriptTemplate'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'template': ('django.db.models.fields.TextField', [], {})
        },
        'charmming.job': {
            'Meta': {'object_name': 'Job'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '45'})
        },
        'charmming.jobtaskparameter': {
            'Meta': {'object_name': 'JobTaskParameter'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'job': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.Job']"}),
            'task_parameter': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.TaskParameter']"}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '25'})
        },
        'charmming.lesson': {
            'Meta': {'object_name': 'Lesson'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '150'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '45'})
        },
        'charmming.lessonstep': {
            'Meta': {'object_name': 'LessonStep'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'lesson': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.Lesson']"})
        },
        'charmming.program': {
            'Meta': {'object_name': 'Program'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '150'}),
            'enabled': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'input_flag': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'output_flag': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'path': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'program_set': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.ProgramSet']"})
        },
        'charmming.programparameter': {
            'Meta': {'object_name': 'ProgramParameter'},
            'flag': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'program': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.Program']"}),
            'use_equals_sign': ('django.db.models.fields.BooleanField', [], {'default': 'False'})
        },
        'charmming.programset': {
            'Meta': {'object_name': 'ProgramSet'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '150'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'program_set_template': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.InputScriptTemplate']"})
        },
        'charmming.session': {
            'Meta': {'object_name': 'Session'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'ip_address': ('django.db.models.fields.CharField', [], {'max_length': '32'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '64'}),
            'last_accessed': ('django.db.models.fields.DateTimeField', [], {}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.User']"}),
            'user_agent': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        'charmming.structure': {
            'Meta': {'object_name': 'Structure'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'directory': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'file_name': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'uploaded_at': ('django.db.models.fields.DateTimeField', [], {}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.User']"})
        },
        'charmming.task': {
            'Meta': {'object_name': 'Task'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'input_script_template': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.InputScriptTemplate']"}),
            'is_multi_step': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '45'}),
            'position': ('django.db.models.fields.IntegerField', [], {}),
            'program_set': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.ProgramSet']"}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'visible': ('django.db.models.fields.BooleanField', [], {'default': 'True'})
        },
        'charmming.taskparameter': {
            'Meta': {'object_name': 'TaskParameter'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '150'}),
            'display_name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'required': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'task': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.Task']"}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '15'})
        },
        'charmming.taskparameteroption': {
            'Meta': {'object_name': 'TaskParameterOption'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'task_parameter': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.TaskParameter']"}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '25'})
        },
        'charmming.user': {
            'Meta': {'object_name': 'User'},
            'email_address': ('django.db.models.fields.CharField', [], {'max_length': '65'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'password_digest': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'permissions': ('django.db.models.fields.IntegerField', [], {}),
            'username': ('django.db.models.fields.CharField', [], {'max_length': '20'})
        },
        'charmming.workingstructure': {
            'Meta': {'object_name': 'WorkingStructure'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'structure_id': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.Structure']"})
        }
    }

    complete_apps = ['charmming']