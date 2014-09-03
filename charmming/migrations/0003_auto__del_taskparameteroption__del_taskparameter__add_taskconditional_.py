# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting model 'TaskParameterOption'
        db.delete_table('charmming_taskparameteroption')

        # Deleting model 'TaskParameter'
        db.delete_table('charmming_taskparameter')

        # Adding model 'TaskConditional'
        db.create_table('charmming_taskconditional', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('object_type', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('object_id', self.gf('django.db.models.fields.IntegerField')()),
            ('group', self.gf('django.db.models.fields.IntegerField')()),
            ('field', self.gf('django.db.models.fields.IntegerField')()),
            ('condition', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('thenDo', self.gf('django.db.models.fields.CharField')(max_length=15)),
        ))
        db.send_create_signal('charmming', ['TaskConditional'])

        # Adding model 'WithConditional'
        db.create_table('charmming_withconditional', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
        ))
        db.send_create_signal('charmming', ['WithConditional'])

        # Adding model 'TaskGroup'
        db.create_table('charmming_taskgroup', (
            ('withconditional_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['charmming.WithConditional'], unique=True, primary_key=True)),
            ('parent_type', self.gf('django.db.models.fields.CharField')(default=None, max_length=15)),
            ('parent_id', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('task', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.Task'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=45)),
            ('order', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal('charmming', ['TaskGroup'])

        # Adding model 'TaskField'
        db.create_table('charmming_taskfield', (
            ('withconditional_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['charmming.WithConditional'], unique=True, primary_key=True)),
            ('group', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.Group'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('display_name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=150)),
            ('order', self.gf('django.db.models.fields.IntegerField')()),
            ('required', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=15)),
        ))
        db.send_create_signal('charmming', ['TaskField'])

        # Adding model 'TaskFieldOption'
        db.create_table('charmming_taskfieldoption', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('task_field', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.TaskField'])),
            ('label', self.gf('django.db.models.fields.CharField')(max_length=25)),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=25)),
        ))
        db.send_create_signal('charmming', ['TaskFieldOption'])

        # Adding model 'TaskType'
        db.create_table('charmming_tasktype', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
        ))
        db.send_create_signal('charmming', ['TaskType'])

        # Deleting field 'Task.position'
        db.delete_column('charmming_task', 'position')

        # Adding field 'Task.order'
        db.add_column('charmming_task', 'order',
                      self.gf('django.db.models.fields.IntegerField')(default=1),
                      keep_default=False)


        # Renaming column for 'Task.type' to match new field type.
        db.rename_column('charmming_task', 'type', 'type_id')
        # Changing field 'Task.type'
        db.alter_column('charmming_task', 'type_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.TaskType']))
        # Adding index on 'Task', fields ['type']
        db.create_index('charmming_task', ['type_id'])


        # Changing field 'JobTaskParameter.task_parameter'
        db.alter_column('charmming_jobtaskparameter', 'task_parameter_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.TaskField']))

    def backwards(self, orm):
        # Removing index on 'Task', fields ['type']
        db.delete_index('charmming_task', ['type_id'])

        # Adding model 'TaskParameterOption'
        db.create_table('charmming_taskparameteroption', (
            ('task_parameter', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.TaskParameter'])),
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=25)),
            ('label', self.gf('django.db.models.fields.CharField')(max_length=25)),
        ))
        db.send_create_signal('charmming', ['TaskParameterOption'])

        # Adding model 'TaskParameter'
        db.create_table('charmming_taskparameter', (
            ('task', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.Task'])),
            ('display_name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=150)),
            ('required', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=15)),
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
        ))
        db.send_create_signal('charmming', ['TaskParameter'])

        # Deleting model 'TaskConditional'
        db.delete_table('charmming_taskconditional')

        # Deleting model 'WithConditional'
        db.delete_table('charmming_withconditional')

        # Deleting model 'TaskGroup'
        db.delete_table('charmming_taskgroup')

        # Deleting model 'TaskField'
        db.delete_table('charmming_taskfield')

        # Deleting model 'TaskFieldOption'
        db.delete_table('charmming_taskfieldoption')

        # Deleting model 'TaskType'
        db.delete_table('charmming_tasktype')


        # User chose to not deal with backwards NULL issues for 'Task.position'
        raise RuntimeError("Cannot reverse this migration. 'Task.position' and its values cannot be restored.")
        
        # The following code is provided here to aid in writing a correct migration        # Adding field 'Task.position'
        db.add_column('charmming_task', 'position',
                      self.gf('django.db.models.fields.IntegerField')(),
                      keep_default=False)

        # Deleting field 'Task.order'
        db.delete_column('charmming_task', 'order')


        # Renaming column for 'Task.type' to match new field type.
        db.rename_column('charmming_task', 'type_id', 'type')
        # Changing field 'Task.type'
        db.alter_column('charmming_task', 'type', self.gf('django.db.models.fields.CharField')(max_length=25))

        # Changing field 'JobTaskParameter.task_parameter'
        db.alter_column('charmming_jobtaskparameter', 'task_parameter_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.TaskParameter']))

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
            'task_parameter': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.TaskField']"}),
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
            'order': ('django.db.models.fields.IntegerField', [], {}),
            'program_set': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.ProgramSet']"}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.TaskType']"}),
            'visible': ('django.db.models.fields.BooleanField', [], {'default': 'True'})
        },
        'charmming.taskconditional': {
            'Meta': {'object_name': 'TaskConditional'},
            'condition': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'field': ('django.db.models.fields.IntegerField', [], {}),
            'group': ('django.db.models.fields.IntegerField', [], {}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'object_id': ('django.db.models.fields.IntegerField', [], {}),
            'object_type': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'thenDo': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '20'})
        },
        'charmming.taskfield': {
            'Meta': {'object_name': 'TaskField', '_ormbases': ['charmming.WithConditional']},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '150'}),
            'display_name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'group': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.Group']"}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'order': ('django.db.models.fields.IntegerField', [], {}),
            'required': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'withconditional_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['charmming.WithConditional']", 'unique': 'True', 'primary_key': 'True'})
        },
        'charmming.taskfieldoption': {
            'Meta': {'object_name': 'TaskFieldOption'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '25'}),
            'task_field': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.TaskField']"}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '25'})
        },
        'charmming.taskgroup': {
            'Meta': {'object_name': 'TaskGroup', '_ormbases': ['charmming.WithConditional']},
            'name': ('django.db.models.fields.CharField', [], {'max_length': '45'}),
            'order': ('django.db.models.fields.IntegerField', [], {}),
            'parent_id': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'parent_type': ('django.db.models.fields.CharField', [], {'default': 'None', 'max_length': '15'}),
            'task': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.Task']"}),
            'withconditional_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['charmming.WithConditional']", 'unique': 'True', 'primary_key': 'True'})
        },
        'charmming.tasktype': {
            'Meta': {'object_name': 'TaskType'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'})
        },
        'charmming.user': {
            'Meta': {'object_name': 'User'},
            'email_address': ('django.db.models.fields.CharField', [], {'max_length': '65'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'password_digest': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'permissions': ('django.db.models.fields.IntegerField', [], {}),
            'username': ('django.db.models.fields.CharField', [], {'max_length': '20'})
        },
        'charmming.withconditional': {
            'Meta': {'object_name': 'WithConditional'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        'charmming.workingstructure': {
            'Meta': {'object_name': 'WorkingStructure'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'structure_id': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.Structure']"})
        }
    }

    complete_apps = ['charmming']