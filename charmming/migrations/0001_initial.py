# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'User'
        db.create_table('charmming_user', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('username', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('password_digest', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('email_address', self.gf('django.db.models.fields.CharField')(max_length=65)),
            ('permissions', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal('charmming', ['User'])

        # Adding model 'Group'
        db.create_table('charmming_group', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=45)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=150)),
        ))
        db.send_create_signal('charmming', ['Group'])

        # Adding model 'Session'
        db.create_table('charmming_session', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.User'])),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=64)),
            ('user_agent', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('ip_address', self.gf('django.db.models.fields.CharField')(max_length=32)),
            ('last_accessed', self.gf('django.db.models.fields.DateTimeField')()),
        ))
        db.send_create_signal('charmming', ['Session'])

        # Adding model 'Structure'
        db.create_table('charmming_structure', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.User'])),
            ('directory', self.gf('django.db.models.fields.CharField')(max_length=25)),
            ('file_name', self.gf('django.db.models.fields.CharField')(max_length=25)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=25)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('uploaded_at', self.gf('django.db.models.fields.DateTimeField')()),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=15)),
        ))
        db.send_create_signal('charmming', ['Structure'])

        # Adding model 'WorkingStructure'
        db.create_table('charmming_workingstructure', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('structure_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.Structure'])),
        ))
        db.send_create_signal('charmming', ['WorkingStructure'])

        # Adding model 'InputScriptTemplate'
        db.create_table('charmming_inputscripttemplate', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('template', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal('charmming', ['InputScriptTemplate'])

        # Adding model 'ProgramSet'
        db.create_table('charmming_programset', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=150)),
            ('program_set_template', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.InputScriptTemplate'])),
        ))
        db.send_create_signal('charmming', ['ProgramSet'])

        # Adding model 'Program'
        db.create_table('charmming_program', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=150)),
            ('program_set', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.ProgramSet'])),
            ('path', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('enabled', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('input_flag', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('output_flag', self.gf('django.db.models.fields.CharField')(max_length=10)),
        ))
        db.send_create_signal('charmming', ['Program'])

        # Adding model 'ProgramParameter'
        db.create_table('charmming_programparameter', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('program', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.Program'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=25)),
            ('flag', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('use_equals_sign', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal('charmming', ['ProgramParameter'])

        # Adding model 'Task'
        db.create_table('charmming_task', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('program_set', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.ProgramSet'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=45)),
            ('input_script_template', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.InputScriptTemplate'])),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=25)),
            ('position', self.gf('django.db.models.fields.IntegerField')()),
            ('visible', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('is_multi_step', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal('charmming', ['Task'])

        # Adding model 'TaskParameter'
        db.create_table('charmming_taskparameter', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('task', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.Task'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('display_name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=150)),
            ('required', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=15)),
        ))
        db.send_create_signal('charmming', ['TaskParameter'])

        # Adding model 'TaskParameterOption'
        db.create_table('charmming_taskparameteroption', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('task_parameters', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.TaskParameter'])),
            ('label', self.gf('django.db.models.fields.CharField')(max_length=25)),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=25)),
        ))
        db.send_create_signal('charmming', ['TaskParameterOption'])

        # Adding model 'Job'
        db.create_table('charmming_job', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=45)),
        ))
        db.send_create_signal('charmming', ['Job'])

        # Adding model 'JobTaskParameter'
        db.create_table('charmming_jobtaskparameter', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('job', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.Job'])),
            ('task_parameter', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.TaskParameter'])),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=25)),
        ))
        db.send_create_signal('charmming', ['JobTaskParameter'])

        # Adding model 'Lesson'
        db.create_table('charmming_lesson', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=45)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=150)),
        ))
        db.send_create_signal('charmming', ['Lesson'])

        # Adding model 'LessonStep'
        db.create_table('charmming_lessonstep', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('lesson', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['charmming.Lesson'])),
        ))
        db.send_create_signal('charmming', ['LessonStep'])


    def backwards(self, orm):
        # Deleting model 'User'
        db.delete_table('charmming_user')

        # Deleting model 'Group'
        db.delete_table('charmming_group')

        # Deleting model 'Session'
        db.delete_table('charmming_session')

        # Deleting model 'Structure'
        db.delete_table('charmming_structure')

        # Deleting model 'WorkingStructure'
        db.delete_table('charmming_workingstructure')

        # Deleting model 'InputScriptTemplate'
        db.delete_table('charmming_inputscripttemplate')

        # Deleting model 'ProgramSet'
        db.delete_table('charmming_programset')

        # Deleting model 'Program'
        db.delete_table('charmming_program')

        # Deleting model 'ProgramParameter'
        db.delete_table('charmming_programparameter')

        # Deleting model 'Task'
        db.delete_table('charmming_task')

        # Deleting model 'TaskParameter'
        db.delete_table('charmming_taskparameter')

        # Deleting model 'TaskParameterOption'
        db.delete_table('charmming_taskparameteroption')

        # Deleting model 'Job'
        db.delete_table('charmming_job')

        # Deleting model 'JobTaskParameter'
        db.delete_table('charmming_jobtaskparameter')

        # Deleting model 'Lesson'
        db.delete_table('charmming_lesson')

        # Deleting model 'LessonStep'
        db.delete_table('charmming_lessonstep')


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
            'task_parameters': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['charmming.TaskParameter']"}),
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