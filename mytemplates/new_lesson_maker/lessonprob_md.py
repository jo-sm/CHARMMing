            if mdt.nstep != {{task.0.nstep}}:
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1}},severity=2,description='NSTEP was not set to {{task.0.nstep}}')
                lessonprob.save()
                return False
#            if mdt.gbmv != {{task.0.gbmv}}:
#                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1}},severity=2,description='gbmv Implicit Solvent was not set to {{task.0.gbmv}}')
#                lessonprob.save()
#                return False
            if mdt.scpism != {{task.0.scpism}}:
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1}},severity=2,description='SCPISM Implicit Solvent was not set to {{task.0.scpism}}')
                lessonprob.save()
                return False
            {%if task.0.make_movie%}blankme
            {%comment%}blankme
            we only check if we want to make a movie, if not
            it doesn't hurt to make one
            {%endcomment%}blankme
            if not mdt.make_movie:
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1}},severity=2,description='make_movie was not set to {{task.0.make_movie}}')
                lessonprob.save()
                return False
            {%endif%}blankme
