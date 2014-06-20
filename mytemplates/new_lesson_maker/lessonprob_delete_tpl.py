        try:
            LessonProblem.objects.get(lesson_type='lesson{{lesson.name}}',lesson_id=self.id).delete()
        except:
            pass
