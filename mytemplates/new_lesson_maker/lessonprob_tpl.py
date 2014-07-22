        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson{{lesson.name}}',lesson_id=self.id)[0]
        except:
            lessonprob = None
        if lessonprob:
            self.curStep = str(float(self.curStep - self.curStep % 1))
            self.save()
            return False
