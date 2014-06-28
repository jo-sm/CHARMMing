# This program will be rewritten every time create_lesson.py is run
# Typically it adds new lesson just created by running create_lesson.py 
{%comment%}blankme
we make a list of dictionaries to pass into here, since I don't trust
django to process objects I randomly create correctly
{%endcomment%}blankme
{%for lesson in lessons%} blankme
import lesson{{lesson.name}}
{%endfor%} blankme

# Lesson number list
name_list = [
{%for lesson in lessons%}blankme
'{{lesson.name}}'{%if not forloop.last%},{%endif%}
{%endfor%}blankme
    ]
#name_lis = ['1','2','3','4','5']

# The file type lessons need to upload
{%comment%}blankme
When we do this, we will fetch the old lesson file types from lesson_config
then update them.
{%endcomment%}blankme
file_type = {
{%for lesson in lessons%}blankme
    '{{lesson.name}}':'{{lesson.file_type}}{%if not forloop.last%},{%endif%}
{%endfor%}blankme
}

{%comment%}blankme
We don't actually need this block of code anymore, but in case we DO,
we leave this in here. THis is a sample!
# Text-format lesson number
lesson_txt = {
    '1':'Lesson One',
    '2':'Lesson Two',
    '3':'Lesson Three',
    '4':'Lesson Four',
    '5':'Lesson Five',
    '6':'Lesson Six',
}

# Lesson description
lesson_desc = {
    '1':'Lesson One will familiarize the user with basic CHARMM commands and basic methods used in computational chemistry.',
    '2':'Lesson Two will use a protein to demonstrate how a user can conduct computational chemistry.',
    '3':'Lesson Three will make the user create their amino acid, replicated 1YJP then running dynamics on it.',
    '4':'Lesson Four introduces custom RTFs and QM/MM.',
    '5':'Lesson Four shows the user how to build a coarse-grained Go model.',
    '6':'Lesson Six will show the user how to calculate the reduction potential of a [4Fe-4S]-containing protein.',
}
{%endcomment%}blankme

#Lesson title - used to simplify skeleton templating
{%comment%}blankme
We do actually need this one. It genericizes skeleton.html.
{%endcomment%}blankme
title = {
{%for lesson in lessons%}blankme
    '{{lesson.name}}':'{{lesson.title}}'{%if not forloop.last%},{%endif%}
{%endfor%}blankme
}
