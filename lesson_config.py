# This program will be rewritten every time lesson-maker_tpl2.py is run
# Typically it adds new lesson just created by running lesson-maker_tpl2.py 

import lesson1
import lesson2
import lesson3
import lesson4
import lesson5
import lesson6

# Lesson number list
lesson_num_lis = ['1','2','3','4','5','6']
#lesson_num_lis = ['1','2','3','4','5']

# The file type lessons need to upload
file_type = {
    '1':'CRD',
    '2':'PDB',
    '3':'CUSTOM/SEQ',
    '4':'CRD',
    '5':'CRD',
    '6':'PDB',
}

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
