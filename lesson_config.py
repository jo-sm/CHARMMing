# This program will be rewritten every time lesson-maker_tpl2.py is run
# Typically it adds new lesson just created by running lesson-maker_tpl2.py 

import lesson1
import lesson2
import lesson3
import lesson4
import lesson5

# Lesson number list
lesson_num_lis = ['1','2','3','4','5']

# The file type lessons need to upload
file_type = {
    '1':'CRD',
    '2':'PDB',
    '3':'CUSTOM/SEQ',
    '4':'CRD',
    '99':'CRD',
    '98':'CRD',
}

# Text-format lesson number
lesson_txt = {
    '1':'Lesson One',
    '2':'Lesson Two',
    '3':'Lesson Three',
    '4':'Lesson Four',
    '5':'Lesson Taco',
    '99':'Lesson 99',
    '98':'Lesson 98'
}

# Lesson description
lesson_desc = {
    '1':'Lesson One will familiarize the user with basic CHARMM commands and basic methods used in computational chemistry.',
    '2':'Lesson Two will use a protein to demonstrate how a user can conduct computational chemistry.',
    '3':'Lesson Three will make the user create their amino acid, replicated 1YJP then running dynamics on it.',
    '4':'Put lesson description here.',
}
