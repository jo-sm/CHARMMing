#
#                            PUBLIC DOMAIN NOTICE
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the authors' official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software is freely available
#  to the public for use.  There is no restriction on its use or
#  reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, NIH and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. NIH, NHLBI, and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any
#  particular purpose.
from django.contrib import auth
from django.contrib import messages
from django.contrib.messages import get_messages
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User,Group
from structure.models import Structure
from account.models import *
import django.forms
import os, re, smtplib, random, string
import lessonaux, charmming_config
from django.core.mail import mail_admins, send_mail

#If the user visits the main page it will redirect them to
#the main app or a login page
def redirectBasedOnLogin(request):
    if(request.user.is_authenticated()):
        return HttpResponseRedirect('/charmming/accounts/profile/')
    return HttpResponseRedirect('/charmming/login/')  

#resets a user's password
def resetUserPassword(request):
    if(request.POST):
        #email-given and user_id are the the variables the user has submitted
        email_given = request.POST['email']
        user_id = request.POST['username']
        #User should not be able to reset the password of the demo account
        if(user_id == 'demo'):
            return HttpResponse("Did you really think we wouldn't test for that?") #I leave this one as is - it's better that way. ~VS
        password = ''
        try:
            #Attempt to get the user_profile based on the user supplied info
            user_profile = User.objects.get(email__iexact = email_given, username=user_id)
        except:
            messages.error(request, "Sorry, no user could be located with the given information.")
            return HttpResponseRedirect("/charmming/accounts/resetpassword/")
#            return HttpResponse("Sorry, no user could be located with the given information")
        if(user_profile.is_superuser or user_profile.is_staff):
            messages.error(request, "Staff and Admins cannot get their password reset. Whoops!")
            return HttpResponseRedirect("/charmming/accounts/resetpassword/")
#            return HttpResponse("Staff and Admins cannot get their password reset. Whoops!")
        chars = string.ascii_letters + string.digits
        for x in range(8):
            password += random.choice(chars)
        mail_message = "Your password for CHARMMing.org has been reset. Your new password is: " + password 
        mail_message += ". You can change the password once you log in. If you did not "
        mail_message += "request this password reset or there are any other "
        mail_message += "problems please E-mail btamiller@gmail.com. Please do not reply to this E-mail.\n\n"
        mail_message += "Thanks for using CHARMMing,\n CHARMMing Development Team"
        send_mail('CHARMMing Password Reset',mail_message,'PasswordReset@charmming.org',[user_profile.email],\
                 fail_silently = False)
        user_profile.set_password(password)
        user_profile.save()
        return HttpResponse("Password changed! You should receive an e-mail shortly.")
    ##TODO: Change this to a dialog box when we have message.tags
    else:
        return render_to_response('html/resetpassword.html')

#will allow a user to change their password
def changeUserPassword(request):
    User = request.user
    #check to see if user has submitted any data
    if(request.POST):
        old_user_password = request.POST['oldpassword']
        new_user_password = request.POST['newpassword']
        new_user_password_confirm = request.POST['newpasswordconfirm']
        #make sure the user supplied the correct old password first
        if(User.check_password(old_user_password)):
            #ensure the new password and the new confirmation password match
            if(new_user_password == new_user_password_confirm):
                User.set_password(new_user_password)
                User.save()
                return HttpResponse("Password changed!")
            else:
                messages.error(request, "The two new passwords you entered did not match! Please go back and re-enter your password.")
                return HttpResponseRedirect("/charmming/accounts/changepassword/")
#                return HttpResponse("The two new passwords you entered did not match! Please go back and re-enter your password.")
        else:
            messages.error(request, "Your old password was not correct. Please go back and re-enter your password.")
            return HttpResponseRedirect("/charmming/accounts/changepassword/")
#            return HttpResponse("Your old password was not correct. Please go back and re-enter your password.")
    #if not, see if they are alteast logged in. If true, direct them to
    #the password change form
    elif(User.is_authenticated()):
        return render_to_response('html/changepassword.html')
    #if they are not logged in, return them to the front page
    else:
        return render_to_response('html/frontpage.html')
    
#pre: Needs a valid teacherProfile
#This will get the list of classrooms based on a teacher
def getClassList(teacher):
    return classroom.objects.filter(teacher=teacher)

#When a students inputs a teacher, this will get the list of classes associated with the teacher
def refreshClassDiv(request,teacher_name):
    slash = re.compile('/')
    teacher_name = slash.sub('',teacher_name)
    try:
        #Stores user object of the student's teacher
        students_teacher_user = User.objects.filter(username=teacher_name)[0]
        #Gets the teacher Profile of the teacher
        students_teacher = teacherProfile.objects.filter(teacher = students_teacher_user)[0]
	class_list = classroom.objects.filter(teacher = students_teacher)
	return render_to_response('html/classlistdiv.html',{'class_list':class_list})
    except:
        return HttpResponse("Sorry, no such teacher exists. Name is case-sensitive.")
 

#registers a student
def registerStudent(request,user):
    #added temporarily to get rid of teacher/student stuff
    ##This may be completely messed up - the indentation was horrible according to Vim
    # So I took actuon and added error messages. If something is broken, that's probably why!
    #So restore the old copy of this file from production. ~VS
    request.user = user
    if user.is_authenticated():
        if(request.POST):
            #the explanation for why the teacher will use the interface
            #is saved in the directory ~/teachers/teacher_user_name/
            new_student = studentProfile(student=request.user)
            new_student.institute = request.POST['institute']
            if new_student.institute.strip() == '':
                user.delete()
                messages.error(request, 'Institute name is required.')
                return HttpResponseRedirect("/charmming/accounts/register/student/")
#                return HttpResponse('Institute name is required.')

            if request.POST['charmmlicense'] == "No":
               new_student.charmm_license = 0
            else:
               new_student.charmm_license = 1
            #Stores user object of the student's teacher
            #students_teacher_user = User.objects.filter(username=request.POST['teacher_name'])[0]
            students_teacher_user = User.objects.filter(id=1)[0]
            #Gets the teacher Profile of the teacher
            #students_teacher = teacherProfile.objects.filter(teacher = students_teacher_user)[0]
            #new_student.teacher = students_teacher
            #Now to set the students classroom
            #teacher.teacher.username is the teacher's username to log into the site
            classname = 'General Public'
            #commented out because getting rid of teacher/student thing
            #classname = request.POST[new_student.teacher.teacher.username]
            #new_student.classroom = classroom.objects.filter(name=classname)[0]
            new_student.student.username = new_student.student.username.lower()
            new_student.save()
            #New student needs a folder
            os.mkdir("%s/%s" % (charmming_config.user_home,new_student.student.username))
            os.chmod("%s/%s" % (charmming_config.user_home,new_student.student.username), 0775)

            preapprove = Group.objects.get(name='preapprove')
            request.user.groups.add(preapprove)
            mail_message = """
New student awaiting approval.
Username:""" + new_student.student.username + """
Name:""" + new_student.student.first_name + " " + new_student.student.last_name + """
E-mail: """ + new_student.student.email + """
Institute:""" + new_student.institute + """
CHARMM License?:""" + request.POST['charmmlicense']
            mail_admins('new student',mail_message,fail_silently=False)
            return HttpResponse('Your account is now awaiting approval. It may take 1-2 days.')
        #return render_to_response('registration/registerstudent.html')
    return HttpResponse('Please go through the registration page first.')

#Will assign user a 'preteacher' group, where the account wont be able to
#do anything until an admin approves their account
def registerTeacher(request):
    if request.user.is_authenticated():
        if(request.POST):
            #the explanation for why the teacher will use the interface
            #is saved in the directory ~/teachers/teacher_user_name/
            new_teacher = teacherProfile(teacher=request.user)
            new_teacher.institute = request.POST['institute']
            new_teacher.phone_number = request.POST['phone_number']
            new_teacher.save()
            try:
                os.mkdir('%s/teachers/%s' % (charmming_config.user_home,new_teacher.teacher.username))
            except:
                pass
            reason_handle = open('%s/teachers/%s/explanation.txt' % (charmming_config.user_home,new_teacher.teacher.username),'w')
            reason_handle.write(request.POST['explanation'])
            reason_handle.close()
            mail_message = """
New Teacher waiting for approval. 
Name:""" + new_teacher.teacher.first_name + " " + new_teacher.teacher.last_name + """
Institute:""" + new_teacher.institute + """
Phone #:""" + new_teacher.phone_number + """
Explanation:""" + request.POST['explanation'] 
            mail_admins('new teacher',mail_message,fail_silently=False)
            return HttpResponse('Your account is now awaiting approval.')
            preapprove = Group.objects.get(name='preapprove')
            request.user.groups.add(preapprove)
            return render_to_response('registration/registerteacher.html')
        return HttpResponse('Please go through the registration page first.')
    else:
        return render_to_response("html/loggedout.html")
    #If not authenticated this should just exit, right? Or what happens?

def register(request):
    if request.user.is_authenticated():
        return HttpResponseRedirect('/charmming/accounts/profile/')
    if request.method == 'POST':
        # create bound form
        form = UserCreationForm(request.POST)
        if form.is_valid():

            username = form.cleaned_data['username'].lower()
            password = form.cleaned_data['password1']
            first_name = request.POST['first_name']
            last_name = request.POST['last_name']
            if not first_name.strip() or not last_name.strip():
                return HttpResponse('A complete name is required.')
            email = request.POST['email']
            if email.strip() == '': 
                return HttpResponse('An E-mail address is required.')
            user = User.objects.create_user(username=username,email=email,password=password)
            user.first_name = first_name
            user.last_name = last_name
            user.save()

            user = auth.authenticate(username=username, password=password)
            auth.login(request, user)

            return registerStudent(request,user)
    else:
        # create unbound form
        form = UserCreationForm()

    return render_to_response("registration/register.html", {'form' : form})


def login(request):
    try:
        username = request.POST['username'].lower()
        password = request.POST['password']
    except:
        return HttpResponseRedirect("/charmming/")

    # sanity check inputs
    if not username or not password:
        return HttpResponseRedirect("/charmming/")
    if ";" in username or "<" in username or ">" in username or "<" in password or ">" in password or ";" in password:
        return HttpResponseRedirect("/charmming/")

    user = auth.authenticate(username=username, password=password)
    if user and user.is_active:
        # Correct password, and the user is marked "active"
        auth.login(request, user)
        #check to see if there is a lesson for the lesson status bar
        try:
            file = Structure.objects.filter(owner=request.user,selected='y')[0]
            lessontype = file.lesson_type
        except:
            lessontype = None
        # Redirect to a success page.
        return HttpResponseRedirect("/charmming/html/skeleton/")
    else:
        # Show an error page
        return render_to_response("html/frontpage.html",{'login_error':'Your username/password did not match'})

def validate(request):
    if request.user.is_authenticated():
        return HttpResponseRedirect("/charmming/html/skeleton/")
    return HttpResponseRedirect("/charmming/")

def skeletonDowntime(request):
    if request.user.is_authenticated():
        return render_to_response("html/downtime.html")
    else:
        return HttpResponse("You can log in now.")

def isUserTrustworthy(user):
    groupList = user.groups.all()
    if(groupList.filter(name='trusted') or user.is_superuser):
        return True
    else:
        return False

def checkPermissions(request):
    # returns lesson, drug gesign
    if request.user.is_authenticated():
        return True, True
    else:
        return False, False

def skeleton(request):
    if request.user.is_authenticated():
        # get a list of all the groups that the current user is in
        groupList = request.user.groups.all()
        #Checks to see if the user still needs to be preapproved and makes sure
	#an admin doesn't fall into this category as in admin is in all groups
        if(groupList.filter(name='preapprove') and not request.user.is_superuser):
            return render_to_response('registration/waitingapproval.html')
        #check to see if there is a lesson for the lesson status bar
        try:
            file = Structure.objects.filter(owner=request.user,selected='y')[0]
            lessontype = file.lesson_type
        except:
            lessontype = None
    else:
        return HttpResponse("Login please")
    #If it is a teacher, send teacher to the skeleton interface
    if request.user.is_superuser:
        superuser = "true"
        trusted = "true"
        return render_to_response("html/skeleton.html", {'superuser':superuser,'trusted':trusted,'lessontype':lessontype, 'lessonuser': 1})
    if groupList.filter(name='trusted'):
        trusted = 1
    else:
        trusted = 0

    if groupList.filter(name='lesson'):
        lessonuser = 1
    else:
        lessonuser = 0
        #If they are only a student load the below page
    return render_to_response("html/skeleton.html",{'lessontype':lessontype, 'trusted': trusted, 'lessonuser': lessonuser})

def logout(request):
    auth.logout(request)
    # Redirect to a success page.
    return HttpResponseRedirect("/charmming", {'logout':1})


def loadFrontPage(request):
    if(request.user.is_authenticated()):
        return HttpResponseRedirect('/charmming/accounts/profile/')

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/frontpage.html', {'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

def loadAboutPage(request):
    lesson_ok, dd_ok = checkPermissions(request)
    messages = get_messages(request)
    return render_to_response('html/about.html', {'lesson_ok': lesson_ok, 'dd_ok': dd_ok, 'messages': messages})
