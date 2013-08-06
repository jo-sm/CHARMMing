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
from django import forms
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from django.contrib.auth.models import User,Group
from django.core.mail import send_mail
from account.models import studentProfile
from statistics.models import DataPoint

def showUnapproved(request):
    #unapproved_list is the list of unapproved users
    #user_profile_list keeps track of the User objects that are selected for approval
    #approved_list keeps track of the users approved
    user_profile_list = []
    unapproved_list = []
    approved_studentProfile_list = []
    disapproved_studentProfile_list = []
    #first checks to see if the user submitted anything
    #if not, display the page to select users
    if(request.POST):
        unapproved_list = User.objects.filter(groups__name='preapprove')
        #obtains student and preapprove group objects
        student = Group.objects.get(name='student')
        preapprove = Group.objects.get(name='preapprove')
        #cycles through each unapproved user testing
        for user_profile in unapproved_list:
            mail_message = ""
            try:
               #check if an E-mail should be sent
                try:
                    if(request.POST['squelch_' + user_profile.username]):
                        send_email = 0
                except:
                    send_email = 1
                #checks to see if the user's profile was checked
                #the name of each checkbox on the front-end webpage is just the user's username
                if(request.POST['approve_or_disapprove_' + user_profile.username] == 'approve_' + user_profile.username):
                    #because users can be in multiple groups, remove them from preapprove first then add to students
                    dict = {}
                    user_profile.groups.remove(preapprove)
                    user_profile.groups.add(student)
                    user_profile.save()
                    student_obj = studentProfile.objects.filter(student = user_profile)[0]
                    #creates a custom mail message
                    #but first check and see if the checkbox NOT to send an E-mail to the user was marked
                    if(send_email == 1):
                        mail_message += "Hi " + user_profile.first_name + ",\n\n"
                        mail_message += "Your account has been approved for using www.CHARMMing.org. If you have any problems "
                        mail_message += "logging in, please E-mail btamiller@gmail.com. Please do not reply to this E-mail\n\n"
                        if(request.POST['comments_' + user_profile.username] != ""):
                            mail_message += "The following comments were appended to your approval:\n"
                            mail_message += request.POST['comments_' + user_profile.username] + "\n\n"
                            dict['comments'] = request.POST['comments_' + user_profile.username] 
                        mail_message += "Thanks for using CHARMMing,\n CHARMMing Development Team"
                        send_mail('CHARMMing Approval',mail_message,'CHARMMingApproval@learn.lobos.nih.gov',[user_profile.email],fail_silently = False)
                    else:
                        dict['comments'] = "No E-mail was sent."
                    dict['username'] = user_profile.username;
                    dict['email'] = user_profile.email
                    dict['institute'] = student_obj.institute
                    dict['charmm_license'] = student_obj.charmm_license
                    dict['date_joined'] = user_profile.date_joined
                    approved_studentProfile_list.append(dict)
                elif(request.POST['approve_or_disapprove_' + user_profile.username] == 'disapprove_' + user_profile.username):
                    #Send user an E-mail, then delete em'! 
                    dict = {}
                    student_obj = studentProfile.objects.filter(student = user_profile)[0]
                    if(send_email == 1):
                        mail_message += "Hi " + user_profile.first_name + ",\n\n"
                        mail_message += "Your account has NOT been approved for using www.CHARMMing.org. If you have any questions "
                        mail_message += "about why you were not approved, please E-mail btamiller@gmail.com. Please do not reply to this E-mail\n\n"
                        if(request.POST['comments_' + user_profile.username] != ""):
                            mail_message += "The following comments were appended to your disapproval:\n"
                            mail_message += request.POST['comments_' + user_profile.username] + "\n\n"
                            dict['comments'] = request.POST['comments_' + user_profile.username]
                        mail_message += "Thanks for using CHARMMing,\n CHARMMing Development Team"
                        send_mail('CHARMMing Disapproval',mail_message,'CHARMMingDisapproval@learn.lobos.nih.gov',[user_profile.email],\
                                      fail_silently = False)
                    else:
                        dict['comments'] = "No E-mail was sent."
                    dict['username'] = user_profile.username;
                    dict['email'] = user_profile.email
                    dict['institute'] = student_obj.institute
                    dict['charmm_license'] = student_obj.charmm_license
                    dict['date_joined'] = user_profile.date_joined 
                    disapproved_studentProfile_list.append(dict)
                    studentProfile.objects.filter(student = user_profile)[0].delete()
                    user_profile.delete()

            except:
                pass
    #get new unapproved list
    unapproved_list = User.objects.filter(groups__name='preapprove')
    unapproved_studentProfile_list = []
    for user_profile in unapproved_list:
        unapproved_studentProfile_list.append(studentProfile.objects.filter(student = user_profile)[0])
    return render_to_response('admin/unapproveduser.html', {'unapproved_list':  unapproved_studentProfile_list, 'approved_list': approved_studentProfile_list,'disapproved_list':disapproved_studentProfile_list})

def showStats(request):
    tdict = {}
    datapoints = DataPoint.objects.all() #Get everything since we don't care about any filtering criteria
    task_actions = {}
    task_action_list = []
    total_tasks = len(datapoints)
    successful_tasks = 0
    failed_tasks = 0
    users = set([])
    for point in datapoints:
        if point.success:
            successful_tasks += 1
        else:
            failed_tasks += 1
        if point.user not in users:
            users = users.union(set([point.user]))
        if point.task_action not in task_actions.keys():
            task_actions[point.task_action] = 1
        else:
            task_actions[point.task_action] += 1
    total_users = len(users)
    #This is silly but it's also the easiest way to get it right. I tried to handle it with a list and the lookups get silly and mess my code up, at least this works
    for key, val in task_actions.iteritems():
        task_action_list.append((key, val))
    task_action_list = sorted(task_action_list,key=lambda x:x[0]) #Make it into a sorted list for easy display later
    tdict['task_list'] = task_action_list
    tdict['successful_tasks'] = successful_tasks
    tdict['failed_tasks'] = failed_tasks
    tdict['total_users'] = total_users
    tdict['total_tasks'] = total_tasks
    return render_to_response('admin/statistics.html', tdict)
