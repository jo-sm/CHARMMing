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
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User,Group
from wiki.models import *
from pdbinfo.models import *
import os
import re
import smtplib
from django.core.mail import mail_admins
from BeautifulSoup import BeautifulSoup, Comment
#The general Wiki will be the generic wiki 
def loadGeneric(request,boardname):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    topics = []
    topics = getRankedTopics(boardname)
    for item in topics:
        item.title = item.title.replace('\\','')
        item.message = item.message.replace('\\','')
    if request.user.is_superuser:
        superuser = True
    else:
        superuser = False
    return render_to_response("html/wikipages/generic.html", {'topics':topics,'boardname':boardname,'superuser':superuser})

#Used in rating wiki topics
def rateTopic(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    checkWikiRequestData(request)
    board_obj = board.objects.filter(name=request.POST['board_name'])[0]
    rtopic = topic.objects.filter(title = re.escape(request.POST['topic_title']))[0]
    user_rating = request.POST['user_rating']
    rtopic.allratings += int(user_rating)
    rtopic.numratings += 1
    rtopic.rating = round(rtopic.allratings/rtopic.numratings)
    rtopic.save()
    return HttpResponse(str(rtopic.rating))

#checks to make sure a topic's title is unique to the board
#if it is not unique, a number is attached to
#the topic title
def checkUniqueName(topic_title,board_name):
    board_obj = board.objects.filter(name=board_name)[0]
    topic_title = re.escape(topic_title)
    try:
        unique = topic.objects.filter(title=topic_title,board=board_obj)[0]
        i = 2
        topic_title = topic_title.replace('\\','')
        tmptopic_title = topic_title + '(' + str(i) + ')'
        tmptopic_title = re.escape(tmptopic_title)
    except:
        unique = None
    while unique != None:
        i += 1
	try:
            unique = topic.objects.filter(title=tmptopic_title,board=board_obj)[0]
	    tmptopic_title = topic_title + '(' + str(i) + ')'
	    tmptopic_title = re.escape(tmptopic_title)
	except:
	    topic_title = tmptopic_title
	    unique = None
    return topic_title

#Creates a topic
def createTopic(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    checkWikiRequestData(request)
    new_topic = topic()
    if request.POST.has_key('topic_title'):
        if request.POST['topic_title'].strip() == '':
            new_topic.title = checkUniqueName('Untitled',request.POST['board_name'])
        else:
            new_topic.title = checkUniqueName(sanitize_html(request.POST['topic_title']),request.POST['board_name'])
    else:
        # see if this works to not clutter up the wiki
        return False

    if request.POST.has_key('message_post'):
        if request.POST['message_post'].strip() == '':
            new_topic.message = "Be the first to respond to this topic."
        else:
            tmpmessage = ''
            for line in sanitize_html(request.POST['message_post']):
                if line in ['\n']:
                    tmpmessage += '<br>'
                else:
                    tmpmessage += line
            new_topic.message = re.escape(sanitize_html(tmpmessage))
    else:
        new_topic.message = "Be the first to respond to this topic."
    board_obj = board.objects.filter(name = request.POST['board_name'])[0]
    new_topic.board = board_obj
    new_topic.rating = 0
    new_topic.numratings = 0
    new_topic.allratings = 0
    new_topic.save()
    return loadGeneric(request,request.POST['board_name'])

#delete topic
def deleteTopic(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    if not request.user.is_superuser:
        return HttpResponse('Access Denied!')
    topic_title = re.escape(request.POST['topictitle'])
    board_name = request.POST['boardname']
    board_obj = board.objects.filter(name = board_name)[0]
    current_topic = topic.objects.filter(title=topic_title,board=board_obj)[0]
    current_topic.delete()
    return HttpResponse('deleted!')

#Uses inline editing, retrieves the data, and then edits it
def editMessage(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    checkWikiRequestData(request)
    #If the user wants to edit the topic title then the request.GET['content'] is the
    #actual topic title whereas edittitle is just the <span> id
    board_name = request.GET['boardname']
    board_obj = board.objects.filter(name = board_name)[0]
    topictitle = request.GET['fieldname']
    realtopictitle = topictitle.split('_',1)
    current_topic = topic.objects.filter(title=re.escape(realtopictitle[1]),board=board_obj)[0]
    if realtopictitle[0] in ['edittitle']:
        current_topic.title = re.escape(sanitize_html((request.GET['content'])))
        current_topic.save()
        return HttpResponse(current_topic.title.replace('\\',''))
    else:
        if request.GET['content'] == '':
            current_topic.message = re.escape("Be the first to respond to this topic.")
	else:
	    newline = re.compile('\n')
	    message = ''
	    for line in request.GET['content']:
	        if newline.search(line):
		    message += '<br>'
		    message += line
		else:
		    message += line
            message = sanitize_html(message)
            current_topic.message = re.escape(message)
        current_topic.save()
        return HttpResponse(current_topic.message.replace('\\',''))

def getRankedTopics(board_name):
    board_obj = board.objects.filter(name = board_name)[0]
    topics = []
    topicrank5 = topic.objects.filter(board = board_obj,rating = 5)
    for item in topicrank5:
        topics.append(item)
    topicrank4 = topic.objects.filter(board = board_obj,rating = 4)
    for item in topicrank4:
        topics.append(item)
    topicrank3 = topic.objects.filter(board = board_obj,rating = 3)
    for item in topicrank3:
        topics.append(item)
    topicrank2 = topic.objects.filter(board = board_obj,rating = 2)
    for item in topicrank2:
        topics.append(item)
    topicrank1 = topic.objects.filter(board = board_obj,rating = 1)
    for item in topicrank1:
        topics.append(item)
    topicrank0 = topic.objects.filter(board = board_obj,rating = 0)
    for item in topicrank0:
        topics.append(item)
    return topics

#check request data for malicious code
def checkWikiRequestData(request):
    for parameter in request.POST:
        checkForWikiMaliciousCode(request.POST[parameter],request)
    for parameter in request.GET:
        checkForWikiMaliciousCode(request.GET[parameter],request)


def sanitize_html(value):
    valid_tags = 'i b u br'.split()
    valid_attrs = ''.split()
    soup = BeautifulSoup(value)
    for comment in soup.findAll(
        text=lambda text: isinstance(text, Comment)):
        comment.extract()
    for tag in soup.findAll(True):
        if tag.name not in valid_tags:
            tag.hidden = True
        tag.attrs = [(attr, val) for attr, val in tag.attrs
                     if attr in valid_attrs]
    return soup.renderContents().decode('utf8').replace('javascript:', '')

def sanitizeHTML(html):
    htmlcode = re.compile('(<[b,i,u]*>)(.*)(</[b,i,u]*>)')
    htmlgroup = htmlcode.match(html)
    validhtmlstart = []
    validhtmlend = []
    validhtml = []
    breakloop = 0
    while not breakloop: 
        if htmlgroup.group(1).lower() in ['<b>','<u>','<i>']:
            validhtmlstart.insert(0,htmlgroup.group(1).strip('<').strip('>'))
        elif htmlgroup.group(3).lower() in ['</b>','</u>','</i>']:
            validhtmlend.insert(0,htmlgroup.group(1).strip('</').strip('>'))
	    if validhtmlend[0] == validhtmlstart[0]:
	        validhtml.append
	        
        htmlgroup = htmlcode.match(htmlgroup.group(2))
	if htmlgroup == None:
	    breakloop = 1
    return validhtml
         

    #checks data for possible malicious code

def checkForWikiMaliciousCode(text,request):
    #syst is how charmm executes system commands
    semicolon = re.compile(';')
    osdot = re.compile('os\.')
    if(semicolon.search(text) or osdot.search(text)):
        msg = "User :" + request.user.username + "\n tried to execute malicous code with:\n " + text + "\n"
        mail_admins('Malcious Code Warning',msg,fail_silently=False)
        sys.exit()
        return "Dangerous Data! Attempt has been logged."
    return text
