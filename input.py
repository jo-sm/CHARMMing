import sys, re
from django.core.mail import mail_admins

#checks data for possible malicious code
def checkForMaliciousCode(text,request):
    #syst is how charmm executes system commands
    syst = re.compile('syst')
    semicolon = re.compile(';')
    osdot = re.compile('os\.')
    if(syst.search(text) or semicolon.search(text) or osdot.search(text)):
        msg = "User :" + request.user.username + "\n tried to execute malicous code with:\n " + text + "\n"
        mail_admins('Malcious Code Warning',msg,fail_silently=False)
        sys.exit()
        return "Dangerous Data! Attempt has been logged."
    return text

