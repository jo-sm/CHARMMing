#!/bin/bash 
# 
# A new, and hopefully improved, installer for the 
# CHARMM interface and graphics (CHARMMing)
# 
# H. Lee Woodcock 8/19     Version 0.1 
#
# To install set the paths and passwords in the 
# "set_install_config" section (directly below)
# and run the installer.... 
#
# ./installer -install        (to install) 
# ./installer -uninstall      (to uninstall) 
#


function set_install_config {
  charmming_path="/var/www/charmming"       # Location of the main CHARMMing installation 
  charmming_utils="/usr/local/charmming"    # Location of CHARMMing util installation 
  charmm="/usr/local/charmm/charmm"         # Location of the CHARMM binary  
  schedd_user="schedd"                      # CHARMMing scheduler (DO NOT MODIFY)!!!  
  charmming_admin="admin"                   # The username for the charmming administrator 
  charmming_passwd="qwerty"                 # The administrator password for charmming 
  database_passwd="qwerty"                  # The mysql password for the charmming database 
  database_root_passwd="qwerty"             # The root password for the MySQL database
  apache_group="www-data"                   # The group that Apache runs under
}

function dependencies { 

echo -e "\nDo you want to use apt to install packages?" 
read answer

if [[ "$answer" == "Y" || "$answer" == "y" || "$answer" = "" || "$answer" = "yes" || "$answer" = "Yes" || "$answer" = "YES" ]];
then
  sudo apt-get install bc subversion python-django python-openbabel openbabel torque-server torque-scheduler torque-client torque-common torque-mom \
               apache2 mysql-server mysql-client mysql-common libapache2-mod-python python-mysqldb python-numpy python-scikits-learn python-rdkit \
               rdkit-data rdkit-doc csh tcsh openmpi-common openmpi-bin libopenmpi-dev
else
  echo    "-----------------------"
  echo -e "Ok, I will skip this..." 
  echo -e "-----------------------\n"
fi
}

function download_charmming {
  mkdir -p /tmp/charmming/
  cd /tmp/charmming/  

  echo -e "\nDo you want to check out and install new versions of charmming and charmminglib?"
  read answer
  if [[ "$answer" == "Y" || "$answer" == "y" || "$answer" = "" || "$answer" = "yes" || "$answer" = "Yes" || "$answer" = "YES" ]];
  then
    if [ -d "charmming" ]; then
	    mv charmming charmming-'date +%Y-%m-%d-%a-%I.%M%p'
    fi 
    
    if [ -d "charmminglib" ]; then
	    mv charmminglib charmminglib-`date +%Y-%m-%d-%a-%I.%M%p`
    fi
    
# checkout charmming and charmminglib from svn 
    svn co http://charmming.googlecode.com/svn/branches/CHARMMING_0_10 charmming
    svn co http://charmming.googlecode.com/svn/charmminglib charmminglib
    
# install charmminglib 
    cd charmminglib 
    sudo python setup.py install
    cd ../
   
# make sure visualization stuff is setup 
    cd /tmp/charmming/charmming/mytemplates/js/
    unzip glmol.zip 
    unzip jsmol-13.1.13.zip
    tmpdir=`dirname  $charmming_path`
    static="$tmpdir/charmming-static"
    sudo mkdir "$static" 
    sudo mv jsmol "$static"
    cd /tmp/charmming/
 
#move charmming and charmming utils to final location.. 
    sudo mv charmming $charmming_path
    sudo cp -a $charmming_path/charmming-private $charmming_utils
  else
    echo    "-----------------------"
    echo -e "Ok, I will skip this..." 
    echo -e "-----------------------\n"
  fi
}

function mysql_config {
  echo -e "\nDo you want to create and setup the default charmming database?" 
  read answer
  if [[ "$answer" == "Y" || "$answer" == "y" || "$answer" = "" || "$answer" = "yes" || "$answer" = "Yes" || "$answer" = "YES" ]];
  then
    cd /tmp/charmming/
    echo "create database charmming;" > create_db.sql
    echo "grant usage on *.* to 'charmming'@'localhost' identified by '$database_passwd';" >> create_db.sql
    echo "grant all on charmming.* to 'charmming'@'localhost';" >> create_db.sql
    echo "flush privileges;" >> create_db.sql

# call mysql to create database and set permissions...
    echo " Creating charmming database:" 
    mysql -u root --password="$database_root_passwd" < /tmp/charmming/create_db.sql

    echo "CREATE TABLE \`job_scheduler\` ("                       >  scheduler.sql
    echo "  \`sched_id\` varchar(64) default NULL,"               >> scheduler.sql
    echo "  \`state\` tinyint(3) unsigned NOT NULL default '0',"  >> scheduler.sql
    echo "  \`started\` datetime default NULL,"                   >> scheduler.sql
    echo "  \`ended\` datetime default NULL,"                     >> scheduler.sql
    echo "  \`userid\` int(11) NOT NULL default '0',"             >> scheduler.sql
    echo "  \`id\` bigint(20) NOT NULL auto_increment,"           >> scheduler.sql
    echo "  \`queued\` datetime default NULL,"                    >> scheduler.sql
    echo "  \`exe\` text,"                                        >> scheduler.sql
    echo "  \`script\` text,"                                     >> scheduler.sql
    echo "  \`batchid\` varchar(64),"                             >> scheduler.sql
    echo "  PRIMARY KEY  (\`id\`)"                                >> scheduler.sql
    echo ") ENGINE=MyISAM DEFAULT CHARSET=latin1;"                >> scheduler.sql

# call mysql to create database and set permissions...
    echo -e " Creating charmming tables:\n" 
    mysql -u charmming --password="$database_passwd" -D charmming < /tmp/charmming/scheduler.sql
   
  else
    echo    "-----------------------"
    echo -e "Ok, I will skip this..." 
    echo -e "-----------------------\n"
  fi
}

function schedd_config {
  echo -e "\nDo you want/need to create the schedd user and modify charmming paths?" 
  read answer
  if [[ "$answer" == "Y" || "$answer" == "y" || "$answer" = "" || "$answer" = "yes" || "$answer" = "Yes" || "$answer" = "YES" ]];
  then
    sudo useradd schedd -c "CHARMMing scheduler daemon" -s /bin/bash -m
    sudo usermod -a -G "$apache_group" schedd
    sudo mkdir -p /home/schedd/$charmming_admin
    sudo chown -R schedd.www-data /home/schedd/ 
    sudo chmod -R g+w /home/schedd/ 
    sudo chown -R schedd.schedd $charmming_path/scheduler/
    sudo sed -i "s/myhost/`hostname`/" $charmming_path/scheduler/pbs_sched.py
  
# set debug options to false for production sites... 
    sudo sed -i "s/DEBUG = True/DEBUG = False/"                                                   $charmming_path/settings.py 
    sudo sed -i "s/TEMPLATE_DEBEUG = True/TEMPLATE_DEBEUG = False/"                               $charmming_path/settings.py 
    sudo sed -i "s/DATABASE_PASSWORD = 'charmming'/DATABASE_PASSWORD = '"$charmming_passwd"'/"    $charmming_path/settings.py
    sudo sed -i "s/MEDIA_ROOT = '\/home\/pdb_uploads'/MEDIA_ROOT = '\/home\/schedd'/"             $charmming_path/settings.py
    sudo sed -i "s/MEDIA_URL = '\/charmming\/pdbuploads\/'/MEDIA_URL = '\/charmming\/schedd\/'/"  $charmming_path/settings.py
    sudo ln -s /home/schedd $charmming_path/schedd

# set charmm binary location and misc.... 
    charmming_path2=`echo $charmming_path | sed -e 's:\/:\\\/:g'`
    sudo sed -i "s/\/var\/www\/charmming/"$charmming_path2"/"                   $charmming_path/charmming_config.py
    sudo sed -i "s/\/home\/pdb_uploads/\/home\/schedd/"                         $charmming_path/charmming_config.py
    charmm2=`echo $charmm | sed -e 's:\/:\\\/:g'`
    sudo sed -i "s/\/usr\/local\/charmming\/gfortran-xxlg-qc.one/"$charmm2"/"   $charmming_path/charmming_config.py
  else
    echo    "-----------------------"
    echo -e "Ok, I will skip this..." 
    echo -e "-----------------------\n"
  fi
}

function torque_config {
  echo -e "\nDo you need to setup torque?" 
  read answer
  if [[ "$answer" == "Y" || "$answer" == "y" || "$answer" = "" || "$answer" = "yes" || "$answer" = "Yes" || "$answer" = "YES" ]];
  then
    sudo killall -9 pbs_server pbs_sched pbs_mom
    sudo /etc/init.d/torque-server stop 
    sudo /etc/init.d/torque-scheduler stop 
    sudo /etc/init.d/torque-mom stop 

    cd /tmp/charmming/ 
    hostname > server_name 
    sudo mv server_name /etc/torque/server_name 
    sudo sed -i 's/127.0.1.1/127.0.0.1/' /etc/hosts
  
    a=`cat /proc/cpuinfo | grep processor | tail -1 | tail -c 2 | head -c 1`
    nproc=`echo $a+1 | bc`
  
    echo "`hostname` np=$nproc" > nodes
    sudo mv nodes /var/spool/torque/server_priv/nodes
  
# setup default torque queues options 
    echo "# step 1 -- create queues: basic set-up with a"         >  conf_torque 
    echo "# routing and execution queue."                         >> conf_torque
    echo "create queue entry"                                     >> conf_torque
    echo "create queue charmming"                                 >> conf_torque
    echo "set queue entry queue_type = route"                     >> conf_torque
    echo "set queue entry route_destinations = charmming"         >> conf_torque
    echo "set queue entry enabled = True"                         >> conf_torque
    echo "set queue entry started = True"                         >> conf_torque
    echo "set queue charmming queue_type = exec"                  >> conf_torque
    echo "set queue charmming from_route_only = True"             >> conf_torque
    echo "set queue charmming resources_default.nodes = 1:ppn=1"  >> conf_torque
    echo "set queue charmming enabled = True"                     >> conf_torque
    echo "set queue charmming started = True"                     >> conf_torque
    echo ""                                                       >> conf_torque
    echo "# Step 2 -- configure server and get things rolling"    >> conf_torque
    echo "set server scheduler_iteration = 60"                    >> conf_torque
    echo "set server default_queue = entry"                       >> conf_torque
    echo "set server scheduling = True"                           >> conf_torque
# create default mom config file 
    echo "\$logevent       0x1ff"                                 >  node_config
    echo "\$clienthost     `hostname`"                            >> node_config
    echo "\$usecp          *:/       /"                           >> node_config
    echo "\$ideal_load      2.25"                                 >> node_config
    echo "\$max_load        3.95"                                 >> node_config
    sudo cp node_config /var/spool/torque/mom_priv/config

# start sever and scheduler 
    sudo pbs_server -t create
    sudo pbs_sched

# create default torque queue and nodes 
    sudo qmgr < conf_torque
    sudo qmgr -c "create node `hostname`"
    sudo qmgr -c "set node `hostname` np = $nproc"
    sudo qmgr -c "set node `hostname` properties = exec"

# start mom 
    sudo pbs_mom
  else
    echo    "-----------------------"
    echo -e "Ok, I will skip this..." 
    echo -e "-----------------------\n"
  fi
}

function config_apache {
#/etc/apache2/sites-enabled/000-default
#  echo "   <Directory "/var/www/charmming">"                    >  conf_apache
#  echo "      SetHandler python-program"                        >> conf_apache
#  echo "      PythonHandler django.core.handlers.modpytho"      >> conf_apache
#  echo "      SetEnv DJANGO_SETTINGS_MODULE settings"           >> conf_apache
#  echo "      PythonDebug On"                                   >> conf_apache
#  echo "      PythonPath "['/var/www/charmming'] + sys.path""   >> conf_apache
#  echo "   </Directory>"                                        >> conf_apache

  echo -e "\nDo you need to setup apache?" 
  read answer
  if [[ "$answer" == "Y" || "$answer" == "y" || "$answer" = "" || "$answer" = "yes" || "$answer" = "Yes" || "$answer" = "YES" ]];
  then
    sudo ln -s /etc/apache2/mods-available/python.load /etc/apache2/mods-enabled/python.load
    cd /tmp/charmming/ 
    echo "/access.log/a\\"                      >  conf_apache
    echo "    \n    <Directory \"/var/www/charmming\">\n      SetHandler python-program\n      PythonHandler django.core.handlers.modpython\n      SetEnv DJANGO_SETTINGS_MODULE settings\n      PythonDebug On\n      PythonPath \"[\'/var/www/charmming\'] + sys.path\"\n    </Directory>" >> conf_apache

    sudo cp /etc/apache2/sites-enabled/000-default /etc/apache2/sites-enabled/000-default.bak
    sed -f conf_apache /etc/apache2/sites-enabled/000-default > 000-default 
    sudo cp 000-default /etc/apache2/sites-enabled/000-default 
    sudo /etc/init.d/apache2 restart
  else
    echo    "-----------------------"
    echo -e "Ok, I will skip this..." 
    echo -e "-----------------------\n"
  fi
}

function config_django {
  echo -e "\nDo you need to setup django and populate the databases?" 
  read answer
  if [[ "$answer" == "Y" || "$answer" == "y" || "$answer" = "" || "$answer" = "yes" || "$answer" = "Yes" || "$answer" = "YES" ]];
  then
    cd $charmming_path
    echo "Answer the following questions this way..." 
    echo "  -> yes" 
    echo "  -> $charmming_admin" 
    echo "  -> email address" 
    echo "  -> $charmming_passwd"
    echo "  -> $charmming_passwd"
    sleep 3 
# sync datbases 
    sudo python manage.py syncdb 

    echo -e " Creating initial groups...\n"
    echo "INSERT INTO auth_group (id,name) VALUES (1,'preapprove');" > /tmp/charmming/groups.sql
    echo "INSERT INTO auth_group (id,name) VALUES (2,'student');" >> /tmp/charmming/groups.sql
    echo "INSERT INTO auth_group (id,name) VALUES (3,'lesson');" >> /tmp/charmming/groups.sql
    mysql -u charmming --password="$database_passwd" -D charmming < /tmp/charmming/groups.sql
    mysql -u charmming --password="$database_passwd" -D charmming < /var/www/charmming/qsar_tables.sql

  else
    echo    "-----------------------"
    echo -e "Ok, I will skip this..." 
    echo -e "-----------------------\n"
  fi
} 

function start_schedd {
  echo -e "\nDo you need to start the charmming scheduler?" 
  read answer
  if [[ "$answer" == "Y" || "$answer" == "y" || "$answer" = "" || "$answer" = "yes" || "$answer" = "Yes" || "$answer" = "YES" ]];
  then
    sudo sed -i "s/passwd\=\"charmming\"/passwd\=\"$charmming_passwd\"/"    $charmming_path/scheduler/schedd.py
    sudo chown -R schedd.www-data /home/schedd/
    sudo chmod -R g+w /home/schedd/
    sudo su - schedd -c "PYTHONPATH=/var/www/charmming $charmming_path/scheduler/schedd.py"
#   echo -e "\n You are almost done, one last thing... Please execute the following commands: \n" 
#   echo -e "  sudo su - schedd" 
#   echo -e "  PYTHONPATH=/var/www/charmming $charmming_path/scheduler/schedd.py\n" 
    echo -e "Congratulations! After this you should have a working version of CHARMMing! :-)\n" 
  else
    echo    "-----------------------"
    echo -e "Ok, I will skip this..." 
    echo -e "-----------------------\n"
  fi
}


# -------------
#  Main Script
# -------------

EXPECTED_ARGS=1
E_BADARGS=-1

ARGV=( $@ )
#echo "ARGV array is: ${ARGV[@]}"
#echo "ARGV array is: ${ARGV[0]}"
#echo $# 

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` -install | -uninstall "
  echo "" 
  echo "Before installing please edit installer.sh and set the "
  echo "following preferences in the set_install_config function :" 
  echo "" 
  echo -e "\t charmming_admin" 
  echo -e "\t charmming_passwd" 
  echo -e "\t database_passwd" 
  echo "" 
  exit $E_BADARGS
fi

if [ "${ARGV[0]}" == "-install" ] 
then 
  echo "" 
  echo -e "    Take Us Up to Warp Speed Mr. Sulu..." 
  echo "" 

  echo "Before proceeding with the installation, please look over the install notes"
  echo "located at http://charmmtutorial.org/index.php/Installation_of_CHARMMing"
  echo ""
  echo "This page provide information on prerequisite software that cannot be"
  echo "distributed woith CHARMMing, as well as known issues with the installation"
  echo "and their work-arounds. So, run don't walk, to:"
  echo ""
  echo "http://charmmtutorial.org/index.php/Installation_of_CHARMMing"
  echo ""
  echo ""

# set config options 
  set_install_config

#install dependencies
  dependencies

#download and install charmming, charmminglib, and charmming-utils 
  download_charmming

#configure mysql 
  mysql_config 

#add schedd user 
  schedd_config

# configure torque 
  torque_config 

# configure apache2
  config_apache

# configure django 
  config_django

# Start scheduler (i.e. schedd.py) 
  start_schedd
  
  echo "Installation is complete ... please review charmming_config.py to make sure that everything is set up correctly."
  echo "At the very least, you will need to provide CHARMM binaries and put their location in charmming_config.py"

elif [ "${ARGV[0]}" == "-uninstall" ]
then 
  echo -e "\nDo you want to uninstall charmming?" 
  read answer
  if [[ "$answer" == "Y" || "$answer" == "y" || "$answer" = "" || "$answer" = "yes" || "$answer" = "Yes" || "$answer" = "YES" ]];
  then
    echo -e "\nOkey dokey Doctor Jones! Hold on to your potato!\n"

# set config options 
    set_install_config

    mkdir -p /tmp/charmming/ 
    sudo /etc/init.d/apache2 stop 
    sudo rm -rf $charmming_path
    sudo rm -rf $charmming_utils
    sudo killall -9 schedd.py 
    sudo /etc/init.d/torque-mom stop
    sudo /etc/init.d/torque-scheduler stop
    sudo /etc/init.d/torque-server stop
    sudo killall -9 pbs_server pbs_sched pbs_mom
    sudo userdel -r schedd 
    sudo cp /etc/apache2/sites-enabled/000-default.bak /etc/apache2/sites-enabled/000-default
    sudo rm -rf /var/spool/torque 
    cd /tmp/charmming/
    echo "drop database charmming;"                     >  cleanup_db.sql 
    echo 'delete from user where User="charmming"; '    >> cleanup_db.sql
    mysql -u root --password="$database_passwd" -D mysql < cleanup_db.sql

    echo -e "\nDo you want to remove all packages installed by charmming?" 
    read answer
      if [[ "$answer" == "Y" || "$answer" == "y" || "$answer" = "" || "$answer" = "yes" || "$answer" = "Yes" || "$answer" = "YES" ]];
      then
        sudo apt-get remove subversion python-django openbabel torque-server torque-scheduler torque-client torque-common torque-mom apache2 mysql-server mysql-client mysql-common libapache2-mod-python python-mysqldb python-numpy
      else
        echo    "-----------------------"
        echo -e "Ok, I will skip this..." 
        echo -e "-----------------------\n"
      fi
    fi 
fi 

exit 


