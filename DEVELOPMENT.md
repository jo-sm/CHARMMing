Development
===========

# Preface 

Development for CHARMMing is split into two areas: backend dev, on Django, and frontend dev, on the javascript and css files, as well as the templates. Both are relatively simple to get started with, but in this document are some guidelines and a guide to get you started.

## Supported environments

CHARMMing is supported officially on Mac and Linux. It may work on Windows but is it not officially supported, so we may be able to help in minimal ways but will not be able to substantially troubleshoot issues.

## Tools and Installation

### Backend

**Virtualenv and virtualenvwrapper**

We use virtualenv to keep deployment simple by isolating the development environment from the rest of the machine. In simpler terms, this means that while we may use an older version of Django in our project, using virtualenv will let you not affect your system's packages. To start, make sure that `virtualenv` is installed by running 

`pip install virtualenv`

You may need to run this as sudo. Once this is installed, you'll need to create a virtual environment. You'll want to keep this separate from the project because the environment is specific to your system and should not be checked into git. The recommended method is to use [`virtualenvwrapper`](http://virtualenvwrapper.readthedocs.org/en/latest/) and following the instructions listed to install and create a virtualenv environment. It simplifies the workflow surrounding using `virtualenv`:

```
$ source /path/to/virtualenvwrapper.sh # this can be added to .bash_profile
$ mkvirtualenv charmming
$ workon charmming # should be automatic but just in case it doesn't show up
(charmming)$ 

... some time later

(charmming)$ git commit -m "new feature"
(charmming)$ deactivate
$ 
```

However, it's also possible to do the steps manually, but they are not recommended and not supported. The manual way to handle this is by creating an external directory, such as `~/.venvs` and then running `virtualenv ~/.venvs/charmming`. You would then activate it by running `source ~/.venvs/charmming/bin/source`. You should then see something prepended to your command line after activation:

```
$ virtualenv ~/.venvs/charmming
$ source ~/.venvs/charmming/bin/source
(charmming)$ ...
(charmming)$ deactivate
$ ...
```

Again this is not officially supported, so if there are issues please use virtualenvwrapper or ask Stack Overflow.

*Package installation*

_Before we go on_: Please ensure you're running a virtual environment by seeing something like the following prompt:

```
(charmming)$ 
```

Without this virtual environment, running `pip install` can mess up your Python system packages. Be careful!

We use a custom `pip` installation script that takes care of development and production deployment differences (i.e. there is no need for the test frameworks or dev database adapter to be installed on production). To install packages for CHARMMing, please run `CHARMMING_ENV=development ./install-deps`. This defaults to `development` but also accepts `test` and `production`. We will keep the requirements stable in the master branch, meaning that you will rarely need to update your dependencies. If we update the dependencies we will make a note of this in the changelog. 

### Frontend

To keep workflow simple for designers and frontend developers, we use the powerful `grunt` and `grunt-watch` to take care of preprocessing our SASS files into CSS and our React JSX files into js, as well as transcompiling our ES6 js into ES5. It also allows us to watch any style changes and update without reloading the page. 

*Node*

To run grunt, you first need `node.js`.

If you're on a Mac, you can install node using [Homebrew](brew.sh). There are other methods but we don't officially support them.

If you're on Linux, you'll need to consult Joyent's documentation on [Installation from a Package Manager](https://github.com/joyent/node/wiki/Installing-Node.js-via-package-manager). Each flavor of Linux has a different preferred method of installation, so we can't cover them individually here.

In either case, this may take some time, so grab a coffee or blueberry muffin while you wait.

Once node has been installed, you can continue with the installation of grunt:

```
$ node -v # should output something like v0.10.36. if not, please make sure it's installed correctly
$ npm install -g grunt-cli # may require sudo
```

Grunt shouldn't take too long to install. Once grunt is installed, you can proceed with installing the rest of the node packages required:

```
$ npm install
```

Once this completes, you're ready for development!

## Creating the first user and other setup

When CHARMMing first runs, it doesn't have an administrator account by default. Additionally, while it can run without asking the questions in the setup, there are some useful packages such as CHARMM input script templates (coming soon!), that are useful to install. To run the setup, run `./manage.py firstrun`. This will ask you for your name, email address, and password to create the account, plus if you'd like to install any additional CHARMMing specific packages. 

## Running locally

To run locally, simply run `start.sh`. It accepts a few environment variables:

`SQLITE_FILE`: The SQLite file that local development uses. Defaults to `local.sqlite` in the top-level directory of the project
`CHARMMING_ENV`: The environment that CHARMMing is currently running in. Defaults to `development` but also accepts `production` and `test`.

If you get any errors, please make sure that you've installed the frontend and backend prerequisites correctly.

## That's it!

Now you should be up and running on development! Check it out at [localhost:3000](http://localhost:3000). Happy hacking!

