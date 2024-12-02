# Intro: Get your computer ready
###### tags: `basic commands`, `terminal`
# Software

We will need at least pieces of software.
the bash terminal and a program to transfer files from and to your computer

what is a terminal? check [here](https://itconnect.uw.edu/learn/workshops/online-tutorials/web-publishing/what-is-a-terminal/)

what is a file transfer software ? check [here](https://en.wikipedia.org/wiki/FileZilla)
terminal
## mac

a full terminal, compaltible with unix system, is pre installed. [Link](https://support.apple.com/guide/terminal/open-or-quit-terminal-apd5265185d-f365-44cb-8b09-71a064a42125/mac).
## windows

You need to activate or install it. [how to install an ubuntu terminal in win 10](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) .
## linux

If you are running linux, you are already familiar with this. But just in case, check [here](https://askubuntu.com/questions/38162/what-is-a-terminal-and-how-do-i-open-and-use-it)
# filezilla

download the appropiate version (Client) [here](https://filezilla-project.org/)
# Firefox

Many outputs use .html format. I know firefox don't have any trouble opening them. Maybe if you are using other browser, you can run into troubles. Let's have it alredy installed, you can still use your favorite browser for your stuff.

# text editor.
## windows

You will be using Notepad a lot.
## mac

you have text edit. Just follow the instruction [here](https://www.techjunkie.com/textedit-plain-text-mode/#:~:text=TextEdit%20opens%20a%20new%20document,in%20the%20TextEdit%20menu%20bar.) to make it plain text
## linux

I use the text edit from gedit. It is a small and lightweight text editor for the GNOME Desktop
# Connect to the MBL BPC computing server
...now the fun begings:
Open the terminal and type


`shh -X 'user'@evol5.mbl.edu`

`password: type here`

![](https://i.imgur.com/bKqGoDR.png)

Actually, make sure you have installed "ssh" in your bash terminal. It's a progran to connect to servers. Why the -X? We'll see it in another time.

When you get in after the password, that is the head node. WE NEVER EVER UNDER NO CIRCUNSTANCE work there.

![](https://i.imgur.com/10EjljO.png)


however. there is a series of commands you can use, like browse among dir and files. But better keep going:

`status` #list the servers and clusters. 

![](https://i.imgur.com/FLWkb97.png)


`news` #guess what...

`usage` #list your quota and storage use.

`handbook` # **check this one out***

now, log to the machine you will be using for running the programs

eg:
`ssh -X luna`

![](https://i.imgur.com/OWeMWhi.png)


Sometimes you have to retype your password at this point.

So now, you are in...but WHERE?

type `pwd` and will show the path to your home dir
`/users/your-name`

but we'll be working mostly in a group directory; so I can see your anaysis. 

`/groups/fr2020`
`/groups/arklab`
