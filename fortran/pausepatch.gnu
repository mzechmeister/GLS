# Version v03  19.04.2012
# by Mathias Zechmeister
# in gnuplot v4.2.3 mouse zooming is not possible with "pause mouse button2" but with "pause mouse key"
# however then bind does not work
# here is a small work around
# disadvantage: double click/press is neccessary to switch between key and zoom mode
# this feels a bit notchy

pause mouse key # "zooming" is active, any keypress switches to the next mode

# another small patch; to prevent at least for the custom bindings the double press
# print MOUSE_CHAR, MOUSE_KEY
BIND_CHAR=(MOUSE_CHAR eq "$") ?"\\$" :\
          (MOUSE_KEY == 1008)?"Left" :\
          (MOUSE_KEY == 1010)?"Right":\
          (MOUSE_KEY == 1009)?"Up"   :\
          (MOUSE_KEY == 1011)?"Down" :\
          (MOUSE_KEY == 1035)?"+" :\
          (MOUSE_KEY == 1037)?"-" :\
           MOUSE_CHAR
BIND_CHAR=(MOUSE_CTRL == 0)?BIND_CHAR :"Ctrl-".BIND_CHAR
BIND_CHAR=(MOUSE_ALT  == 0)?BIND_CHAR :"Alt-".BIND_CHAR
# print BIND_CHAR
set macro
# sed script to search the bind
# combine the muliple lines to one line
# and grep the binding
# requires that bindkey and bindexpression are enclosed by double quotes
#cmd='sed -n ":again /\\\\$/{N;s/\\\\\\n//; t again}; /bind *\"'.BIND_CHAR.'\"/{s/.*\"'.BIND_CHAR.'\" *\"//; s/\"//p}" ~/zoom.gnu'
zoomfile='/home/raid0/zechmeister/programs/GLS/zoom.gnu'
cmd='sed -n ":again /\\\\$/{N;s/\\\\\\n//; t again}; /bind *\"'.BIND_CHAR.'\"/{s/.*\"'.BIND_CHAR.'. *.//; s/\"//p}" '.zoomfile
# print cmd
a=system(cmd);
@a

#a=system('sed -n "/\\\\$/{N;s/\\\\\\n//}; /.*\"'.MOUSE_CHAR.'/{s/.*\"'.MOUSE_CHAR.'\" *\"//; s/\"//p}" ~/zoom.gnu')

if(MOUSE_KEY!=113) pause mouse button1,button2,button3; # now binding active; mouse click will return to key mode
if(MOUSE_KEY!=113) reread # 'q' to quit
# if(MOUSE_KEY!=113) pause mouse button1,button2,button3; reread # not good when hasty clicking

