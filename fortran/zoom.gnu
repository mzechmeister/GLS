# Version v03  02.04.2012
# by Mathias Zechmeister
# this script allows to scroll and zoom with keyboard keys
# (wxt does not support KP_Add and KP_Subtract, i.e Num + and Num -, and ctrl+char bindings)
# bindings can be reset with:  bind! or reset bind

DELTA_X(a)=a*(GPVAL_X_MAX-GPVAL_X_MIN)
DELTA_Y(a)=a*(GPVAL_Y_MAX-GPVAL_Y_MIN) 
SCR_X=.25     # x scroll fraction
SCR_Y=.25     # y scroll fraction
ZOOM_X=0.05   # x zoom factor
ZOOM_Y=0.05   # y zoom factor
#bind "XXX" "print 'aaa'"
# scrolling:
bind "Right"       "set xrange [GPVAL_X_MIN+DELTA_X(SCR_X):GPVAL_X_MAX+DELTA_X(SCR_X)]; replot"
bind "Left"        "set xrange [GPVAL_X_MIN-DELTA_X(SCR_X):GPVAL_X_MAX-DELTA_X(SCR_X)]; replot"
bind "Up"          "set yrange [GPVAL_Y_MIN+DELTA_Y(SCR_Y):GPVAL_Y_MAX+DELTA_Y(SCR_Y)]; replot"
bind "Down"        "set yrange [GPVAL_Y_MIN-DELTA_Y(SCR_Y):GPVAL_Y_MAX-DELTA_Y(SCR_Y)]; replot"
# scroll page
bind "Alt-Right"   "set xrange [GPVAL_X_MIN+DELTA_X(1.)   :GPVAL_X_MAX+DELTA_X(1.)];    replot"
bind "Alt-Left"    "set xrange [GPVAL_X_MIN-DELTA_X(1.)   :GPVAL_X_MAX-DELTA_X(1.)];    replot"
bind "Alt-Up"      "set xrange [GPVAL_Y_MIN+DELTA_X(1.)   :GPVAL_Y_MAX+DELTA_X(1.)];    replot"
bind "Alt-Down"    "set xrange [GPVAL_Y_MIN-DELTA_X(1.)   :GPVAL_Y_MAX-DELTA_X(1.)];    replot"

# zooming:
bind "Ctrl-Right"  "set xrange [GPVAL_X_MIN+DELTA_X(ZOOM_X):GPVAL_X_MAX-DELTA_X(ZOOM_X)]; replot" # zoom in x
bind "Ctrl-Left"   "set xrange [GPVAL_X_MIN-DELTA_X(ZOOM_X):GPVAL_X_MAX+DELTA_X(ZOOM_X)]; replot" # zoom out x
bind "Ctrl-Up"     "set yrange [GPVAL_Y_MIN-DELTA_Y(ZOOM_Y):GPVAL_Y_MAX+DELTA_Y(ZOOM_Y)]; replot" # zoom out y
bind "Ctrl-Down"   "set yrange [GPVAL_Y_MIN+DELTA_Y(ZOOM_Y):GPVAL_Y_MAX-DELTA_Y(ZOOM_Y)]; replot" # zoom in y

bind "KP_Add"      "set xrange [GPVAL_X_MIN+DELTA_X(ZOOM_X):GPVAL_X_MAX-DELTA_X(ZOOM_X)];\
                    set yrange [GPVAL_Y_MIN+DELTA_Y(ZOOM_Y):GPVAL_Y_MAX-DELTA_Y(ZOOM_Y)]; replot" # zoom in both
bind "+"           "set xrange [GPVAL_X_MIN+DELTA_X(ZOOM_X):GPVAL_X_MAX-DELTA_X(ZOOM_X)];\
                    set yrange [GPVAL_Y_MIN+DELTA_Y(ZOOM_Y):GPVAL_Y_MAX-DELTA_Y(ZOOM_Y)]; replot" # zoom in both
bind "KP_Subtract" "set xrange [GPVAL_X_MIN-DELTA_X(ZOOM_X):GPVAL_X_MAX+DELTA_X(ZOOM_X)];\
                    set yrange [GPVAL_Y_MIN-DELTA_Y(ZOOM_Y):GPVAL_Y_MAX+DELTA_Y(ZOOM_Y)]; replot" # zoom out both
bind "-"           "set xrange [GPVAL_X_MIN-DELTA_X(ZOOM_X):GPVAL_X_MAX+DELTA_X(ZOOM_X)];\
                    set yrange [GPVAL_Y_MIN-DELTA_Y(ZOOM_Y):GPVAL_Y_MAX+DELTA_Y(ZOOM_Y)]; replot" # zoom out both

bind "C"      "set cbrange [*:*]; replot"                   # reset color range
bind "Alt-C"  "set cbrange [GPVAL_CB_MIN:GPVAL_CB_MAX];"    # save  color range (resetted by 'u', use 'U')
bind "Alt-c"  "set cbrange [GPVAL_CB_MIN:GPVAL_CB_MAX];"    # save  color range (resetted by 'u', use 'U')

# resize all
if (GPVAL_VERSION>4.2)\
bind "Ctrl-u" "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX];\
               set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]; replot";\
bind "Alt-u"      "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX];\
               set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]; replot;"  # reset ranges to max/min data
bind "U"  "set xrange [*:*]; set yrange [*:*]; replot"  # total reset
bind "x"      "set xrange [*:*]; replot"
bind "y"      "set yrange [*:*]; replot"

if (GPVAL_VERSION>4.2)\
bind "Ctrl-x" "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]; replot";\
bind "X"      "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]; replot";
bind "Ctrl-y" "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]; replot";\
bind "Y"      "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]; replot";

 # reset ranges to min/max y data

# pan with mouse (press a, then )
bind "a"      "pause mouse; DELTAX=MOUSE_X;        DELTAY=MOUSE_Y;\
               pause mouse; DELTAX=DELTAX-MOUSE_X; DELTAY=DELTAY-MOUSE_Y;\
               set xrange [GPVAL_X_MIN+DELTAX:GPVAL_X_MAX+DELTAX];\
               set yrange [GPVAL_Y_MIN+DELTAY:GPVAL_Y_MAX+DELTAY]; replot"

# center to mouse
bind "c"      "pause mouse; DELTAX=0.5*(GPVAL_X_MAX-GPVAL_X_MIN);\
                            DELTAY=0.5*(GPVAL_Y_MAX-GPVAL_Y_MIN);\
               set xrange [MOUSE_X-DELTAX:MOUSE_X+DELTAX];\
               set yrange [MOUSE_Y-DELTAY:MOUSE_Y+DELTAY]; replot"

bind "G"      "set grid front lt GR_COL=exists('GR_COL')?(GR_COL+1):1; replot"  # toggle various grid colors
bind "Alt-G"  "set grid front lt GR_COL=exists('GR_COL')?(GR_COL-1):1; replot"  # toggle previous grid colors
bind "Alt-g"  "set grid front lt GR_COL=exists('GR_COL')?(GR_COL-1):1; replot"  # toggle previous grid colors

if (GPVAL_VERSION>4.2)\
bind "P"      'PAL_MOD=exists("PAL_MOD")?(PAL_MOD+1)%5:1;\
               eval "set palette col model ".(PAL_MOD==0?"RGB":PAL_MOD==1?"HSV":PAL_MOD==2?"CMY":PAL_MOD==3?"YIQ":"XYZ"); replot'
               # toggling various palette models

bind "B"      'B_FUNC=exists("B_FUNC")?(B_FUNC+1)%36:1; set palette rgb 3,2,B_FUNC; replot'  # toggle blue
bind "Alt-b"  'B_FUNC=exists("B_FUNC")?(B_FUNC-1)%36:1; set palette rgb 3,2,B_FUNC; replot'  # toggle blue
bind "R"      'R_FUNC=exists("R_FUNC")?(R_FUNC+1)%36:1; set palette rgb R_FUNC,3,2; replot'  # toggle red
bind "Alt-r"  'R_FUNC=exists("R_FUNC")?(R_FUNC-1)%36:1; set palette rgb R_FUNC,3,2; replot'  # toggle red

# bind "7"     is built-in toggle ratio
bind "Alt-7"  "set size square; GP_Y_CEN=0.5*(GPVAL_Y_MAX+GPVAL_Y_MIN);\
               set yrange [GP_Y_CEN-DELTA_X(.5):GP_Y_CEN+DELTA_X(.5)]; replot" # set square; keep x fixed
bind "/"      "set size square; GP_X_CEN=0.5*(GPVAL_X_MAX+GPVAL_X_MIN);\
               set xrange [GP_X_CEN-DELTA_Y(.5):GP_X_CEN+DELTA_Y(.5)]; replot" # set square; keep y fixed
               # assumes "Shift-7" == "/" German keyboard

bind "M"      'GP_CBOX=exists("GP_CBOX")?(GP_CBOX+1)%2:1; if (GP_CBOX) unset colorbox; replot;\
               else set colorbox; replot'   # toggle colorbox

bind "Alt-R"  "reset; replot" # reset with "Shift-Alt-r"

