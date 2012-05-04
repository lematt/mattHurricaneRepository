#!/bin/csh

set OLDDISPLAY=$DISPLAY
unsetenv DISPLAY

nohup matlab <runPotentialIntensity.m> output &
setenv DISPLAY $OLDDISPLAY
