# CalVIC

MATLAB functions for calibrating the Variable Infiltration Capacity model. Currently implemented: Shuffled Complex Evolution (SCE) algorithm. There are two branches:

* master - for calibrating the distributed IRB model
* ucrb - for calibrating the lumped UCRB VIC model

The core of this code is taken from Q. Duan's 2004 MATLAB SCE-UA implementation, available via the MATLAB file exchange.

Revision to make: rewrite CalVIC to be more flexible with the kinds of inputs it accepts. Make it so the user does not need to modify the source code.