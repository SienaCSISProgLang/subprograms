| To force emacs to compile use: -*- compile-command: "make" -*-
| Starter for assembly-only code.  main() called from Go button in simple.c
|
| Simple main program to call gcd function with 2 args, result in d0
|
| Jim Teresco, Williams College, CS237
| Based on code from Duane Bailey
| $Id: main.s 124 2006-10-03 21:27:40Z terescoj $
|
	.include "cs237.h"
	.globl main
	.text
| no locals or params being used
main:
	link	%a6, #0
	| first just trace a simple example in the debugger
	move.w	#24, -(%sp)
	move.w	#36, -(%sp)
	brk
	bsr	gcd
	add.l	#4, %sp
	brk
	| now draw pixels on the screen where gcd(row,col) = 1
	clr.w	%d6  | row in d6
rloop:	cmp.w	#160,%d6
	bge	endrloop
	clr.w	%d7  | row in d7
cloop:	cmp.w	#160,%d7
	bge	endcloop
	move.w	%d7,-(%sp)
	move.w	%d6,-(%sp)
	bsr	gcd
	add.l	#4,%sp
	cmp.w	#1,%d0
	bne	overif
	move.w	%d7,-(%sp)
	move.w	%d6,-(%sp)
	oscall	WinDrawPixel
	add.l	#4,%sp
overif:	
	add.w	#1,%d7
	bra	cloop
endcloop:
	add.w	#1,%d6
	bra	rloop
endrloop:
	unlk	%a6
	rts
