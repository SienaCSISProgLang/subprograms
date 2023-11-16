|
| gcd function expect 2 word args, a and b, result in d0
|
| implements this C function:
| int gcd(int a, int b)
| // pre: a <= b
| {
|    if (a == 0) return b;
|    if (a > b) return gcd(b,a);
|    return gcd(b%a,a);
| }
|
| Jim Teresco, Williams College, CS237
| Based on code from Duane Bailey
| $Id: gcd.s 124 2006-10-03 21:27:40Z terescoj $
|
	.include "cs237.h"
	.globl gcd
	.text
| Assume 0 <= a,b < 255
| Call frame for gcd: symbols assigned offsets from a6
| oldfp   =  0  -- always 0
| retaddr =  4  -- always 4
a       =  8   | always 8
b       = 10   | 8 + sizeof(a)
| end of frame

gcd:    
	link	%a6,#0
	move.w	a(%a6),%d0	| get param a into d0
	tst.w	%d0		| if (a == 0) return b
	bne	overif1
	move.w	b(%a6),%d0	| put param b into d0 for return
	bra	gcdend
overif1:			| if (a > b) return gcd(b,a)
	cmp.w	b(%a6),%d0	| d0 still contains a
	ble	overif2
	move.w	a(%a6),-(%sp)	| push params for
	move.w	b(%a6),-(%sp)	| recursive gcd call
	bsr	gcd		| recursion!
	add.l	#4,%sp
	bra	gcdend		| result in d0, ready for return
overif2:			| return gcd(b%a,a)
	move.w	a(%a6),%d0
	move.w	b(%a6),%d1
	ext.l	%d0		| not strictly nec.; good
	ext.l	%d1		| habit, though
	divu	%d0,%d1
	swap	%d1		| remainder into lowest word
	move.w	%d0,-(%sp)	| push a
	move.w	%d1,-(%sp)	| push result of b%a
	bsr	gcd
	add.l	#4,%sp
gcdend:
	unlk	%a6
	rts    
