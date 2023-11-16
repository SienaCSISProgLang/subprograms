# factorial.s -- a SPIM-capable factorial computing program
#
# Jim Teresco, CSCI 2500, RPI, Spring 2009
#              CSIS 220, Siena College, Fall 2010
#
        .data
prompt:	.asciiz "Please enter an integer: "
anslabel:	.asciiz "!="
	

        .align 2
        .text
        .globl main
	
factorial:
        # make space for 2 words on the stack
        addi    $sp, $sp, -8

        # save $ra and $a0 on the stack
        sw     $a0, 4($sp)
        sw     $ra, 0($sp)

        slti   $t0, $a0, 1     # is x < 1?
        beq    $t0, $zero, L1  # if no (t0==0), skip over if part

        # x < 1 here, so just return 1
        addi   $v0, $zero, 1   # return value is 1

        # we could restore $a0 and $ra but we know they haven't
        # changed when we take this quick exit, so let's not
        # but we still need to put $sp back the way we found it
        addi   $sp, $sp, 8
        jr     $ra            # return to caller

        # here, x>=1 so we need to make a recursive call
L1:     addi   $a0, $a0, -1   # get x-1 for the parameter
        jal    factorial      # recursive call, will put answer in $v0

        # We now want to set up for the multiply, but we destroyed $a0
        # above, but have it on the stack, so let's load it
        lw     $a0, 4($sp)
        add    $a1, $v0, $zero # put factorial(x-1) result into $a1
        jal    multiply       # multiply $a0*$a1, put result in $v0

        # $v0 is where we want our answer, so no work there
        # but multiply could have changed $a0 and did change $ra
        lw     $ra, 0($sp)    # restore $ra
        lw     $a0, 4($sp)    # restore $a0
        addi   $sp, $sp, 8    # restore $sp
        jr     $ra            # return to caller

# Multiply with bit shifts and logical ops
#  product = 0;
#  while (y) {
#    if (y&1) product += x;
#    y >>= 1;
#    x <<= 1;
#  }
#
# x is in $a0, y in $a1, product in $v0
multiply:	
	addi	$sp, $sp, -8
	sw	$a0, 0($sp)
	sw	$a1, 4($sp)
	add	$v0, $zero, $zero	# product = 0
Loop:	beq	$a1, $zero, End		# escape loop if y==0
	andi	$t0, $a1, 1		# get y&1 into t0
	beq	$t0, $zero, Overif	# branch over if when y&1==0
	add	$v0, $v0, $a0		# product += x
Overif:	srl	$a1, $a1, 1		# y >>= 1
	sll	$a0, $a0, 1		# x <<= 1
	j	Loop
End:
	lw	$a0, 0($sp)
	lw	$a1, 4($sp)
	addi	$sp, $sp, 8
	jr	$ra

main:                           # main() {
	addi $sp, $sp, -12	# make space for 3 words on stack
	sw $s0, 8($sp)		# push $s0
	sw $s1, 4($sp)		# push $s1
	sw $ra, 0($sp)		# push $ra
	
	la $a0, prompt		# load address of prompt string
	addi $v0, $0, 4		# set up for print_string
	syscall			# to print_string

	addi $v0, $0, 5		# set up for read_int
	syscall			# read int into $v0

	add $a0, $0, $v0	# put input number into $a0
	add $s0, $0, $v0	# also put input number into $s0
	jal factorial		# call factorial subroutine, answer in $v0
	add $s1, $0, $v0	# save answer in $s1
	
	addi $v0, $0, 1		# set up for print_int
	add $a0, $0, $s0	# original number into $a0 for printing
	syscall
	la $a0, anslabel	# load address of anslabel string
	addi $v0, $0, 4		# set up for print_string
	syscall			# to print_string
	addi $v0, $0, 1		# set up for print_int
	add $a0, $0, $s1	# answer into $a0 for printing
	syscall
	addi $v0, $0, 11	# set up for print_char
	add $a0, $0, 13		# new line character
	syscall
	addi $v0, $0, 11	# set up for print_char
	add $a0, $0, 10		# line feed character
	syscall
	
	lw $ra, 0($sp)		# restore $ra from stack
	lw $s1, 4($sp)		# restore $s1 from stack
	lw $s0, 8($sp)		# restore $s0 from stack
	addi $sp, $sp, 12	# restore stack pointer
	jr $ra
