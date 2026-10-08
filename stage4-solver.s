	.section .rodata
input_state:
	.asciz "21345671111111"

	.section .text
	.globl main

main:
	# Stage 4 initial RV32I skeleton.
	# The reference input is embedded temporarily for reproducible testing.

	la      a0, input_state

	# Exit with code 0.
	li      a0, 0
	li      a7, 93
	ecall
