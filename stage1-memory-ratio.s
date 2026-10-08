.text
.globl main

main:
    lui  t0, 0x10000
    li   t1, 2097152
    li   t2, 1

loop:
    sb   t2, 0(t0)
    addi t0, t0, 1
    addi t1, t1, -1
    bne  t1, zero, loop

hold:
    j hold
