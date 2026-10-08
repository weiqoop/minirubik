#include <stdio.h>
#include <stdint.h>

#include "stage2-table.h"
#include "stage3-transition-table.h"
#include "stage3-pdb.h"

static void emit_u32(const uint32_t *values, uint32_t count)
{
    for (uint32_t index = 0; index < count; ++index) {
        printf("    .word %u\n", values[index]);
    }
}

int main(void)
{
    puts(".section .rodata");
    puts(".align 2");
    puts(".global stage2_ranks");
    puts("stage2_ranks:");
    emit_u32(stage2_ranks, STAGE2_TABLE_COUNT);

    return 0;
}