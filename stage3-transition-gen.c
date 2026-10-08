#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum { CUBIES = 7, PERMUTATIONS = 5040, ORIENTATIONS = 729 };

typedef struct { uint8_t p[CUBIES], o[CUBIES]; } state_t;

static const uint8_t source[3][CUBIES] = {
    {1, 4, 2, 0, 3, 5, 6}, {0, 1, 2, 4, 5, 6, 3},
    {0, 2, 5, 3, 1, 4, 6}
};
static const uint8_t twist[3][CUBIES] = {
    {1, 2, 0, 2, 1, 0, 0}, {0, 0, 0, 1, 2, 1, 2},
    {0, 0, 0, 0, 0, 0, 0}
};

static state_t quarter_turn(state_t state, uint8_t face)
{
    state_t result;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t from = source[face][i];
        result.p[i] = state.p[from];
        result.o[i] = (uint8_t) ((state.o[from] + twist[face][i]) % 3U);
    }
    return result;
}

static uint32_t rank_state(const state_t *state)
{
    uint32_t p = 0, o = 0;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t smaller = 0;
        for (uint8_t j = (uint8_t) (i + 1U); j < CUBIES; ++j)
            if (state->p[j] < state->p[i])
                ++smaller;
        p = p * (CUBIES - i) + smaller;
    }
    for (uint8_t i = 0; i < 6; ++i)
        o = o * 3U + state->o[i];
    return p * ORIENTATIONS + o;
}

static void unrank_state(uint32_t rank, state_t *state)
{
    uint8_t available[CUBIES] = {0, 1, 2, 3, 4, 5, 6};
    uint32_t p = rank / ORIENTATIONS, o = rank % ORIENTATIONS, f = 720;
    uint8_t sum = 0;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t q = (uint8_t) (p / f);
        p %= f;
        state->p[i] = available[q];
        for (uint8_t j = q; j + 1U < CUBIES - i; ++j)
            available[j] = available[j + 1U];
        if (i < 5)
            f /= 6U - i;
    }
    for (uint8_t i = 6; i-- > 0;) {
        state->o[i] = (uint8_t) (o % 3U);
        sum = (uint8_t) (sum + state->o[i]);
        o /= 3U;
    }
    state->o[6] = (uint8_t) ((3U - sum % 3U) % 3U);
}

static void build(uint16_t permutation[3][PERMUTATIONS],
                  uint16_t orientation[3][ORIENTATIONS])
{
    state_t state;
    for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
        unrank_state((uint32_t) rank * ORIENTATIONS, &state);
        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);
            permutation[face][rank] =
                (uint16_t) (rank_state(&next) / ORIENTATIONS);
        }
    }
    for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
        unrank_state(rank, &state);
        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);
            orientation[face][rank] =
                (uint16_t) (rank_state(&next) % ORIENTATIONS);
        }
    }
}

static int write_table(const char *path,
                       const uint16_t permutation[3][PERMUTATIONS],
                       const uint16_t orientation[3][ORIENTATIONS])
{
    FILE *file = fopen(path, "w");
    if (!file)
        return 0;
    fputs("#ifndef STAGE3_TRANSITION_TABLE_H\n"
          "#define STAGE3_TRANSITION_TABLE_H\n\n"
          "#include <stdint.h>\n\n"
          "enum { STAGE3_PERMUTATIONS = 5040, "
          "STAGE3_ORIENTATIONS = 729 };\n\n", file);
    fputs("static const uint16_t stage3_permutation_transition[3]"
          "[STAGE3_PERMUTATIONS] = {\n", file);
    for (uint8_t face = 0; face < 3; ++face) {
        fputs("    {\n", file);
        for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
            if (rank % 12U == 0)
                fputs("        ", file);
            fprintf(file, "%u%s", permutation[face][rank],
                    rank + 1U == PERMUTATIONS ? "" : ", ");
            if (rank % 12U == 11U || rank + 1U == PERMUTATIONS)
                fputc('\n', file);
        }
        fprintf(file, "    }%s\n", face == 2 ? "" : ",");
    }
    fputs("};\n\nstatic const uint16_t stage3_orientation_transition[3]"
          "[STAGE3_ORIENTATIONS] = {\n", file);
    for (uint8_t face = 0; face < 3; ++face) {
        fputs("    {\n", file);
        for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
            if (rank % 12U == 0)
                fputs("        ", file);
            fprintf(file, "%u%s", orientation[face][rank],
                    rank + 1U == ORIENTATIONS ? "" : ", ");
            if (rank % 12U == 11U || rank + 1U == ORIENTATIONS)
                fputc('\n', file);
        }
        fprintf(file, "    }%s\n", face == 2 ? "" : ",");
    }
    fputs("};\n\n#endif\n", file);
    return fclose(file) == 0;
}

int main(int argc, char **argv)
{
    const char *path = argc == 2 ? argv[1] : "stage3-transition-table.h";
    static uint16_t permutation[3][PERMUTATIONS];
    static uint16_t orientation[3][ORIENTATIONS];
    build(permutation, orientation);
    if (!write_table(path, permutation, orientation)) {
        fprintf(stderr, "could not write %s\n", path);
        return EXIT_FAILURE;
    }
    printf("generated %zu bytes of transition data in %s\n",
           sizeof permutation + sizeof orientation, path);
    return EXIT_SUCCESS;
}
