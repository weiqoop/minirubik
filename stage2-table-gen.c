#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    CUBIES = 7,
    PERMUTATIONS = 5040,
    ORIENTATIONS = 729,
    STATES = PERMUTATIONS * ORIENTATIONS,
    MOVES = 9,
    MAX_DEPTH = 5,
    TABLE_STATES = 12224
};

typedef struct {
    uint8_t p[CUBIES], o[CUBIES];
} state_t;

static const uint8_t source[3][CUBIES] = {
    {1, 4, 2, 0, 3, 5, 6},
    {0, 1, 2, 4, 5, 6, 3},
    {0, 2, 5, 3, 1, 4, 6},
};
static const uint8_t twist[3][CUBIES] = {
    {1, 2, 0, 2, 1, 0, 0},
    {0, 0, 0, 1, 2, 1, 2},
    {0, 0, 0, 0, 0, 0, 0},
};
static const uint8_t inverse_move[MOVES] = {
    2, 1, 0, 5, 4, 3, 8, 7, 6
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

static state_t apply_move(state_t state, uint8_t move)
{
    uint8_t turns = (uint8_t) (move % 3U + 1U);
    for (uint8_t i = 0; i < turns; ++i)
        state = quarter_turn(state, (uint8_t) (move / 3U));
    return state;
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

static void build_transitions(uint16_t permutation[3][PERMUTATIONS],
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

static int build_table(const uint16_t permutation[3][PERMUTATIONS],
                       const uint16_t orientation[3][ORIENTATIONS],
                       uint8_t *distance, uint8_t *toward_solved,
                       uint32_t *queue, uint32_t *count)
{
    uint32_t head = 0, tail = 1;
    distance[0] = 1;
    toward_solved[0] = 0;
    queue[0] = 0;
    while (head < tail) {
        uint32_t here = queue[head++];
        uint8_t here_depth = (uint8_t) (distance[here] - 1U);
        if (here_depth == MAX_DEPTH)
            continue;
        uint16_t p = (uint16_t) (here / ORIENTATIONS);
        uint16_t o = (uint16_t) (here % ORIENTATIONS);
        for (uint8_t face = 0; face < 3; ++face) {
            uint16_t next_p = p, next_o = o;
            for (uint8_t turn = 0; turn < 3; ++turn) {
                next_p = permutation[face][next_p];
                next_o = orientation[face][next_o];
                uint32_t there =
                    (uint32_t) next_p * ORIENTATIONS + next_o;
                if (distance[there] == 0) {
                    uint8_t move = (uint8_t) (face * 3U + turn);
                    distance[there] = (uint8_t) (here_depth + 2U);
                    toward_solved[there] = inverse_move[move];
                    queue[tail++] = there;
                }
            }
        }
    }
    *count = tail;
    return tail == TABLE_STATES;
}

static int verify_table(const uint32_t *ranks, const uint8_t *info,
                        uint32_t count)
{
    state_t state;
    if (count != TABLE_STATES || ranks[0] != 0 || info[0] != 0)
        return 0;
    for (uint32_t i = 1; i < count; ++i) {
        if (ranks[i - 1] >= ranks[i])
            return 0;
    }
    for (uint32_t i = 0; i < count; ++i) {
        uint32_t rank = ranks[i];
        uint8_t depth = (uint8_t) (info[i] >> 4);
        uint8_t move = (uint8_t) (info[i] & 0x0fU);
        if (depth > MAX_DEPTH || move >= MOVES)
            return 0;
        unrank_state(rank, &state);
        for (uint8_t step = depth; step > 0; --step) {
            state = apply_move(state, move);
            rank = rank_state(&state);
            uint32_t j;
            for (j = 0; j < count && ranks[j] != rank; ++j)
                ;
            if (j == count || (uint8_t) (info[j] >> 4) != step - 1U)
                return 0;
            move = (uint8_t) (info[j] & 0x0fU);
        }
        if (rank != 0)
            return 0;
    }
    return 1;
}

static int write_array_file(const char *path, const uint32_t *ranks,
                            const uint8_t *info, uint32_t count)
{
    FILE *file = fopen(path, "w");
    if (!file)
        return 0;
    fprintf(file, "#ifndef STAGE2_TABLE_H\n#define STAGE2_TABLE_H\n\n");
    fprintf(file, "#include <stdint.h>\n\n");
    fprintf(file, "enum { STAGE2_TABLE_COUNT = %u };\n\n", count);
    fprintf(file, "/* info: high nibble is distance; low nibble is move. */\n");
    fprintf(file, "#define STAGE2_DISTANCE(value) ((uint8_t) ((value) >> 4))\n");
    fprintf(file, "#define STAGE2_MOVE(value) ((uint8_t) ((value) & 0x0fU))\n\n");
    fprintf(file, "static const uint32_t stage2_ranks[STAGE2_TABLE_COUNT] = {\n");
    for (uint32_t i = 0; i < count; ++i) {
        if (i % 8U == 0)
            fputs("    ", file);
        fprintf(file, "%u%s", ranks[i], i + 1U == count ? "" : ", ");
        if (i % 8U == 7U || i + 1U == count)
            fputc('\n', file);
    }
    fprintf(file, "};\n\n");
    fprintf(file, "static const uint8_t stage2_info[STAGE2_TABLE_COUNT] = {\n");
    for (uint32_t i = 0; i < count; ++i) {
        if (i % 16U == 0)
            fputs("    ", file);
        fprintf(file, "0x%02x%s", info[i], i + 1U == count ? "" : ", ");
        if (i % 16U == 15U || i + 1U == count)
            fputc('\n', file);
    }
    fprintf(file, "};\n\n#endif\n");
    return fclose(file) == 0;
}

int main(int argc, char **argv)
{
    const char *output = argc == 2 ? argv[1] : "stage2-table.h";
    uint16_t permutation[3][PERMUTATIONS], orientation[3][ORIENTATIONS];
    uint8_t *distance = calloc(STATES, sizeof *distance);
    uint8_t *toward_solved = malloc(STATES);
    uint32_t *queue = malloc((size_t) STATES * sizeof *queue);
    uint32_t *ranks = malloc((size_t) TABLE_STATES * sizeof *ranks);
    uint8_t *info = malloc(TABLE_STATES);
    uint32_t count = 0, output_count = 0;
    int result = EXIT_FAILURE;

    if (!distance || !toward_solved || !queue || !ranks || !info) {
        fputs("allocation failed\n", stderr);
        goto cleanup;
    }
    build_transitions(permutation, orientation);
    if (!build_table(permutation, orientation, distance, toward_solved,
                     queue, &count)) {
        fprintf(stderr, "expected %d states, got %u\n", TABLE_STATES, count);
        goto cleanup;
    }
    for (uint32_t rank = 0; rank < STATES; ++rank) {
        if (distance[rank] != 0) {
            ranks[output_count] = rank;
            info[output_count++] = (uint8_t) (((distance[rank] - 1U) << 4) |
                                               toward_solved[rank]);
        }
    }
    if (!verify_table(ranks, info, output_count)) {
        fputs("table verification failed\n", stderr);
        goto cleanup;
    }
    if (!write_array_file(output, ranks, info, output_count)) {
        fprintf(stderr, "could not write %s\n", output);
        goto cleanup;
    }
    printf("generated %u states in %s\n", output_count, output);
    printf("rank data: %zu bytes; info data: %u bytes\n",
           (size_t) output_count * sizeof *ranks, output_count);
    puts("verification: sorted ranks, distances, and toward-solved moves OK");
    result = EXIT_SUCCESS;

cleanup:
    free(distance);
    free(toward_solved);
    free(queue);
    free(ranks);
    free(info);
    return result;
}
