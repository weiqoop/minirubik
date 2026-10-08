#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stage2-table.h"
#include "stage3-transition-table.h"

enum {
    CUBIES = 7,
    PERMUTATIONS = 5040,
    ORIENTATIONS = 729,
    STATES = PERMUTATIONS * ORIENTATIONS,
    FACES = 3,
    MOVES = 9,
    STAGE2_DEPTH = 5
};

typedef struct {
    uint8_t p[CUBIES], o[CUBIES];
} state_t;

static const uint8_t source[FACES][CUBIES] = {
    {1, 4, 2, 0, 3, 5, 6},
    {0, 1, 2, 4, 5, 6, 3},
    {0, 2, 5, 3, 1, 4, 6}
};
static const uint8_t twist[FACES][CUBIES] = {
    {1, 2, 0, 2, 1, 0, 0},
    {0, 0, 0, 1, 2, 1, 2},
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
    uint32_t p = rank / ORIENTATIONS, o = rank % ORIENTATIONS, factor = 720;
    uint8_t sum = 0;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t choice = (uint8_t) (p / factor);
        p %= factor;
        state->p[i] = available[choice];
        uint8_t limit = (uint8_t) (CUBIES - i);
        for (uint8_t j = choice; j + 1U < limit; ++j)
            available[j] = available[j + 1U];
        if (i < 5)
            factor /= 6U - i;
    }
    for (uint8_t i = 6; i-- > 0;) {
        state->o[i] = (uint8_t) (o % 3U);
        sum = (uint8_t) (sum + state->o[i]);
        o /= 3U;
    }
    state->o[6] = (uint8_t) ((3U - sum % 3U) % 3U);
}

static void build_direct_transitions(uint16_t permutation[FACES][PERMUTATIONS],
                                     uint16_t orientation[FACES][ORIENTATIONS])
{
    state_t state;
    for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
        unrank_state((uint32_t) rank * ORIENTATIONS, &state);
        for (uint8_t face = 0; face < FACES; ++face) {
            state_t next = quarter_turn(state, face);
            permutation[face][rank] =
                (uint16_t) (rank_state(&next) / ORIENTATIONS);
        }
    }
    for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
        unrank_state(rank, &state);
        for (uint8_t face = 0; face < FACES; ++face) {
            state_t next = quarter_turn(state, face);
            orientation[face][rank] =
                (uint16_t) (rank_state(&next) % ORIENTATIONS);
        }
    }
}

static void build_oracle(uint16_t permutation[FACES][PERMUTATIONS],
                         uint16_t orientation[FACES][ORIENTATIONS],
                         uint8_t *distance, uint32_t *queue)
{
    uint32_t head = 0, tail = 1;
    distance[0] = 0;
    queue[0] = 0;
    while (head < tail) {
        uint32_t here = queue[head++];
        uint16_t permutation_rank = (uint16_t) (here / ORIENTATIONS);
        uint16_t orientation_rank = (uint16_t) (here % ORIENTATIONS);
        uint8_t next_distance = (uint8_t) (distance[here] + 1U);
        for (uint8_t face = 0; face < FACES; ++face) {
            uint16_t next_permutation = permutation_rank;
            uint16_t next_orientation = orientation_rank;
            for (uint8_t turn = 0; turn < 3; ++turn) {
                next_permutation = permutation[face][next_permutation];
                next_orientation = orientation[face][next_orientation];
                uint32_t there = (uint32_t) next_permutation * ORIENTATIONS +
                                  next_orientation;
                if (distance[there] == UINT8_MAX) {
                    distance[there] = next_distance;
                    queue[tail++] = there;
                }
            }
        }
    }
    if (tail != STATES)
        fprintf(stderr, "oracle visited %u of %u states\n", tail, STATES);
}

static int table_find(uint32_t rank, uint32_t *index)
{
    uint32_t low = 0, high = STAGE2_TABLE_COUNT;
    while (low < high) {
        uint32_t middle = low + (high - low) / 2U;
        if (stage2_ranks[middle] < rank)
            low = middle + 1U;
        else
            high = middle;
    }
    if (low == STAGE2_TABLE_COUNT || stage2_ranks[low] != rank)
        return 0;
    *index = low;
    return 1;
}

static void report(const char *name, uint32_t checked, uint32_t failures)
{
    printf("%s: %s (checked=%u, failures=%u)\n", name,
           failures == 0 ? "PASS" : "FAIL", checked, failures);
}

int main(void)
{
    static uint16_t direct_permutation[FACES][PERMUTATIONS];
    static uint16_t direct_orientation[FACES][ORIENTATIONS];
    uint8_t *distance = malloc(STATES);
    uint32_t *queue = malloc((size_t) STATES * sizeof *queue);
    uint32_t failures = 0;
    uint32_t total_failures = 0;
    uint32_t checked;
    uint8_t max_oracle_distance = 0;
    uint32_t oracle_count = 0;
    if (!distance || !queue) {
        fputs("allocation failed\n", stderr);
        free(distance);
        free(queue);
        return EXIT_FAILURE;
    }

    memset(distance, UINT8_MAX, STATES);
    build_direct_transitions(direct_permutation, direct_orientation);
    build_oracle(direct_permutation, direct_orientation, distance, queue);
    for (uint32_t rank = 0; rank < STATES; ++rank) {
        if (distance[rank] == UINT8_MAX)
            ++failures;
        else {
            if (distance[rank] > max_oracle_distance)
                max_oracle_distance = distance[rank];
            if (distance[rank] <= STAGE2_DEPTH)
                ++oracle_count;
        }
    }
    checked = STATES;
    report("oracle complete", checked, failures);
    total_failures += failures;
    printf("oracle_max_distance: %u\n", max_oracle_distance);
    printf("oracle_depth_0_to_5_count: %u\n", oracle_count);

    failures = 0;
    for (uint32_t i = 1; i < STAGE2_TABLE_COUNT; ++i)
        if (stage2_ranks[i - 1U] >= stage2_ranks[i])
            ++failures;
    report("stage2 ranks strictly increasing", STAGE2_TABLE_COUNT, failures);
    total_failures += failures;

    failures = 0;
    uint32_t covered = 0;
    for (uint32_t rank = 0; rank < STATES; ++rank) {
        if (distance[rank] <= STAGE2_DEPTH) {
            ++covered;
            uint32_t index;
            if (!table_find(rank, &index))
                ++failures;
        }
    }
    report("stage2 covers oracle depths 0..5", covered, failures);
    total_failures += failures;

    failures = 0;
    uint8_t max_table_distance = 0;
    for (uint32_t i = 0; i < STAGE2_TABLE_COUNT; ++i) {
        if (stage2_ranks[i] >= STATES) {
            ++failures;
            continue;
        }
        uint8_t table_distance = STAGE2_DISTANCE(stage2_info[i]);
        if (table_distance > max_table_distance)
            max_table_distance = table_distance;
        if (distance[stage2_ranks[i]] != table_distance)
            ++failures;
    }
    report("stage2 distances equal oracle", STAGE2_TABLE_COUNT, failures);
    total_failures += failures;
    printf("stage2_max_distance: %u\n", max_table_distance);
    printf("stage2_rank0_distance: %u\n", STAGE2_DISTANCE(stage2_info[0]));
    printf("stage2_rank0_move: %u\n", STAGE2_MOVE(stage2_info[0]));

    failures = 0;
    checked = 0;
    for (uint32_t i = 0; i < STAGE2_TABLE_COUNT; ++i) {
        uint32_t rank = stage2_ranks[i];
        if (rank == 0)
            continue;
        ++checked;
        uint8_t move = STAGE2_MOVE(stage2_info[i]);
        if (rank >= STATES || move >= MOVES)
            ++failures;
        else {
            state_t state;
            unrank_state(rank, &state);
            state_t next_state = apply_move(state, move);
            uint32_t next = rank_state(&next_state);
            if (distance[next] + 1U != distance[rank])
                ++failures;
        }
    }
    report("toward-solved moves reduce distance", checked, failures);
    total_failures += failures;

    failures = 0;
    checked = 0;
    uint16_t max_permutation = 0, max_orientation = 0;
    for (uint8_t face = 0; face < FACES; ++face) {
        for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
            ++checked;
            uint16_t value = direct_permutation[face][rank];
            if (value >= PERMUTATIONS || value !=
                stage3_permutation_transition[face][rank])
                ++failures;
            if (value > max_permutation)
                max_permutation = value;
        }
    }
    report("permutation transitions", checked, failures);
    total_failures += failures;
    printf("permutation_max_value: %u\n", max_permutation);
    printf("permutation_solved_entries: %u, %u, %u\n",
           direct_permutation[0][0], direct_permutation[1][0],
           direct_permutation[2][0]);

    failures = 0;
    checked = 0;
    for (uint8_t face = 0; face < FACES; ++face) {
        for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
            ++checked;
            uint16_t value = direct_orientation[face][rank];
            if (value >= ORIENTATIONS || value !=
                stage3_orientation_transition[face][rank])
                ++failures;
            if (value > max_orientation)
                max_orientation = value;
        }
    }
    report("orientation transitions", checked, failures);
    total_failures += failures;
    printf("orientation_max_value: %u\n", max_orientation);
    printf("orientation_solved_entries: %u, %u, %u\n",
           direct_orientation[0][0], direct_orientation[1][0],
           direct_orientation[2][0]);

    free(distance);
    free(queue);
    return total_failures == 0 && max_oracle_distance == 11 &&
           oracle_count == STAGE2_TABLE_COUNT ? EXIT_SUCCESS : EXIT_FAILURE;
}
