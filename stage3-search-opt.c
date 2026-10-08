#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "stage2-table.h"

enum { CUBIES = 7, PERMUTATIONS = 5040, ORIENTATIONS = 729,
       MOVES = 9, MAX_DEPTH = 6 };

typedef struct { uint8_t p[CUBIES], o[CUBIES]; } state_t;
typedef struct {
    uint8_t best_total, best_depth, best_forward[MAX_DEPTH];
    uint16_t best_permutation, best_orientation;
} search_result_t;

static const char *const move_names[MOVES] = {
    "R", "R2", "R'", "B", "B2", "B'", "D", "D2", "D'"
};
static const uint8_t source[3][CUBIES] = {
    {1, 4, 2, 0, 3, 5, 6}, {0, 1, 2, 4, 5, 6, 3},
    {0, 2, 5, 3, 1, 4, 6}
};
static const uint8_t twist[3][CUBIES] = {
    {1, 2, 0, 2, 1, 0, 0}, {0, 0, 0, 1, 2, 1, 2},
    {0, 0, 0, 0, 0, 0, 0}
};
static uint16_t permutation_transition[3][PERMUTATIONS];
static uint16_t orientation_transition[3][ORIENTATIONS];

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

static void build_transitions(void)
{
    state_t state;
    for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
        unrank_state((uint32_t) rank * ORIENTATIONS, &state);
        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);
            permutation_transition[face][rank] =
                (uint16_t) (rank_state(&next) / ORIENTATIONS);
        }
    }
    for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
        unrank_state(rank, &state);
        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);
            orientation_transition[face][rank] =
                (uint16_t) (rank_state(&next) % ORIENTATIONS);
        }
    }
}

static int parse_state(const char *input, state_t *state)
{
    for (int i = 0; i < 14; ++i) {
        int limit = i < 7 ? 7 : 3;
        if (input[i] < '1' || input[i] > '0' + limit)
            return 0;
        (i < 7 ? state->p : state->o)[i % 7] =
            (uint8_t) (input[i] - '1');
    }
    if (input[14] != '\0')
        return 0;
    uint8_t sum = 0;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        if (state->p[i] >= CUBIES || state->o[i] >= 3)
            return 0;
        for (uint8_t j = 0; j < i; ++j)
            if (state->p[j] == state->p[i])
                return 0;
        sum = (uint8_t) (sum + state->o[i]);
    }
    return sum % 3U == 0;
}

/* Lower-bound search with an early equality return. */
static int table_find(uint32_t rank, uint32_t *index)
{
    uint32_t lo = 0;
    uint32_t hi = STAGE2_TABLE_COUNT - 1U;
    if (stage2_ranks[lo] > rank || stage2_ranks[hi] < rank)
        return 0;
    while (lo <= hi) {
        uint32_t middle = lo + (hi - lo) / 2U;
        uint32_t value = stage2_ranks[middle];
        if (value == rank) {
            *index = middle;
            return 1;
        }
        if (value < rank)
            lo = middle + 1U;
        else {
            if (middle == 0)
                break;
            hi = middle - 1U;
        }
    }
    return 0;
}

static void search(uint16_t permutation, uint16_t orientation,
                   uint8_t depth, int8_t previous_face, uint8_t *forward,
                   search_result_t *result)
{
    uint32_t rank = (uint32_t) permutation * ORIENTATIONS + orientation;
    uint32_t index;
    if (table_find(rank, &index)) {
        uint8_t reverse_depth = STAGE2_DISTANCE(stage2_info[index]);
        uint8_t total = (uint8_t) (depth + reverse_depth);
        if (total < result->best_total) {
            result->best_total = total;
            result->best_depth = depth;
            memcpy(result->best_forward, forward, depth);
            result->best_permutation = permutation;
            result->best_orientation = orientation;
        }
    }
    if (depth == MAX_DEPTH || depth >= result->best_total)
        return;
    for (uint8_t move = 0; move < MOVES; ++move) {
        uint8_t face = (uint8_t) (move / 3U);
        if ((int8_t) face == previous_face)
            continue;
        uint16_t next_permutation = permutation;
        uint16_t next_orientation = orientation;
        uint8_t turns = (uint8_t) (move % 3U + 1U);
        for (uint8_t turn = 0; turn < turns; ++turn) {
            next_permutation = permutation_transition[face][next_permutation];
            next_orientation = orientation_transition[face][next_orientation];
        }
        forward[depth] = move;
        search(next_permutation, next_orientation, (uint8_t) (depth + 1U),
               (int8_t) face, forward, result);
    }
}

static int make_solution(const state_t *input, const search_result_t *result,
                         uint8_t *solution, uint8_t *length)
{
    uint32_t rank = (uint32_t) result->best_permutation * ORIENTATIONS +
                    result->best_orientation;
    uint8_t count = result->best_depth;
    memcpy(solution, result->best_forward, count);
    while (rank != 0) {
        uint32_t index;
        if (!table_find(rank, &index))
            return 0;
        solution[count++] = STAGE2_MOVE(stage2_info[index]);
        uint8_t move = solution[count - 1U];
        uint8_t face = (uint8_t) (move / 3U);
        uint8_t turns = (uint8_t) (move % 3U + 1U);
        for (uint8_t turn = 0; turn < turns; ++turn)
            rank = (uint32_t) permutation_transition[face][rank / ORIENTATIONS] *
                   ORIENTATIONS + orientation_transition[face][rank % ORIENTATIONS];
    }
    state_t state = *input;
    for (uint8_t i = 0; i < count; ++i)
        state = apply_move(state, solution[i]);
    if (rank_state(&state) != 0)
        return 0;
    *length = count;
    return 1;
}

int main(void)
{
    static const char input_text[] = "21345671111111";
    state_t input;
    uint8_t forward[MAX_DEPTH], solution[MAX_DEPTH + 5], length;
    search_result_t result = {UINT8_MAX, 0, {0}, 0, 0};
    build_transitions();
    if (!parse_state(input_text, &input))
        return 1;
    uint32_t input_rank = rank_state(&input);
    search((uint16_t) (input_rank / ORIENTATIONS),
           (uint16_t) (input_rank % ORIENTATIONS), 0, -1, forward, &result);
    if (result.best_total == UINT8_MAX ||
        !make_solution(&input, &result, solution, &length))
        return 1;
    printf("input: %s\nsolution:", input_text);
    for (uint8_t i = 0; i < length; ++i)
        printf("%s%s", i == 0 ? "" : " ", move_names[solution[i]]);
    printf("\nlength: %u\n", length);
    puts("verification: solved (rank 0)");
    return 0;
}
