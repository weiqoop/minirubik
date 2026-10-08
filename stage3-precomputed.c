#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "stage2-table.h"
#include "stage3-transition-table.h"

enum { CUBIES = 7, ORIENTATIONS = 729, MOVES = 9, MAX_DEPTH = 6 };

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
            next_permutation =
                stage3_permutation_transition[face][next_permutation];
            next_orientation =
                stage3_orientation_transition[face][next_orientation];
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
            rank = (uint32_t)
                       stage3_permutation_transition[face][rank / ORIENTATIONS] *
                   ORIENTATIONS +
                   stage3_orientation_transition[face][rank % ORIENTATIONS];
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
    printf("precomputed_transition_bytes: %zu\n",
           sizeof stage3_permutation_transition +
           sizeof stage3_orientation_transition);
    puts("runtime_transition_build: none");
    return 0;
}
