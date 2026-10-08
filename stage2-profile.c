#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stage2-table.h"

enum {
    CUBIES = 7,
    ORIENTATIONS = 729,
    MOVES = 9,
    MAX_DEPTH = 6,
    MAX_CANDIDATES = 64
};

typedef struct {
    uint8_t p[CUBIES], o[CUBIES];
} state_t;

typedef enum {
    PHASE_SEARCH,
    PHASE_REBUILD,
    PHASE_VERIFY
} phase_t;

typedef struct {
    uint64_t dfs_nodes[MAX_DEPTH + 1];
    uint64_t table_calls[3];
    uint64_t table_comparisons[3];
    uint64_t rank_calls[3];
    uint64_t apply_moves[3];
    uint64_t quarter_turns[3];
    uint8_t candidate_count;
    struct {
        uint8_t forward_depth;
        uint8_t reverse_depth;
        uint8_t total;
    } candidates[MAX_CANDIDATES];
    phase_t phase;
} profile_t;

static profile_t profile;

static const char *const move_names[MOVES] = {
    "R", "R2", "R'", "B", "B2", "B'", "D", "D2", "D'"
};
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

static state_t quarter_turn(state_t state, uint8_t face)
{
    profile.quarter_turns[profile.phase]++;
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
    profile.apply_moves[profile.phase]++;
    uint8_t turns = (uint8_t) (move % 3U + 1U);
    for (uint8_t i = 0; i < turns; ++i)
        state = quarter_turn(state, (uint8_t) (move / 3U));
    return state;
}

static uint32_t rank_state(const state_t *state)
{
    profile.rank_calls[profile.phase]++;
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
    profile.table_calls[profile.phase]++;
    while (low < high) {
        uint32_t middle = low + (high - low) / 2U;
        profile.table_comparisons[profile.phase]++;
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

typedef struct {
    uint8_t best_total;
    uint8_t best_depth;
    uint8_t best_forward[MAX_DEPTH];
    state_t best_meeting;
} search_result_t;

static void record_candidate(uint8_t forward_depth, uint8_t reverse_depth,
                             uint8_t total)
{
    if (profile.candidate_count < MAX_CANDIDATES) {
        uint8_t i = profile.candidate_count++;
        profile.candidates[i].forward_depth = forward_depth;
        profile.candidates[i].reverse_depth = reverse_depth;
        profile.candidates[i].total = total;
    }
}

static void search(state_t state, uint8_t depth, int8_t previous_face,
                   uint8_t *forward, search_result_t *result)
{
    profile.dfs_nodes[depth]++;
    uint32_t index;
    if (table_find(rank_state(&state), &index)) {
        uint8_t reverse_depth = STAGE2_DISTANCE(stage2_info[index]);
        uint8_t total = (uint8_t) (depth + reverse_depth);
        if (total < result->best_total) {
            record_candidate(depth, reverse_depth, total);
            result->best_total = total;
            result->best_depth = depth;
            memcpy(result->best_forward, forward, depth);
            result->best_meeting = state;
        }
    }
    if (depth == MAX_DEPTH || depth >= result->best_total)
        return;
    for (uint8_t move = 0; move < MOVES; ++move) {
        uint8_t face = (uint8_t) (move / 3U);
        if ((int8_t) face == previous_face)
            continue;
        forward[depth] = move;
        search(apply_move(state, move), (uint8_t) (depth + 1U),
               (int8_t) face, forward, result);
    }
}

static int make_solution(const state_t *input, const search_result_t *result,
                         uint8_t *solution, uint8_t *length)
{
    state_t state = result->best_meeting;
    uint8_t count = result->best_depth;
    memcpy(solution, result->best_forward, count);
    while (rank_state(&state) != 0) {
        uint32_t index;
        if (!table_find(rank_state(&state), &index))
            return 0;
        solution[count] = STAGE2_MOVE(stage2_info[index]);
        state = apply_move(state, solution[count]);
        ++count;
    }
    state = *input;
    for (uint8_t i = 0; i < count; ++i)
        state = apply_move(state, solution[i]);
    if (rank_state(&state) != 0)
        return 0;
    *length = count;
    return 1;
}

static void print_phase(const char *name, phase_t phase)
{
    printf("%s: table_calls=%llu table_comparisons=%llu rank_calls=%llu "
           "apply_moves=%llu quarter_turns=%llu\n", name,
           (unsigned long long) profile.table_calls[phase],
           (unsigned long long) profile.table_comparisons[phase],
           (unsigned long long) profile.rank_calls[phase],
           (unsigned long long) profile.apply_moves[phase],
           (unsigned long long) profile.quarter_turns[phase]);
}

int main(void)
{
    static const char input_text[] = "21345671111111";
    state_t input;
    uint8_t forward[MAX_DEPTH], solution[MAX_DEPTH + 5], length;
    search_result_t result = {UINT8_MAX, 0, {0}, {{0}, {0}}};
    profile.phase = PHASE_SEARCH;
    if (!parse_state(input_text, &input)) {
        fputs("input parse failed\n", stderr);
        return 1;
    }
    search(input, 0, -1, forward, &result);
    if (result.best_total == UINT8_MAX) {
        fputs("no solution found\n", stderr);
        return 1;
    }
    for (uint8_t i = 0; i < profile.candidate_count; ++i)
        printf("candidate %u: forward_depth=%u reverse_depth=%u total=%u\n",
               (unsigned) i + 1U, profile.candidates[i].forward_depth,
               profile.candidates[i].reverse_depth,
               profile.candidates[i].total);
    profile.phase = PHASE_REBUILD;
    if (!make_solution(&input, &result, solution, &length)) {
        fputs("solution verification failed\n", stderr);
        return 1;
    }
    profile.phase = PHASE_VERIFY;
    state_t check = input;
    for (uint8_t i = 0; i < length; ++i)
        check = apply_move(check, solution[i]);
    if (rank_state(&check) != 0) {
        fputs("final rank is not zero\n", stderr);
        return 1;
    }
    printf("solution:");
    for (uint8_t i = 0; i < length; ++i)
        printf("%s%s", i == 0 ? "" : " ", move_names[solution[i]]);
    printf("\nlength: %u\n", length);
    puts("verification: solved (rank 0)");
    for (uint8_t depth = 0; depth <= MAX_DEPTH; ++depth)
        printf("dfs_depth_%u_nodes: %llu\n", depth,
               (unsigned long long) profile.dfs_nodes[depth]);
    print_phase("search", PHASE_SEARCH);
    print_phase("rebuild", PHASE_REBUILD);
    print_phase("verify", PHASE_VERIFY);
    return 0;
}
