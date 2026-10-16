#include <stdint.h>
#include <stdio.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#endif

#include "stage2-table.h"
#include "stage3-transition-table.h"
#include "stage3-pdb.h"

enum { CUBIES = 7, ORIENTATIONS = 729, MOVES = 9, MAX_DEPTH = 6 };

typedef struct { uint8_t p[CUBIES], o[CUBIES]; } state_t;
typedef struct {
    uint8_t best_total, best_depth, best_forward[MAX_DEPTH];
    uint16_t best_permutation, best_orientation;
} search_result_t;

typedef struct {
    uint64_t nodes[MAX_DEPTH + 1];
    uint64_t table_queries;
    uint64_t heuristic_queries;
    uint64_t heuristic_prunes;
} search_stats_t;

typedef struct {
    uint16_t permutation, orientation;
    uint8_t depth, next_move, entered;
    int8_t previous_face;
} search_frame_t;

enum { SEARCH_STACK_CAPACITY = MAX_DEPTH + 1 };

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

static double wall_seconds(void)
{
#ifdef _WIN32
    LARGE_INTEGER counter, frequency;
    QueryPerformanceCounter(&counter);
    QueryPerformanceFrequency(&frequency);
    return (double) counter.QuadPart / (double) frequency.QuadPart;
#else
    return 0.0;
#endif
}

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

static uint8_t heuristic(uint16_t permutation, uint16_t orientation,
                         search_stats_t *stats)
{
    ++stats->heuristic_queries;
    uint8_t permutation_distance =
        stage3_permutation_pdb[permutation];
    uint8_t orientation_distance =
        stage3_orientation_pdb[orientation];
    return permutation_distance > orientation_distance
               ? permutation_distance : orientation_distance;
}

static int stack_push(search_frame_t *stack, uint8_t *stack_size,
                      uint16_t permutation, uint16_t orientation,
                      uint8_t depth, int8_t previous_face)
{
    if (*stack_size == SEARCH_STACK_CAPACITY)
        return 0;
    stack[*stack_size] = (search_frame_t) {
        permutation, orientation, depth, 0, 0, previous_face
    };
    ++*stack_size;
    return 1;
}

static int search_iterative(uint16_t permutation, uint16_t orientation,
                            uint8_t *forward, search_result_t *result,
                            search_stats_t *stats)
{
    search_frame_t stack[SEARCH_STACK_CAPACITY];
    uint8_t stack_size = 0;
    if (!stack_push(stack, &stack_size, permutation, orientation, 0, -1))
        return 0;
    while (stack_size != 0) {
        search_frame_t *frame = &stack[stack_size - 1U];
        if (!frame->entered) {
            ++stats->nodes[frame->depth];
            uint32_t rank = (uint32_t) frame->permutation * ORIENTATIONS +
                            frame->orientation;
            uint32_t index;
            ++stats->table_queries;
            if (table_find(rank, &index)) {
                uint8_t reverse_depth = STAGE2_DISTANCE(stage2_info[index]);
                uint8_t total = (uint8_t) (frame->depth + reverse_depth);
                if (total < result->best_total) {
                    result->best_total = total;
                    result->best_depth = frame->depth;
                    memcpy(result->best_forward, forward, frame->depth);
                    result->best_permutation = frame->permutation;
                    result->best_orientation = frame->orientation;
                }
            }
            if (frame->depth == MAX_DEPTH ||
                frame->depth >= result->best_total) {
                --stack_size;
                continue;
            }
            uint8_t lower_bound = (uint8_t) (frame->depth +
                                              heuristic(frame->permutation,
                                                        frame->orientation,
                                                        stats));
            if (lower_bound >= result->best_total) {
                ++stats->heuristic_prunes;
                --stack_size;
                continue;
            }
            frame->entered = 1;
            frame->next_move = 0;
        }
        if (frame->next_move >= MOVES) {
            --stack_size;
            continue;
        }
        uint8_t move = frame->next_move++;
        uint8_t face = (uint8_t) (move / 3U);
        if ((int8_t) face == frame->previous_face)
            continue;
        uint16_t next_permutation = frame->permutation;
        uint16_t next_orientation = frame->orientation;
        uint8_t turns = (uint8_t) (move % 3U + 1U);
        for (uint8_t turn = 0; turn < turns; ++turn) {
            next_permutation =
                stage3_permutation_transition[face][next_permutation];
            next_orientation =
                stage3_orientation_transition[face][next_orientation];
        }
        forward[frame->depth] = move;
        if (!stack_push(stack, &stack_size, next_permutation, next_orientation,
                        (uint8_t) (frame->depth + 1U), (int8_t) face))
            return 0;
        continue;
    }
    return 1;
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

int main(int argc, char **argv)
{
    static const char input_text[] = "21345671111111";
    const char *input_text_arg = input_text;
    state_t input;
    uint8_t forward[MAX_DEPTH], solution[MAX_DEPTH + 5], length;
    search_result_t result = {UINT8_MAX, 0, {0}, 0, 0};
    search_stats_t stats = {{0}, 0, 0, 0};
    if (argc > 2 || (argc == 2 && !parse_state(argv[1], &input)))
        return 1;
    if (argc == 2)
        input_text_arg = argv[1];
    if (argc != 2 && !parse_state(input_text_arg, &input))
        return 1;
    uint32_t input_rank = rank_state(&input);
    double start = wall_seconds();
    if (!search_iterative((uint16_t) (input_rank / ORIENTATIONS),
                          (uint16_t) (input_rank % ORIENTATIONS), forward,
                          &result, &stats))
        return 1;
    if (result.best_total == UINT8_MAX ||
        !make_solution(&input, &result, solution, &length))
        return 1;
    double finish = wall_seconds();
    printf("input: %s\nsolution:", input_text_arg);
    for (uint8_t i = 0; i < length; ++i)
        printf("%s%s", i == 0 ? "" : " ", move_names[solution[i]]);
    printf("\nlength: %u\n", length);
    puts("verification: solved (rank 0)");
    for (uint8_t depth = 0; depth <= MAX_DEPTH; ++depth)
        printf("nodes_depth_%u: %llu\n", depth,
               (unsigned long long) stats.nodes[depth]);
    printf("nodes_total: %llu\n",
           (unsigned long long) stats.nodes[0] + stats.nodes[1] +
           stats.nodes[2] + stats.nodes[3] + stats.nodes[4] +
           stats.nodes[5] + stats.nodes[6]);
    printf("table_queries: %llu\n", (unsigned long long) stats.table_queries);
    printf("heuristic_queries: %llu\n",
           (unsigned long long) stats.heuristic_queries);
    printf("heuristic_prunes: %llu\n",
           (unsigned long long) stats.heuristic_prunes);
    printf("host_search_seconds: %.6f\n", finish - start);
    printf("pdb_data_bytes: %zu\n",
           sizeof stage3_permutation_pdb + sizeof stage3_orientation_pdb);
    return 0;
}