#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif

/* 引入Stage 3 預先計算版本的搜尋程式*/
#define main stage3_precomputed_original_main
#include "stage3-precomputed.c"
#undef main

enum { H3_STATES = 5040 * 729, H3_FACES = 3 };

static double h3_wall_seconds(void)
{
#ifdef _WIN32
    LARGE_INTEGER counter, frequency;
    QueryPerformanceCounter(&counter);
    QueryPerformanceFrequency(&frequency);
    return (double) counter.QuadPart / (double) frequency.QuadPart;
#else
    return (double) clock() / (double) CLOCKS_PER_SEC;
#endif
}

static void h3_unrank_state(uint32_t rank, state_t *state)
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

static uint32_t h3_build_oracle(uint8_t *distance, uint32_t *queue)
{
    uint32_t head = 0, tail = 1;
    memset(distance, UINT8_MAX, H3_STATES);
    distance[0] = 0;
    queue[0] = 0;
    while (head < tail) {
        uint32_t here = queue[head++];
        uint16_t permutation = (uint16_t) (here / ORIENTATIONS);
        uint16_t orientation = (uint16_t) (here % ORIENTATIONS);
        uint8_t next_distance = (uint8_t) (distance[here] + 1U);
        for (uint8_t face = 0; face < H3_FACES; ++face) {
            uint16_t next_permutation = permutation;
            uint16_t next_orientation = orientation;
            for (uint8_t turn = 0; turn < 3; ++turn) {
                next_permutation =
                    stage3_permutation_transition[face][next_permutation];
                next_orientation =
                    stage3_orientation_transition[face][next_orientation];
                uint32_t there = (uint32_t) next_permutation * ORIENTATIONS +
                                  next_orientation;
                if (distance[there] == UINT8_MAX) {
                    distance[there] = next_distance;
                    queue[tail++] = there;
                }
            }
        }
    }
    return tail;
}

static int h3_verify_solution(uint32_t rank, const uint8_t *solution,
                              uint8_t length)
{
    state_t state;
    h3_unrank_state(rank, &state);
    for (uint8_t i = 0; i < length; ++i)
        state = apply_move(state, solution[i]);
    return rank_state(&state) == 0;
}

static int h3_check_state(uint32_t rank, uint8_t expected,
                          uint32_t *successes, uint32_t *failures)
{
    uint8_t forward[MAX_DEPTH], solution[MAX_DEPTH + 5], length;
    search_result_t result = {UINT8_MAX, 0, {0}, 0, 0};
    state_t input;
    uint16_t permutation = (uint16_t) (rank / ORIENTATIONS);
    uint16_t orientation = (uint16_t) (rank % ORIENTATIONS);
    h3_unrank_state(rank, &input);
    search(permutation, orientation, 0, -1, forward, &result);
    if (result.best_total == UINT8_MAX ||
        !make_solution(&input, &result, solution, &length)) {
        fprintf(stderr, "H3 failure: rank=%u expected=%u actual=none reason=no solution\n",
                rank, expected);
        ++*failures;
        return 0;
    }
    if (length != expected) {
        fprintf(stderr, "H3 failure: rank=%u expected=%u actual=%u reason=non-optimal\n",
                rank, expected, length);
        ++*failures;
        return 0;
    }
    if (!h3_verify_solution(rank, solution, length)) {
        fprintf(stderr, "H3 failure: rank=%u expected=%u actual=%u reason=does not solve\n",
                rank, expected, length);
        ++*failures;
        return 0;
    }
    ++*successes;
    return 1;
}

int main(int argc, char **argv)
{
    uint32_t limit = 1000;
    int full = 0;
    uint8_t *oracle;
    uint32_t *queue;
    uint32_t successes = 0, failures = 0, checked = 0;
    uint32_t oracle_count;
    double start, finish;
    if (argc == 2 && strcmp(argv[1], "--all") == 0)
        full = 1;
    else if (argc == 3 && strcmp(argv[1], "--limit") == 0) {
        char *end;
        unsigned long value = strtoul(argv[2], &end, 10);
        if (*argv[2] == '\0' || *end != '\0' || value == 0 || value > H3_STATES) {
            fprintf(stderr, "usage: %s [--limit N|--all]\n", argv[0]);
            return 2;
        }
        limit = (uint32_t) value;
    } else if (argc != 1) {
        fprintf(stderr, "usage: %s [--limit N|--all]\n", argv[0]);
        return 2;
    }
    if (full)
        limit = H3_STATES;

    oracle = malloc(H3_STATES);
    queue = malloc((size_t) H3_STATES * sizeof *queue);
    if (!oracle || !queue) {
        fputs("allocation failed\n", stderr);
        free(oracle);
        free(queue);
        return 1;
    }
    start = h3_wall_seconds();
    oracle_count = h3_build_oracle(oracle, queue);
    if (oracle_count != H3_STATES) {
        fprintf(stderr, "oracle incomplete: visited=%u expected=%u\n",
                oracle_count, H3_STATES);
        free(oracle);
        free(queue);
        return 1;
    }
    for (uint32_t rank = 0; rank < limit; ++rank) {
        ++checked;
        h3_check_state(rank, oracle[rank], &successes, &failures);
    }
    finish = h3_wall_seconds();
    printf("mode: %s\nchecked_states: %u\nsuccesses: %u\nfailures: %u\n",
           full ? "all" : "limit", checked, successes, failures);
    printf("oracle_states: %u\n", oracle_count);
    printf("wall_clock_seconds: %.6f\n", finish - start);
    printf("estimated_full_seconds_from_sample: %.6f\n",
           checked == 0 ? 0.0 : (finish - start) * H3_STATES / checked);
    free(oracle);
    free(queue);
    return failures == 0 ? 0 : 1;
}
