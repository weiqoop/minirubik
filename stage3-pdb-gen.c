#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif

#include "stage3-transition-table.h"

enum {
    FACES = 3,
    MOVES = 9,
    PERMUTATIONS = 5040,
    ORIENTATIONS = 729,
    STATES = PERMUTATIONS * ORIENTATIONS
};

static uint8_t *bfs_pdb(uint16_t table[3][PERMUTATIONS],
                        uint8_t *max_distance)
{
    uint8_t *distance = malloc(PERMUTATIONS);
    uint16_t *queue = malloc((size_t) PERMUTATIONS * sizeof *queue);
    uint32_t head = 0, tail = 1;
    if (!distance || !queue) {
        free(distance);
        free(queue);
        return NULL;
    }
    memset(distance, UINT8_MAX, PERMUTATIONS);
    distance[0] = 0;
    queue[0] = 0;
    *max_distance = 0;
    while (head < tail) {
        uint16_t here = queue[head++];
        uint8_t next_distance = (uint8_t) (distance[here] + 1U);
        for (uint8_t face = 0; face < FACES; ++face) {
            uint16_t next = here;
            for (uint8_t turn = 0; turn < 3; ++turn) {
                next = table[face][next];
                if (distance[next] == UINT8_MAX) {
                    distance[next] = next_distance;
                    queue[tail++] = next;
                    if (next_distance > *max_distance)
                        *max_distance = next_distance;
                }
            }
        }
    }
    free(queue);
    if (tail != PERMUTATIONS) {
        free(distance);
        return NULL;
    }
    return distance;
}

static uint8_t *bfs_orientation(uint16_t table[3][ORIENTATIONS],
                                uint8_t *max_distance)
{
    uint8_t *distance = malloc(ORIENTATIONS);
    uint16_t *queue = malloc((size_t) ORIENTATIONS * sizeof *queue);
    uint32_t head = 0, tail = 1;
    if (!distance || !queue) {
        free(distance);
        free(queue);
        return NULL;
    }
    memset(distance, UINT8_MAX, ORIENTATIONS);
    distance[0] = 0;
    queue[0] = 0;
    *max_distance = 0;
    while (head < tail) {
        uint16_t here = queue[head++];
        uint8_t next_distance = (uint8_t) (distance[here] + 1U);
        for (uint8_t face = 0; face < FACES; ++face) {
            uint16_t next = here;
            for (uint8_t turn = 0; turn < 3; ++turn) {
                next = table[face][next];
                if (distance[next] == UINT8_MAX) {
                    distance[next] = next_distance;
                    queue[tail++] = next;
                    if (next_distance > *max_distance)
                        *max_distance = next_distance;
                }
            }
        }
    }
    free(queue);
    if (tail != ORIENTATIONS) {
        free(distance);
        return NULL;
    }
    return distance;
}

static double wall_seconds(void)
{
#ifdef _WIN32
    LARGE_INTEGER counter, frequency;
    QueryPerformanceCounter(&counter);
    QueryPerformanceFrequency(&frequency);
    return (double) counter.QuadPart / (double) frequency.QuadPart;
#else
    return (double) time(NULL);
#endif
}

static uint8_t *build_oracle(uint16_t permutation[3][PERMUTATIONS],
                             uint16_t orientation[3][ORIENTATIONS])
{
    uint8_t *distance = malloc(STATES);
    uint32_t *queue = malloc((size_t) STATES * sizeof *queue);
    uint32_t head = 0, tail = 1;
    if (!distance || !queue) {
        free(distance);
        free(queue);
        return NULL;
    }
    memset(distance, UINT8_MAX, STATES);
    distance[0] = 0;
    queue[0] = 0;
    while (head < tail) {
        uint32_t here = queue[head++];
        uint16_t p = (uint16_t) (here / ORIENTATIONS);
        uint16_t o = (uint16_t) (here % ORIENTATIONS);
        uint8_t next_distance = (uint8_t) (distance[here] + 1U);
        for (uint8_t face = 0; face < FACES; ++face) {
            uint16_t next_p = p, next_o = o;
            for (uint8_t turn = 0; turn < 3; ++turn) {
                next_p = permutation[face][next_p];
                next_o = orientation[face][next_o];
                uint32_t there = (uint32_t) next_p * ORIENTATIONS + next_o;
                if (distance[there] == UINT8_MAX) {
                    distance[there] = next_distance;
                    queue[tail++] = there;
                }
            }
        }
    }
    free(queue);
    if (tail != STATES) {
        free(distance);
        return NULL;
    }
    return distance;
}

static int write_pdb(const char *path, const uint8_t *permutation,
                     const uint8_t *orientation)
{
    FILE *file = fopen(path, "w");
    if (!file)
        return 0;
    fputs("#ifndef STAGE3_PDB_H\n#define STAGE3_PDB_H\n\n"
          "#include <stdint.h>\n\n"
          "enum { STAGE3_PDB_PERMUTATIONS = 5040, "
          "STAGE3_PDB_ORIENTATIONS = 729 };\n\n", file);
    fputs("static const uint8_t stage3_permutation_pdb["
          "STAGE3_PDB_PERMUTATIONS] = {\n", file);
    for (uint16_t i = 0; i < PERMUTATIONS; ++i) {
        if (i % 24U == 0)
            fputs("    ", file);
        fprintf(file, "%u%s", permutation[i], i + 1U == PERMUTATIONS ? "" : ", ");
        if (i % 24U == 23U || i + 1U == PERMUTATIONS)
            fputc('\n', file);
    }
    fputs("};\n\nstatic const uint8_t stage3_orientation_pdb["
          "STAGE3_PDB_ORIENTATIONS] = {\n", file);
    for (uint16_t i = 0; i < ORIENTATIONS; ++i) {
        if (i % 24U == 0)
            fputs("    ", file);
        fprintf(file, "%u%s", orientation[i], i + 1U == ORIENTATIONS ? "" : ", ");
        if (i % 24U == 23U || i + 1U == ORIENTATIONS)
            fputc('\n', file);
    }
    fputs("};\n\n#endif\n", file);
    return fclose(file) == 0;
}

static int verify_h1(const uint8_t *permutation, const uint8_t *orientation,
                     const uint8_t *oracle, uint32_t *failures)
{
    *failures = 0;
    for (uint32_t rank = 0; rank < STATES; ++rank) {
        uint16_t p = (uint16_t) (rank / ORIENTATIONS);
        uint16_t o = (uint16_t) (rank % ORIENTATIONS);
        uint8_t h = permutation[p] > orientation[o] ? permutation[p] : orientation[o];
        if (h > oracle[rank])
            ++*failures;
    }
    return *failures == 0;
}

static int verify_complete(const uint8_t *distance, uint32_t count,
                           uint8_t *maximum)
{
    *maximum = 0;
    for (uint32_t i = 0; i < count; ++i) {
        if (distance[i] == UINT8_MAX)
            return 0;
        if (distance[i] > *maximum)
            *maximum = distance[i];
    }
    return 1;
}

int main(int argc, char **argv)
{
    const char *output = argc == 2 ? argv[1] : "stage3-pdb.h";
    uint16_t permutation[3][PERMUTATIONS];
    uint16_t orientation[3][ORIENTATIONS];
    uint8_t *permutation_pdb;
    uint8_t *orientation_pdb;
    uint8_t *oracle;
    uint8_t permutation_max, orientation_max, oracle_max;
    uint32_t failures;
    int h1_pass;
    double start, finish;

    memcpy(permutation, stage3_permutation_transition, sizeof permutation);
    memcpy(orientation, stage3_orientation_transition, sizeof orientation);
    permutation_pdb = bfs_pdb(permutation, &permutation_max);
    orientation_pdb = bfs_orientation(orientation, &orientation_max);
    if (!permutation_pdb || !orientation_pdb) {
        fputs("PDB BFS failed or incomplete\n", stderr);
        free(permutation_pdb);
        free(orientation_pdb);
        return EXIT_FAILURE;
    }
    if (!verify_complete(permutation_pdb, PERMUTATIONS, &permutation_max) ||
        !verify_complete(orientation_pdb, ORIENTATIONS, &orientation_max)) {
        fputs("PDB completeness: FAIL\n", stderr);
        free(permutation_pdb);
        free(orientation_pdb);
        return EXIT_FAILURE;
    }
    if (!write_pdb(output, permutation_pdb, orientation_pdb)) {
        fprintf(stderr, "could not write %s\n", output);
        free(permutation_pdb);
        free(orientation_pdb);
        return EXIT_FAILURE;
    }
    printf("PDB completeness: PASS (permutation=%u, orientation=%u)\n",
           PERMUTATIONS, ORIENTATIONS);
    printf("permutation_max_distance: %u\n", permutation_max);
    printf("permutation_solved_entry: %u\n", permutation_pdb[0]);
    printf("orientation_max_distance: %u\n", orientation_max);
    printf("orientation_solved_entry: %u\n", orientation_pdb[0]);
    printf("pdb_data_bytes: %zu\n",
           (size_t) PERMUTATIONS + ORIENTATIONS);

    start = wall_seconds();
    oracle = build_oracle(permutation, orientation);
    if (!oracle) {
        fputs("oracle BFS failed or incomplete\n", stderr);
        free(permutation_pdb);
        free(orientation_pdb);
        return EXIT_FAILURE;
    }
    if (!verify_complete(oracle, STATES, &oracle_max)) {
        fputs("oracle completeness: FAIL\n", stderr);
        free(oracle);
        free(permutation_pdb);
        free(orientation_pdb);
        return EXIT_FAILURE;
    }
    printf("oracle completeness: PASS (checked=%u, max_distance=%u)\n",
           STATES, oracle_max);
        h1_pass = verify_h1(permutation_pdb, orientation_pdb, oracle, &failures);
        finish = wall_seconds();
        printf("H1: %s (checked=%u, failures=%u)\n",
            h1_pass ? "PASS" : "FAIL", STATES, failures);
    printf("H1_wall_clock_seconds: %.6f\n", finish - start);

    free(oracle);
    free(permutation_pdb);
    free(orientation_pdb);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
