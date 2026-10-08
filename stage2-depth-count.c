/* Stage 2：統計 HTM 深度 0 到 5 的狀態數量*/
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
	MAX_DEPTH = 6
};

typedef struct {
	uint8_t p[CUBIES], o[CUBIES];
} state_t;

/* 每個目標位置的方塊來自 source[face][destination] */
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

static int build_transitions(uint16_t permutation[3][PERMUTATIONS],
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
	return 0;
}

int main(void)
{
	uint16_t permutation[3][PERMUTATIONS], orientation[3][ORIENTATIONS];
	uint8_t *visited = malloc(STATES);
	uint32_t *frontier = malloc((size_t) STATES * sizeof *frontier);
	uint32_t *next_frontier = malloc((size_t) STATES * sizeof *next_frontier);
	uint32_t frontier_count = 1;
	uint32_t cumulative = 1;

	if (!visited || !frontier || !next_frontier) {
		fprintf(stderr, "allocation failed\n");
		free(visited);
		free(frontier);
		free(next_frontier);
		return EXIT_FAILURE;
	}

	build_transitions(permutation, orientation);
	memset(visited, 0, STATES);
	visited[0] = 1;
	frontier[0] = 0;

	printf("depth  new_states  cumulative_states\n");
	printf("%5d  %11u  %17u\n", 0, 1U, cumulative);
	for (uint8_t depth = 1; depth <= MAX_DEPTH; ++depth) {
		uint32_t next_count = 0;
		for (uint32_t i = 0; i < frontier_count; ++i) {
			uint32_t here = frontier[i];
			uint16_t p = (uint16_t) (here / ORIENTATIONS);
			uint16_t o = (uint16_t) (here % ORIENTATIONS);
			for (uint8_t face = 0; face < 3; ++face) {
				uint16_t next_p = p, next_o = o;
				for (uint8_t turn = 0; turn < 3; ++turn) {
					next_p = permutation[face][next_p];
					next_o = orientation[face][next_o];
					uint32_t there =
						(uint32_t) next_p * ORIENTATIONS + next_o;
					if (!visited[there]) {
						visited[there] = 1;
						next_frontier[next_count++] = there;
					}
				}
			}
		}
		cumulative += next_count;
		printf("%5u  %11u  %17u\n", depth, next_count, cumulative);
		uint32_t *temporary = frontier;
		frontier = next_frontier;
		next_frontier = temporary;
		frontier_count = next_count;
	}

	free(visited);
	free(frontier);
	free(next_frontier);
	return EXIT_SUCCESS;
}
