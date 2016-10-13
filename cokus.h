
#pragma once
typedef unsigned long uint32;
#define N              (624)                 // length of state vector
#define M              (397)                 // a period parameter

extern uint32   state[N + 1];     // state vector + 1 extra to not violate ANSI C
extern uint32   *Nnext;          // next random value is computed from here
extern int      Lleft;


void seedMT(uint32 seed);
uint32 reloadMT(void);
uint32 randomMT(void);