#include <iostream>

#include <thread>
#include <chrono>
#include <atomic>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <SDL2/SDL.h>

#include "cell.h"


int main() {
  CellList cl;

  Cont<Vect> pp = {Vect(0)};

  cl.Reinit(pp, MIdx(3), Vect(0), Vect(1));

  return 0;
}
