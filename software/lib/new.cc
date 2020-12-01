#include "new.h"

void* operator new(const size_t n, void* p)
{
  return p;
}
