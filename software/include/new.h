#ifndef __NEW_H
#define __NEW_H

typedef unsigned long	size_t;

void* operator new(const size_t n, void* p);

#endif
