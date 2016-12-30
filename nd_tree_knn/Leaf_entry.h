#include "config.h"

#ifndef LEAF_ENTRY_H
#define LEAF_ENTRY_H

struct Leaf_entry
{
    unsigned char key[DIM]; // Key: a string of DIM characters, represented as 0 to (alphabet_size-1)
    ND_tree_record record; // Actual data, usually a pointer
    void printKey() {
        for ( int i = 0 ; i < DIM ; i++ ) 
            printf("%c ",key[i]);
        puts("");
    }
};

#endif
