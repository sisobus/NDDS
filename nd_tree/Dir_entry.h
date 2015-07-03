#ifndef DIR_ENTRY_H
#define DIR_ENTRY_H

#include "Leaf_node.h"

struct Dir_entry{
   unsigned char DMBR[DMBR_SIZE]; // Digital Minimum Bounding Rectangle
   unsigned int child; // Block number of the child of the current entry
};

#endif
