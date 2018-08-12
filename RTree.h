#ifndef RTREE_H
#define RTREE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PGSIZE 512
#define M 4 //Max. entries per node
#define m 2;  //Min. entries per node
int level = 0;  //Begin at root node

typedef int bool;
#define true 1
#define false 0

struct Rect{
int x1;
int x2;
int y1;
int y2;
};

struct Entry{
struct Rect entryMBR;   ///Minimum bounding rectangle of entry
int *childNode;      ///Ptr to child node
};

struct Node {
bool isSplit;       ///Flag to indicate whether node has been split by insertion
int level;      /// Level of node; root is 0
int nodeID;     /// Same as page number
int numEntries;     ///Number of entries e[m,M]
struct Rect nodeMBR;        ///Minimum bounding rectangle of entire node ie. max of all entries MBRs
struct Entry entries[M];        ///Each node holds up to M entries,empty entry is NULL
//struct Node *child;
};

struct Node* Root(int x1,int x2, int y1, int y2)
{
  /// Allocate memory for new root
  struct Node* Root = (struct Node*)malloc(sizeof(struct Node));

  Root->nodeMBR.x1 = x1;
  Root->nodeMBR.x2 = x2;
  Root->nodeMBR.y1 = y1;
  Root->nodeMBR.y2 = y2;

  Root->level = 0;
  Root->isSplit = 0;
  Root->nodeID = 0;
  Root->numEntries = 0;

  /// Initialize child as NULL
  //Node->child = NULL;

  return(Root);
}

#endif
