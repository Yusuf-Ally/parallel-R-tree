#ifndef RTREE_H
#define RTREE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PGSIZE 512
#define NUMRECTS 12
#define M 3 //Max. entries per node
#define m 1;  //Min. entries per node
int level = 0;  //Begin at root node

typedef int bool;
#define true 1
#define false 0


/// Get number of leaf nodes
    int numLeafs(){
        return ceil(NUMRECTS/M);
        }

/// Get number of slices
    int numSlices(){
        return ceil(sqrt(numLeafs()));
        }

/// Get number of rects per slice
    int rectSlices(){
        return numSlices()*M;
        }


struct Rect{
int x1;
int x2;
int y1;
int y2;
int rectID;
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

void insert(struct Entry newEntry,struct Node *Root){

    /* First invoke chooseLeaf - Find leaf requiring min perimeter increase to accomodate new entry.*/

    if (level == 0){  //If node is a leaf
        if (Root->numEntries < M){     //If node can accomodate entry
            Root->entries[Root->numEntries] = newEntry;
            Root->numEntries++;
        }

    }
}

/// Sort values of rect array by ascending x value
void sortX(struct Rect *arr){
struct Rect temp;

        for (int i=0;i<NUMRECTS-1;i++){
            for (int j=i+1;j<NUMRECTS;j++){
                if (arr[i].x1 > arr[j].x2){
                    temp = arr[i];
                    arr[i] = arr[j];
                    arr[j] = temp;}
                }
        }
    }

void sortY(struct Rect *arr){
struct Rect temp;
    for (int k=0;k<numSlices();k++){
        for (int i=k*rectSlices();i<((k+1)*rectSlices())-1;i++){
            for (int j=i+1;j<(k+1)*rectSlices();j++){
                if (arr[i].y1 > arr[j].y1){
                    temp = arr[i];
                    arr[i] = arr[j];
                    arr[j] = temp;}
                }
        }
    }
    }

#endif

