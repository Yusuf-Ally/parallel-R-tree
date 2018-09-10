#ifndef RTREE_H
#define RTREE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include "queue.h"

#define PGSIZE 512
#define GRID_MAX 1000000
#define NUMRECTS 1001
#define M 50   //Max. entries per node
#define THREADS 2
#define null -1
int level = 0;  //Begin at root node

typedef int bool;
#define true 1
#define false 0

int numberSlices;

/* Helper functions*/

/// Get number of leaf nodes
    int numLeafs(){
        int result = NUMRECTS/M;
        if (NUMRECTS%M != 0)
            result++;

        return result;
        }

/// Get number of slices
    int numSlices(){
        double result = ceil(sqrt(numLeafs()));
        return result;
        }

/// Get number of rects per slice
    int rectSlices(){
        return numSlices()*M;
        }

/// Get height of tree, with the height of the root being 1
    int treeHeight(){
        return ceil(log(NUMRECTS)/log(M));
        }

/// Calculate max. number of nodes
    int numNodes(){
        double N = NUMRECTS;
        double capacity = M;
        int result = 0;

        for (int i=1;i<=treeHeight();i++){
            result += ceil(N/pow(capacity,i));
        }
        return result++;
    }

/*Structure definitions*/

struct Rect{
int x1;
int x2;
int y1;
int y2;
int rectID;
struct Node *childNode;
};

struct Node {
int level;      /// Level of node; root is 1
int nodeID;     /// Same as page number
int numEntries;
int sliceNumber;
struct Rect nodeMBR;        ///Minimum bounding rectangle of entire node ie. max of all entries MBRs
struct Rect entries[M];        ///Each node holds up to M entries,empty entry is NULL
};

struct Rect * genRectList(int numRects){
    struct Rect *arrPtr;
    arrPtr = (struct Rect*)malloc(numRects*sizeof(struct Rect));
        if (arrPtr == NULL){
            printf("Error allocating memory.");
            exit(0);}
    return arrPtr;
        }


int randnorm()
{
    int r = 0;

    for (int i = 0; i < 20; i++)
        r += rand();

    return r;
}

int __r_srand_done = 0;

struct Rect * genClusteredRectList(int numRects, int numClusters, int gridsize){
    struct Rect *arrPtr;
    arrPtr = (struct Rect*)malloc(numRects*sizeof(struct Rect));
        if (arrPtr == NULL){
            printf("Error allocating memory.");
            exit(0);}

    if (!__r_srand_done)
    {
        srand((unsigned int)time(NULL));
        __r_srand_done = 1;
    }

    int numofdims = 2;
    int rectsize = gridsize * 5.0 / 100.0;
    int valuespercluster = numRects / numClusters;

    double dv = (double)RAND_MAX / gridsize;

    int ccentre[2];

    int nval = 0;

    for (int d = 0; d < numClusters; d++)
    {
        ccentre[0] = rand() % gridsize;
        ccentre[1] = rand() % gridsize;

        for (int i = 0; i < valuespercluster; i++)
        {
            
            arrPtr[nval].x1 = ccentre[0] + (randnorm() / dv - gridsize / 2);

            if (arrPtr[nval].x1 < 0.0 )
                 arrPtr[nval].x1 = 0.0;
            
            arrPtr[nval].x2 = arrPtr[nval].x1 +  rand() % rectsize;
            
            arrPtr[nval].y1 = ccentre[0] + (randnorm() / dv - gridsize / 2);

            if (arrPtr[nval].y1 < 0.0 )
                 arrPtr[nval].y1 = 0.0;
            
            arrPtr[nval].y2 = arrPtr[nval].y1 +  rand() % rectsize;
            
            arrPtr[nval].rectID = nval;
            arrPtr[nval].childNode = NULL;
            nval += 1;
        }
    }

    for (int i = nval; i < numRects; i++)
    {
        for (int j = 0; j < numofdims; j++)
        {
            arrPtr[nval].x1 = rand() % gridsize;
            arrPtr[nval].x2 = arrPtr[nval].x1 + rand() % rectsize;
            arrPtr[nval].y1 = rand() % gridsize;
            arrPtr[nval].y2 = arrPtr[nval].y1 + rand() % rectsize;
            arrPtr[nval].rectID = i;
            arrPtr[nval].childNode = NULL;
        }
    }


    return arrPtr;
}


/// Sort values of rect array by ascending x value
void sortX(struct Rect *arr){
struct Rect temp;

        for (int i=0;i<NUMRECTS-1;i++){
            for (int j=i+1;j<NUMRECTS;j++){
                if (arr[i].x1 > arr[j].x1){
                    temp = arr[i];
                    arr[i] = arr[j];
                    arr[j] = temp;}
                }
        }
    }

// Radix sort implementation, adapted from https://www.geeksforgeeks.com/radix-sort/
int getMax(struct Rect *arr, int n)
{
    int mx = arr[0].x1;
    for (int i = 1; i < n; i++)
        if (arr[i].x1 > mx)
            mx = arr[i].x1;
    return mx;
}
 
// A function to do counting sort of arr[] according to
// the digit represented by exp.
void countSort(struct Rect *arr, int n, int exp)
{
    int *output = (int*)malloc(n*sizeof(int)); // output array
    int i, count[10] = {0};
 
    // Store count of occurrences in count[]
    for (i = 0; i < n; i++)
        count[ (arr[i].x1/exp)%10 ]++;
 
    // Change count[i] so that count[i] now contains actual
    //  position of this digit in output[]
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];
 
    // Build the output array
    for (i = n - 1; i >= 0; i--)
    {
        output[count[ (arr[i].x1/exp)%10 ] - 1] = arr[i].x1;
        count[ (arr[i].x1/exp)%10 ]--;
    }
 
    // Copy the output array to arr[], so that arr[] now
    // contains sorted numbers according to current digit
    for (i = 0; i < n; i++)
        arr[i].x1 = output[i];
}
 
// The main function to that sorts arr[] of size n using 
// Radix Sort
void radixSort(struct Rect *arr, int n)
{
    // Find the maximum number to know number of digits
    int m = getMax(arr, n);
 
    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    for (int exp = 1; m/exp > 0; exp *= 10)
        countSort(arr, n, exp);
}

void sortY(struct Rect *arr){
struct Rect temp;
int numberSlices = numSlices();
int lastEntries = NUMRECTS%rectSlices(); /// Number of entries in  last slice

if (rectSlices()*numberSlices == NUMRECTS){
    for (int k=0;k<numberSlices ;k++){
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
else{

    if (rectSlices()*numberSlices > (NUMRECTS + rectSlices())){
        numberSlices--;
    }

    for (int k=0;k<numberSlices ;k++){

            if (k < (numberSlices-1)){
        for (int i=k*rectSlices();i<((k+1)*rectSlices())-1;i++){
            for (int j=i+1;j<(k+1)*rectSlices();j++){
                if (arr[i].y1 > arr[j].y1){
                    temp = arr[i];
                    arr[i] = arr[j];
                    arr[j] = temp;}
                }
        }
            }
            if (k == (numberSlices-1)){
            for (int i=k*rectSlices();i< (k*rectSlices()) + (lastEntries)-1;i++){
                for (int j=i+1;j<(k*rectSlices()) + (lastEntries);j++){
                    if (arr[i].y1 > arr[j].y1){
                        temp = arr[i];
                        arr[i] = arr[j];
                        arr[j] = temp;}
                }
        }
            }
    }
    }

}

/// Build RTree and return pointer to root node
struct Node* buildTree(struct Rect* rectArray){
    int currentNumSlices = numSlices();
    double leaves = numLeafs();
    int rectIndex = 0;
    int nodeIndex = 0;
    int nodes = 0;

        if (rectSlices()*currentNumSlices > (NUMRECTS + rectSlices())){
        currentNumSlices--;
    }

/// Allocate memory for all nodes
    struct Node* nodeArray  = (struct Node*)malloc(numNodes()*sizeof(struct Node));
    if (nodeArray == NULL){
        printf("Error allocating memory.");
        exit(0);}

/// First create leaf nodes, at level being equal to tree height
for (int Level=treeHeight();Level > 0;Level--){

int nodesInLevel = ceil(leaves/pow(M,treeHeight()-Level));

/// No extra entries case
        /// If at the leaf level
        if (Level == treeHeight()){
            for (int leaf = 0; leaf < numLeafs(); leaf++){

                if (leaf != (numLeafs()-1) ){
                for (int entry = 0; entry < M; entry++){
                    nodeArray[leaf].level = Level;
                    nodeArray[leaf].nodeID = nodes;
                    nodeArray[leaf].sliceNumber = rectIndex/rectSlices();
                    nodeArray[leaf].entries[entry] = rectArray[rectIndex];
                    nodeArray[leaf].numEntries = entry+1;

                    /// Set node MBR
                    if (entry == 0){
                        nodeArray[leaf].nodeMBR.x1 = rectArray[rectIndex].x1;
                        nodeArray[leaf].nodeMBR.x2 = rectArray[rectIndex].x2;
                        nodeArray[leaf].nodeMBR.y1 = rectArray[rectIndex].y1;
                        nodeArray[leaf].nodeMBR.y2 = rectArray[rectIndex].y2;
;                    }
                    else{
                        if (rectArray[rectIndex].x1 < nodeArray[leaf].nodeMBR.x1)
                            nodeArray[leaf].nodeMBR.x1 = rectArray[rectIndex].x1;
                        if (rectArray[rectIndex].x2 > nodeArray[leaf].nodeMBR.x2)
                            nodeArray[leaf].nodeMBR.x2 = rectArray[rectIndex].x2;
                        if (rectArray[rectIndex].y1 < nodeArray[leaf].nodeMBR.y1)
                            nodeArray[leaf].nodeMBR.y1 = rectArray[rectIndex].y1;
                        if (rectArray[rectIndex].y2 > nodeArray[leaf].nodeMBR.y2)
                            nodeArray[leaf].nodeMBR.y2 = rectArray[rectIndex].y2;
                    }
                    rectIndex++;
                }
                nodes++;
            }

                if (leaf == (numLeafs()-1) ){
                int lastEntries = NUMRECTS % M;
                if (lastEntries == 0){
                    lastEntries = M;
                }

                for (int entry = 0; entry < lastEntries; entry++){
                    nodeArray[leaf].level = Level;
                    nodeArray[leaf].nodeID = nodes;
                    nodeArray[leaf].sliceNumber = rectIndex/rectSlices();
                    nodeArray[leaf].entries[entry] = rectArray[rectIndex];
                    nodeArray[leaf].numEntries = entry+1;

                    /// Set node MBR
                    if (entry == 0){
                        nodeArray[leaf].nodeMBR.x1 = rectArray[rectIndex].x1;
                        nodeArray[leaf].nodeMBR.x2 = rectArray[rectIndex].x2;
                        nodeArray[leaf].nodeMBR.y1 = rectArray[rectIndex].y1;
                        nodeArray[leaf].nodeMBR.y2 = rectArray[rectIndex].y2;
;                    }
                    else{
                        if (rectArray[rectIndex].x1 < nodeArray[leaf].nodeMBR.x1)
                            nodeArray[leaf].nodeMBR.x1 = rectArray[rectIndex].x1;
                        if (rectArray[rectIndex].x2 > nodeArray[leaf].nodeMBR.x2)
                            nodeArray[leaf].nodeMBR.x2 = rectArray[rectIndex].x2;
                        if (rectArray[rectIndex].y1 < nodeArray[leaf].nodeMBR.y1)
                            nodeArray[leaf].nodeMBR.y1 = rectArray[rectIndex].y1;
                        if (rectArray[rectIndex].y2 > nodeArray[leaf].nodeMBR.y2)
                            nodeArray[leaf].nodeMBR.y2 = rectArray[rectIndex].y2;
                    }
                    rectIndex++;
                }
                nodes++;
            }
            }

        }

        /// If not at leaf level

        if ( (Level < treeHeight()) && (Level > 1) ){
            int nodesInPrevLevel = ceil(leaves/pow(M,treeHeight()-Level-1));
            int remainingNodes = nodesInPrevLevel % M;

            for (int curr_node = 0;curr_node < nodesInLevel; curr_node++){
                if (curr_node != (nodesInLevel-1) ){
                    for (int entry = 0;entry < M;entry++){
                        //printf("entry %d of node %d points to node %d ===\n",entry,nodes,nodeIndex);
                        nodeArray[nodes].level = Level;
                        nodeArray[nodes].nodeID = nodes;
                        nodeArray[nodes].entries[entry].rectID = rectIndex;
                        nodeArray[nodes].numEntries = entry+1;
                        nodeArray[nodes].entries[entry].childNode = &nodeArray[nodeIndex];
                        nodeArray[nodes].entries[entry].x1 = nodeArray[nodeIndex].nodeMBR.x1;
                        nodeArray[nodes].entries[entry].x2 = nodeArray[nodeIndex].nodeMBR.x2;
                        nodeArray[nodes].entries[entry].y1 = nodeArray[nodeIndex].nodeMBR.y1;
                        nodeArray[nodes].entries[entry].y2 = nodeArray[nodeIndex].nodeMBR.y2;

                    /// Set node MBR
                    if (entry == 0){
                        nodeArray[nodes].nodeMBR.x1 = nodeArray[nodes].entries[entry].x1;
                        nodeArray[nodes].nodeMBR.x2 = nodeArray[nodes].entries[entry].x2;
                        nodeArray[nodes].nodeMBR.y1 = nodeArray[nodes].entries[entry].y1;
                        nodeArray[nodes].nodeMBR.y2 = nodeArray[nodes].entries[entry].y2;
;                    }
                    else{
                        if (nodeArray[nodes].entries[entry].x1 < nodeArray[nodes].nodeMBR.x1)
                            nodeArray[nodes].nodeMBR.x1 = nodeArray[nodes].entries[entry].x1;
                        if (nodeArray[nodes].entries[entry].x2 > nodeArray[nodes].nodeMBR.x2)
                            nodeArray[nodes].nodeMBR.x2 = nodeArray[nodes].entries[entry].x2;
                        if (nodeArray[nodes].entries[entry].y1 < nodeArray[nodes].nodeMBR.y1)
                            nodeArray[nodes].nodeMBR.y1 = nodeArray[nodes].entries[entry].y1;
                        if (nodeArray[nodes].entries[entry].y2 > nodeArray[nodes].nodeMBR.y2)
                            nodeArray[nodes].nodeMBR.y2 = nodeArray[nodes].entries[entry].y2;
                    }

                        rectIndex++;
                        nodeIndex++;

                    }
                    nodes++;
                }

                if (curr_node == (nodesInLevel-1) ){
                    for (int entry = 0;entry < remainingNodes; entry++){
                        //printf("entry %d of node %d points to node %d ===\n",entry,nodes,nodeIndex);
                        nodeArray[nodes].level = Level;
                        nodeArray[nodes].nodeID = nodes;
                        nodeArray[nodes].entries[entry].rectID = rectIndex;
                        nodeArray[nodes].numEntries = entry+1;
                        nodeArray[nodes].entries[entry].childNode = &nodeArray[nodeIndex];
                        nodeArray[nodes].entries[entry].x1 = nodeArray[nodeIndex].nodeMBR.x1;
                        nodeArray[nodes].entries[entry].x2 = nodeArray[nodeIndex].nodeMBR.x2;
                        nodeArray[nodes].entries[entry].y1 = nodeArray[nodeIndex].nodeMBR.y1;
                        nodeArray[nodes].entries[entry].y2 = nodeArray[nodeIndex].nodeMBR.y2;

                    /// Set node MBR
                    if (entry == 0){
                        nodeArray[nodes].nodeMBR.x1 = nodeArray[nodes].entries[entry].x1;
                        nodeArray[nodes].nodeMBR.x2 = nodeArray[nodes].entries[entry].x2;
                        nodeArray[nodes].nodeMBR.y1 = nodeArray[nodes].entries[entry].y1;
                        nodeArray[nodes].nodeMBR.y2 = nodeArray[nodes].entries[entry].y2;
;                    }
                    else{
                        if (nodeArray[nodes].entries[entry].x1 < nodeArray[nodes].nodeMBR.x1)
                            nodeArray[nodes].nodeMBR.x1 = nodeArray[nodes].entries[entry].x1;
                        if (nodeArray[nodes].entries[entry].x2 > nodeArray[nodes].nodeMBR.x2)
                            nodeArray[nodes].nodeMBR.x2 = nodeArray[nodes].entries[entry].x2;
                        if (nodeArray[nodes].entries[entry].y1 < nodeArray[nodes].nodeMBR.y1)
                            nodeArray[nodes].nodeMBR.y1 = nodeArray[nodes].entries[entry].y1;
                        if (nodeArray[nodes].entries[entry].y2 > nodeArray[nodes].nodeMBR.y2)
                            nodeArray[nodes].nodeMBR.y2 = nodeArray[nodes].entries[entry].y2;
                    }

                        rectIndex++;
                        nodeIndex++;

                    }
                    nodes++;
                }
            }


        }
        /// If at root level ie. height = 1
        if (Level == 1){
            int nodesInPrevLevel = ceil(leaves/pow(M,treeHeight()-Level-1));
            int remainingNodes = nodesInPrevLevel % M;

                    for (int entry = 0;entry < remainingNodes; entry++){
                        //printf("entry %d of node %d points to node %d ===\n",entry,nodes,nodeIndex);
                        nodeArray[nodes].level = Level;
                        nodeArray[nodes].nodeID = nodes;
                        nodeArray[nodes].entries[entry].rectID = rectIndex;
                        nodeArray[nodes].numEntries = entry+1;
                        nodeArray[nodes].entries[entry].childNode = &nodeArray[nodeIndex];
                        nodeArray[nodes].entries[entry].x1 = nodeArray[nodeIndex].nodeMBR.x1;
                        nodeArray[nodes].entries[entry].x2 = nodeArray[nodeIndex].nodeMBR.x2;
                        nodeArray[nodes].entries[entry].y1 = nodeArray[nodeIndex].nodeMBR.y1;
                        nodeArray[nodes].entries[entry].y2 = nodeArray[nodeIndex].nodeMBR.y2;

                    /// Set node MBR
                    if (entry == 0){
                        nodeArray[nodes].nodeMBR.x1 = nodeArray[nodes].entries[entry].x1;
                        nodeArray[nodes].nodeMBR.x2 = nodeArray[nodes].entries[entry].x2;
                        nodeArray[nodes].nodeMBR.y1 = nodeArray[nodes].entries[entry].y1;
                        nodeArray[nodes].nodeMBR.y2 = nodeArray[nodes].entries[entry].y2;
;                    }
                    else{
                        if (nodeArray[nodes].entries[entry].x1 < nodeArray[nodes].nodeMBR.x1)
                            nodeArray[nodes].nodeMBR.x1 = nodeArray[nodes].entries[entry].x1;
                        if (nodeArray[nodes].entries[entry].x2 > nodeArray[nodes].nodeMBR.x2)
                            nodeArray[nodes].nodeMBR.x2 = nodeArray[nodes].entries[entry].x2;
                        if (nodeArray[nodes].entries[entry].y1 < nodeArray[nodes].nodeMBR.y1)
                            nodeArray[nodes].nodeMBR.y1 = nodeArray[nodes].entries[entry].y1;
                        if (nodeArray[nodes].entries[entry].y2 > nodeArray[nodes].nodeMBR.y2)
                            nodeArray[nodes].nodeMBR.y2 = nodeArray[nodes].entries[entry].y2;
                    }

                        rectIndex++;
                        nodeIndex++;

                    }
                    nodes++;
        }

}

    /// Return root node
    return (nodeArray+nodes-1);
}

int hitCount = 0;

int RangeSearch(struct Node* Root,struct Rect QueryWindow){
int entry;

    /// If not at the leaf level

    if (Root->level != treeHeight()){
        for ( entry = 0;entry < Root->numEntries; entry++){
             if ( QueryWindow.x2 >= Root->entries[entry].x1 &&
                  QueryWindow.x1 <= Root->entries[entry].x2 &&
                  QueryWindow.y2 >= Root->entries[entry].y1 &&
                  QueryWindow.y1 <= Root->entries[entry].y2){
                    // printf("\nNode %d : %dth entry points to node %d\n",Root->nodeID,entry,Root->entries[entry].childNode->nodeID);

                    RangeSearch(Root->entries[entry].childNode,QueryWindow);


                  }
        }
    }


    if (Root->level == treeHeight()){   /// If at leaf level ie. Level treeHeight

int ret = 0;
        for (int entry = 0;entry < Root->numEntries; entry++){
             if ( QueryWindow.x2 >= Root->entries[entry].x1 &&
                  QueryWindow.x1 <= Root->entries[entry].x2 &&
                  QueryWindow.y2 >= Root->entries[entry].y1 &&
                  QueryWindow.y1 <= Root->entries[entry].y2){
                    //printf("R%d ",Root->entries[entry].rectID);
                    ret++;

                  }
        }
        hitCount += ret;
        }

    return hitCount;


}

int RangeSearchParallel(struct Node* Root,struct Rect QueryWindow){
    #pragma omp parallel num_threads(4)
    {
        #pragma omp single
        RangeSearch(Root,QueryWindow);
    }
    return hitCount;
}

int RangeSearchQ(struct Node* Root,struct Rect QueryWindow){
    int nCount = 0;
    int tCount = 0; /// Number of idle threads

    queue *nodeQ = queue_create();

    queue_push(nodeQ,Root);

    #pragma omp parallel num_threads(THREADS) shared(nCount,tCount)
    {

        struct Node *node;
        #pragma omp critical
        node = queue_pop(nodeQ);

        if (node == NULL){
            #pragma omp ATOMIC
            tCount++;
            }

        int threadCount = omp_get_num_threads();
        while (tCount < threadCount){
            if (node != NULL){

                 if (node->level != treeHeight()){
        for (int entry = 0;entry < node->numEntries; entry++){
             if ( QueryWindow.x2 >= node->entries[entry].x1 &&
                  QueryWindow.x1 <= node->entries[entry].x2 &&
                  QueryWindow.y2 >= node->entries[entry].y1 &&
                  QueryWindow.y1 <= node->entries[entry].y2){
                    // printf("\nNode %d : %dth entry points to node %d\n",Root->nodeID,entry,Root->entries[entry].childNode->nodeID);
                    #pragma omp critical
                    queue_push(nodeQ,node->entries[entry].childNode);
                    //printf("pushing node %d\n",node->entries[entry].childNode->nodeID);


                  }
        }
    }

int ret = 0;
    if (node->level == treeHeight()){   /// If at leaf level ie. Level treeHeight

        for (int entry = 0;entry < node->numEntries; entry++){
             if ( QueryWindow.x2 >= node->entries[entry].x1 &&
                  QueryWindow.x1 <= node->entries[entry].x2 &&
                  QueryWindow.y2 >= node->entries[entry].y1 &&
                  QueryWindow.y1 <= node->entries[entry].y2){
                    //printf("R%d ",Root->entries[entry].rectID);
                    ret++;

                  }
        }
        #pragma omp ATOMIC
        nCount += ret;
        }

                } // End of node != NULL

            struct Node *nextNode;
            #pragma omp critical
            nextNode = queue_pop(nodeQ);

            if (nextNode == NULL && node != NULL){
                #pragma omp ATOMIC
                tCount++;
            }

            if (nextNode != NULL && node == NULL){
                #pragma omp ATOMIC
                tCount--;
            }

            node = nextNode;

            }

    }
    return nCount;
}



#endif

