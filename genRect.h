#ifndef GENRECT_H
#define GENRECT_H

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "RTree.h"

struct Rect generateRect(int grid,int ID){

int maxLength = grid*5/100;    ///5% of grid length
int allowed = 0;
int rX;
int rY;
struct Rect rect;

int rWidth= 1+rand()%maxLength;
int rHeight = 1+rand()%maxLength;

while (allowed == 0){
    rX = 1+rand()%grid;
    rY = 1+rand()%grid;

    if ( (rX+rWidth < grid) && (rY+rHeight < grid) )
        {allowed = 1;}
}

rect.x2 = rX+rWidth;
rect.x1 = rX;

rect.y2 = rY+rHeight;
rect.y1 = rY;

rect.rectID = ID;
rect.childNode = NULL;

return (rect);

}

struct Rect ownRect(int x1,int x2, int y1, int y2,int ID){
    struct Rect myRect;

    myRect.x1 = x1;
    myRect.x2 = x2;
    myRect.y1 = y1;
    myRect.y2 = y2;
    myRect.rectID = ID;
    myRect.childNode = NULL;

    return (myRect);
}

#endif

