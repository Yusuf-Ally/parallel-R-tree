#include <stdio.h>
#include <stdlib.h>
#include "lsdtree.h"

int main()
{
    int numofdims = 4;
    int upper = 10000;
    srand((unsigned int)time(NULL));

    struct lsdtree * tree;

    tree = lsd_create(numofdims, 100, "testtree.lsdt");
    //tree = lsd_open("testtree.lsdt");

    printf("RNG started\n");

    int numofvalues = 1000000;
    num_t * arrvalues = malloc(sizeof(num_t)*numofvalues*numofdims);
    for (int i = 0; i < numofvalues*numofdims; i++)
    {
        arrvalues[i] = rand() % upper;
        //printf("%f\n", arrvalues[i]);
    }

    printf("inserted started\n");

    for (int i = 0; i < numofvalues; i++)
    {
        num_t * val = &arrvalues[i * numofdims];
        lsd_insert(tree,val);
    }

    printf("%d values inserted\n", numofvalues);

    lsd_free(tree);

    free(arrvalues);

    //lsd_print(tree);

    //lsd_sort(values, 1, numofvalues);
    //for (int i = 0; i < numofvalues; i++)
    //    printf("%f ",values[i][1]);

    return 0;
}
