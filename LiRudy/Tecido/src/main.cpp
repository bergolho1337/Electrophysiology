/*
    Author: Lucas Berg
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "../include/lirudy.h"
#include "../include/timer.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 < 3)
    {
        cout << "-------- LiRudy 2011 -------------" << endl;
        cout << argv[0] << " <dt> <tmax> <nthreads>" << endl;
        cout << "----------------------------------" << endl;
        exit (EXIT_FAILURE);
    }
    LiRudy *lirudy = new LiRudy(argc,argv);
    // DEBUG
    //cout << *lirudy << endl;
    
    double start, finish, elapsed;
    GET_TIME(start);

    lirudy->solve();
    
    GET_TIME(finish);
    elapsed = finish - start;
    printf("Time elapsed = %lf\n",elapsed);

    return 0;
}