
// This is a direct port of the C++ RFTL test program.
//

#include "stdafx.h"
#include <stdio.h>
#include "RTree.h"
#include "Records.h"
#include <iostream>
#include <time.h>
using namespace std;

int nrects = sizeof(rects) / sizeof(rects[0]);

RTree<int, int, 2, float> tree;

bool MySearchCallback(int id, void* arg) 
{
  printf("Hit data rect %d\n", id);
  return true; // keep going
}

void test1();
void test2();

static double diffclock(clock_t clock1,clock_t clock2)
{
    double diffticks=clock1-clock2;
    double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
    return diffms;
}
void main()
{
 
  ofstream myoutput;
  myoutput.open("RFTL_Output");   //write RFTL results 

  ofstream Rtreeoutput;
  Rtreeoutput.open("RTree_Output");  //write original Rtree results

  //printf("nrects = %d\n", nrects);
  clock_t t1 = 0;

  printf("Rtree is bulding...\n");

  test1();
 
  int itIndex = 0;
  RTree<int, int, 2, float>::Iterator it;
  for( tree.GetFirst(it); 
       !tree.IsNull(it); 
       tree.GetNext(it) )
  {
    int value = tree.GetAt(it);
    
    int boundsMin[2] = {0,0};
    int boundsMax[2] = {0,0};
    it.GetBounds(boundsMin, boundsMax);
   // printf("it[%d] %d = (%d,%d,%d,%d)\n", itIndex++, value, boundsMin[0], boundsMin[1], boundsMax[0], boundsMax[1]);
	Rtreeoutput << value <<" = (" << boundsMin[0] << ", " << boundsMin[1]<< ", " << boundsMax[0] <<", " << boundsMax[1] << " ) \n";
	
  }
  
 
  tree.Save("rtreefile");
  clock_t t2 = clock();

  printf("time for Rtree constraction in ms is %f\n",diffclock(t2,t1));
  Rtreeoutput.close();


  RTree<int, int, 2, float> tree1;

  tree1.Load("rtreefile");
  for( tree1.GetFirst(it); 
       !tree1.IsNull(it);
       tree1.GetNext(it) )
  {
    int value = tree.GetAt(it);
    
    int boundsMin[2] = {0,0};
    int boundsMax[2] = {0,0};
    it.GetBounds(boundsMin, boundsMax);
  //  printf("it[%d] %d = (%d,%d,%d,%d)\n", itIndex++, value, boundsMin[0], boundsMin[1], boundsMax[0], boundsMax[1]);
	myoutput << value <<" = (" << boundsMin[0] << ", " << boundsMin[1]<< ", " << boundsMax[0] <<", " << boundsMax[1] << " ) \n";
  }
  

  myoutput.close();


  printf("press a button to exit\n");
  getchar(); 


  }

// Build Rtree with 7000 elements
void test1()
{

	  for(int i=0; i<nrects; i++)
  {
    tree.Insert(rects[i].min, rects[i].max, i); 
  }

}

// Make 1090 insertions and 400 deletions
void test2()
{

 for(int i=0; i<1000; i++)
  {
	 tree.Insert(rects[i].min, rects[i].max, i); 
  }
	  for(int i=0; i<10; i++)
  {
	 tree.Remove(rects[i].min, rects[i].max, i); 
  }
  
 for(int i=20; i<30; i++)
  {
	 tree.Remove(rects[i].min, rects[i].max, i); 
  }

 for(int i=200; i<300; i++)
  {
	 tree.Remove(rects[i].min, rects[i].max, i); 
  }

 for(int i=600; i<800; i++)
  {
	 tree.Remove(rects[i].min, rects[i].max, i); 
  }
	
 for(int i=0; i<10; i++)
  {
    tree.Insert(rects[i].min, rects[i].max, i); 
  }

 for(int i=820; i<900; i++)
  {
    tree.Remove(rects[i].min, rects[i].max, i); 
  }
 for(int i=700; i<780; i++)
  {
    tree.Insert(rects[i].min, rects[i].max, i); 
  }

}
