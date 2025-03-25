// this file is a copy of the graph generator written by the previous person to work on this project
// I've slightly modified it to work with my infrastructure here (output a graph object instead of writing to a file)
// I've also made it a header file instead of a regular program which command-line arguments are passed to

#ifndef __GRAPH_GENERATOR_HXX__
#define __GRAPH_GENERATOR_HXX__

#include <stdio.h>
#include <stdlib.h>
#include "graph.hxx"

//params: nC lC nK lK three_edges seed
graph generate_graph(long nC, long lC, long nK, long lK, long three_edges, long seed = -1) {
   if (seed == -1) seed = time(0);
   srand(seed);

   long n=nC*lC+nK*lK;
   long m=nC*lC+nK*(lK*(lK-1))/2+(2+three_edges)*(nC+nK-1);
  
   //shuffle nodes
   long* nodes = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){nodes[i]=i;}
   for(long i=0;i<n;i++)
   {
      long j=i+rand()%(n-i);
      long temp=nodes[i];
      nodes[i]=nodes[j];
      nodes[j]=temp;
   }   

   //shuffle types
   char* graph_type = (char*)malloc(sizeof(char)*(nC+nK));
   long indx=0;
   for(long i=0;i<nC;i++){graph_type[indx++]=0;}
   for(long i=0;i<nK;i++){graph_type[indx++]=1;}
   for(long i=0;i<nC+nK;i++)
   {
      long j=i+rand()%(nC+nK-i);
      char temp=graph_type[i];
      graph_type[i]=graph_type[j];
      graph_type[j]=temp;
   }

   //create subgraph edges
   long* edges = (long*)malloc(sizeof(long)*m*2);
   long edge_indx=0;
   long* startNode = (long*)malloc(sizeof(long)*(nC+nK));
   long currentNode=0;
   for(long i=0;i<nC+nK;i++)
   {
      startNode[i]=currentNode;
      if(graph_type[i]==0)
      {
         for(long j=0;j<lC;j++)
         {
            edges[edge_indx++]=nodes[currentNode+j];
            edges[edge_indx++]=nodes[currentNode+(j+1)%lC];
         }
         currentNode+=lC;  
      }
      else
      {
         for(long j=0;j<lK;j++)
         {
            for(long k=j+1;k<lK;k++)
            {
               edges[edge_indx++]=nodes[currentNode+j];
               edges[edge_indx++]=nodes[currentNode+k];
            }
         }
         currentNode+=lK; 
      }
   }

   //connect the subgraphs in a tree structure
   for(long i=1;i<nC+nK;i++)
   {
      long j=rand()%i;
      long mod1=lC, mod2=lC;
      if(graph_type[i]==1){mod1=lK;}
      if(graph_type[j]==1){mod2=lK;}
      if(!three_edges)
      {
         long x1,y1,x2,y2;
         x1=rand()%mod1;
         x2=(x1+(1+rand()%(mod1-2)))%mod1;
         y1=rand()%mod2;
         y2=(y1+(1+rand()%(mod2-2)))%mod2;
         edges[edge_indx++]=nodes[startNode[i]+x1]; edges[edge_indx++]=nodes[startNode[j]+y1]; 
         edges[edge_indx++]=nodes[startNode[i]+x2]; edges[edge_indx++]=nodes[startNode[j]+y2]; 
      }
      else
      {
         long x1,y1,x2,y2,x3,y3;
         if(mod1==3){x1=0;x2=1;x3=2;}
         else
         {
            x1=rand()%mod1;
            x2=(x1+(2+rand()%(mod1-3)))%mod1;
            x3=(x1+(1+rand()%((mod1+x2-x1-1)%mod1)))%mod1;
         }
         if(mod2==3){y1=0;y2=1;y3=2;}
         else
         {
            y1=rand()%mod2;
            y2=(y1+(2+rand()%(mod2-3)))%mod2;
            y3=(y1+(1+rand()%((mod2+y2-y1-1)%mod2)))%mod2;
         }
         edges[edge_indx++]=nodes[startNode[i]+x1]; edges[edge_indx++]=nodes[startNode[j]+y1]; 
         edges[edge_indx++]=nodes[startNode[i]+x2]; edges[edge_indx++]=nodes[startNode[j]+y2]; 
         edges[edge_indx++]=nodes[startNode[i]+x3]; edges[edge_indx++]=nodes[startNode[j]+y3]; 
      }
   }

   //shuffle edges
   for(long i=0;i<m;i++)
   {
      long j=i+rand()%(m-i);
      long tempx=edges[2*i]; long tempy=edges[2*i+1];
      edges[2*i]=edges[2*j]; edges[2*i+1]=edges[2*j+1];
      if(rand()%2==0)
      {
         edges[2*j]=tempx; edges[2*j+1]=tempy;   
      }
      else
      {
         edges[2*j]=tempy; edges[2*j+1]=tempx;   
      }
   }

   // create a graph object (this is the main difference between the previous version and my version)
   // my version uses ints instead of longs (we wouldn't have enough memory to run the implementation on a graph with more than two billion vertices anyway, we use ~40 bytes an edge without even including the certificates), so I cast to ints
   graph retval;
   retval.n = n;
   retval.e = m;
   for (long i=0;i<n;i++)
   {
   	  retval.adjLists.emplace_back();
   }

   for(long i=0;i<m;i++)
   {
      retval.add_edge((int)(edges[2*i]), (int)(edges[2*i+1]));
   }

   for (long i=0;i<n;i++)
   {
   	  retval.adjLists[i].shrink_to_fit();
   }

   free(nodes);
   free(graph_type);
   free(edges);
   free(startNode);
   return retval;
}

#endif