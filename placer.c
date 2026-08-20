/////////////////////////////////////////////////////////////////////////////////////
// File       : placer.c
// Language   : C/C++
// Description: 8x8 (QP) Quadratic Placer
// Author     : Sandesh Prabhakar
// Created    : April 13, 2013
// Modified   : August 19, 2026
// Version    : 1.0
//
// Copyright (c) 2026 Sandesh Prabhakar
// All rights reserved
//////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

////////////////////////////////////////////////////////////////////////
// circuit class
////////////////////////////////////////////////////////////////////////

class circuit
{
    // circuit information

    // you may wish to add more private variables here

public:
    circuit(char *);            // constructor

    // you may wish to add more function declarations (or public variables) here
	int matrixSolver (float** , float** , float*  , int , int );	
	int verticalCut (float**, float**, float**, int**, int**, float*, int**, int, int*,int*, int*,float, float, float, float );
	int horizontalCut (float**, float**, float**, int**, int**, float*, int**, int, int*,int*, int*, float xmin, float ymin, float xmax, float ymax );
};
//void strcpy (char* , char* );
//char* strcat(char* , char const* );
////////////////////////////////////////////////////////////////////////
// constructor: reads in the *.lev file and builds basic data structure
// 		for the gate-level ckt
////////////////////////////////////////////////////////////////////////
circuit::circuit(char *cktName)
{
    // circuit information
    int **netlist;	// gate and net connectivity information
    int **fullnetlist;	// gate and net + pad and net connectivity information
    int **pads;		// pad and net connectivity information
	int **net;
    float **gate_xy;
    float **full_gate_xy;
	float **C_matrix;
	float **A_matrix;
	float *BX_matrix;
	float *BY_matrix;
	float *weight;

    FILE *inFile;
    FILE *outFile;
    char fName[40];
	char foName[40];
    int i, j, gcount, ncount, pcount;
    char c;
    int gatenum, padnum, junk;
    int f1, f2, f3;

    strcpy(fName, cktName);
    strcpy(foName, cktName);
    strcat(foName, ".3qp.output");

    inFile = fopen(fName, "r");
    if (inFile == NULL)
    {
	fprintf(stderr, "Can't open input file\n");
	exit(-1);
    }

    outFile = fopen(foName, "w");
    if (outFile == NULL)
    {
	fprintf(stderr, "Can't open 3qp output file to write\n");
	exit(-1);
    }

    fscanf(inFile, "%d", &gcount);	// number of gates
    fscanf(inFile, "%d", &ncount);	// number of nets 

    // allocate space for gates data structure
    netlist = new int * [gcount+1];
	netlist[0] = new int[2];
	netlist[0][0] = gcount;
	netlist[0][1] = ncount;

	net = new int * [15000];
	for (i=0; i<15000; i++)
    {
	net[i] = new int[15000];
	net[i][0]=0;
	}
    // now read in the netlist connectivity 
    for (i=1; i<=netlist[0][0]; i++)
    {
	fscanf(inFile, "%d", &gatenum);
	fscanf(inFile, "%d", &f1);

	f1 = (int) f1;

	netlist[gatenum] = new int[f1+1];
	netlist[gatenum][0] = f1;

		for (j=1; j<=netlist[gatenum][0]; j++)
		{
	    fscanf(inFile, "%d", &f1);
	    netlist[gatenum][j] = (int) f1;
		net[f1][0]=net[f1][0]+1;
		net[f1][net[f1][0]]=gatenum;
		}
	}

    fscanf(inFile, "%d", &pcount);	// number of pads
    pads = new int * [pcount+1];
	pads[0] = new int[1];
	pads[0][0] = pcount;

    // now read in the pad connectivity 
    for (i=1; i<=pads[0][0]; i++)
    {
	fscanf(inFile, "%d", &padnum);
	fscanf(inFile, "%d", &f1);
	fscanf(inFile, "%d", &f2);
	fscanf(inFile, "%d", &f3);

	f1 = (int) f1;
	f2 = (int) f2;
	f3 = (int) f3;

	pads[padnum] = new int[3];
	pads[padnum][0] = f1;
	net[f1][0]=net[f1][0]+1;
	net[f1][net[f1][0]]=gcount+padnum;

	pads[padnum][1] = f2; //X coordinate
	pads[padnum][2] = f3; //Y coordinate
	}

	// skip till end of line
	while ((c = getc(inFile)) != '\n' && c != EOF)
	    ;
//    }	// for (i...)
    fclose(inFile);

    // allocate space for fullnetlist data structure
    fullnetlist = new int * [15000];
	for (i=0; i<15000; i++)
    {
	fullnetlist[i] = new int[15000];
	}

    fullnetlist[0][0]=gcount;
    fullnetlist[0][1]=ncount;
    fullnetlist[0][2]=pcount;

	for (i=1; i<=gcount; i++)
    {
	fullnetlist[i][0]=netlist[i][0];
	for (j=1; j<=netlist[i][0]; j++)
		{
		fullnetlist[i][j] = netlist[i][j];
		}
	}

	for (i=1; i<=pcount; i++)
    {
	fullnetlist[gcount+i][0]=1;
	fullnetlist[gcount+i][1] = pads[i][0];
	}

    printf("Successfully read in netlist:\n");
    printf("\t%d Number of gates.\n", netlist[0][0]);
    printf("\t%d Number of nets. \n", netlist[0][1]);
    printf("\t%d Number of pads. \n", pads[0][0]);

    printf("Successfully read in FULL netlist:\n");
    printf("\t%d Number of gates.\n", fullnetlist[0][0]);
    printf("\t%d Number of nets. \n", fullnetlist[0][1]);
    printf("\t%d Number of pads. \n", fullnetlist[0][2]);

    printf("%d %d\n", fullnetlist[0][0], fullnetlist[0][1]);

//Clique formation
    weight = new float[15000];
    for (i=1; i<=netlist[0][1]; i++)
    {
	   weight[i]=1;
	}

float wt;
int p,q, r, s,x,y,z, a,b;
//int new_ncount=ncount;
	for (i=1;i<=netlist[0][1];i++)
	{
		if (net[i][0]>2)
		{
		wt=1/(float(net[i][0])-1);
		weight[i]=wt;

		}
	}

    // allocate space for Matrix data structure
    A_matrix = new float * [netlist[0][0]+1];
    C_matrix = new float * [netlist[0][0]+1];
    BX_matrix = new float[netlist[0][0]+1];
    BY_matrix = new float[netlist[0][0]+1];
	
    for (i=0; i<=netlist[0][0]; i++)
    {
	A_matrix[i] = new float[netlist[0][0]+1];
	C_matrix[i] = new float[netlist[0][0]+1];
	}

    for (i=0; i<=netlist[0][0]; i++)
    {
		BX_matrix[i] = 0;
		BY_matrix[i] = 0;
		for (j=0; j<=netlist[0][0]; j++)
		  {
			A_matrix[i][j] = 0;
			C_matrix[i][j] = 0;
		  }
	}


int k, m, flag;
flag=0;
	//Forming the C_matrix
    for (i=1; i<=fullnetlist[0][0]; i++)
    {
		for (j=1; j<=fullnetlist[0][0]; j++)
		  {
			flag=0;
			if(i!=j)
			{
				for (k=1; k<=fullnetlist[i][0]; k++)
					{
					for (m=1; m<=fullnetlist[j][0]; m++)
						{
						if ((fullnetlist[i][k]!=0)&&(fullnetlist[j][m]!=0)&&(fullnetlist[i][k]==fullnetlist[j][m]))
							{
							flag=1;
							C_matrix[i][j] = weight[fullnetlist[i][k]];
							}
						}
					}
			}
		  }
	}

//Forming the A_matrix and B_matrix

    for (i=1; i<=fullnetlist[0][0]; i++)
    {
		for (j=1; j<=fullnetlist[0][0]; j++)
		  {
			if(i!=j)
				{
				A_matrix[i][j]=-C_matrix[i][j];
				}
			if(i==j)
				{
	//			 A_matrix[i][j]=0; 
				for(k=1;k<=fullnetlist[0][0];k++)
					{
						A_matrix[i][j]=A_matrix[i][j]+C_matrix[i][k];
					}
				for(k=1;k<=fullnetlist[i][0];k++)
					{
					if(fullnetlist[i][k]!=0)
						{
						for(x=1;x<=net[fullnetlist[i][k]][0];x++)
							{
							if(net[fullnetlist[i][k]][x]>gcount)
								{
								  	A_matrix[i][j]=A_matrix[i][j]+weight[fullnetlist[i][k]];
									BX_matrix[i]=BX_matrix[i]+(weight[fullnetlist[i][k]]*(pads[net[fullnetlist[i][k]][x]-gcount][1]));
									BY_matrix[i]=BY_matrix[i]+(weight[fullnetlist[i][k]]*(pads[net[fullnetlist[i][k]][x]-gcount][2]));
								}
							}
						}
					}
				}
		  }
	}

    gate_xy = new float * [gcount+1];
	gate_xy[0] = new float[2];
    for (i=1;i<=gcount;i++)
	{
	gate_xy[i] = new float[2];
	}

    full_gate_xy = new float * [gcount+pads[0][0]+1];
	full_gate_xy[0] = new float[2];
    for (i=1;i<=gcount+pads[0][0];i++)
	{
	full_gate_xy[i] = new float[2];
	}

	matrixSolver (gate_xy, A_matrix, BX_matrix, gcount, 0 );
	matrixSolver (gate_xy, A_matrix, BY_matrix, gcount, 1 );

for(i=1;i<=gcount;i++)
{
full_gate_xy[i][0]=gate_xy[i][0];
full_gate_xy[i][1]=gate_xy[i][1];
}

for(i=1;i<=pads[0][0];i++)
{
full_gate_xy[gcount+i][0]=float(pads[i][1]);
full_gate_xy[gcount+i][1]=float(pads[i][2]);
}

//Print the random XY coordinates
   printf("\nOriginal XY coordinates:\n");
   printf("\n");
   for (i=1;i<=gcount;i++)
	{
	printf("%f %f ",gate_xy[i][0],gate_xy[i][1]);
	printf("\n");
	}


//Calculating original Quadratic Wirelength QL
float QL=0;
for (i=1;i<=fullnetlist[0][1];i++)
{
	QL = QL + weight[i]*(((full_gate_xy[net[i][1]][0]-full_gate_xy[net[i][2]][0])*(full_gate_xy[net[i][1]][0]-full_gate_xy[net[i][2]][0]))+((full_gate_xy[net[i][1]][1]-full_gate_xy[net[i][2]][1])*(full_gate_xy[net[i][1]][1]-full_gate_xy[net[i][2]][1])));
}
printf ("\n Original Quadratic Wirelength=%f\n",QL);

int *num_regions;
int *num_cuts;
    num_regions = new int[1];
    num_cuts = new int[1];
num_regions[0]=0;
num_cuts[0]=0;

	int *sort_gates;

//Assignment
    sort_gates = new int[gcount+1];
    for (i=1; i<=gcount; i++)
    {
	sort_gates[i]=i;
	}

float xmin=0;
float ymin=0;
float xmax=100;
float ymax=100;

verticalCut(gate_xy, full_gate_xy, A_matrix, fullnetlist, pads, weight, net, gcount, num_regions, num_cuts, sort_gates, xmin, ymin, xmax, ymax );

//Print the 3qp XY coordinates
   printf("\n3qp placement XY coordinates:\n\n");
    for (i=1; i<=gcount; i++)
    {
	fprintf(outFile, "%d ", i);
    fprintf(outFile, "%f %f ", gate_xy[i][0],gate_xy[i][1]);	
	fprintf(outFile, "\n");
	printf("%f %f ",gate_xy[i][0],gate_xy[i][1]);
	printf("\n");
	}
	fprintf(outFile, "\n");

    fclose(outFile);

for(i=1;i<=gcount;i++)
{
full_gate_xy[i][0]=gate_xy[i][0];
full_gate_xy[i][1]=gate_xy[i][1];
}

//Calculating 3QP Quadratic Wirelength QL
QL=0;
for (i=1;i<=fullnetlist[0][1];i++)
{
	QL = QL + weight[i]*(((full_gate_xy[net[i][1]][0]-full_gate_xy[net[i][2]][0])*(full_gate_xy[net[i][1]][0]-full_gate_xy[net[i][2]][0]))+((full_gate_xy[net[i][1]][1]-full_gate_xy[net[i][2]][1])*(full_gate_xy[net[i][1]][1]-full_gate_xy[net[i][2]][1])));
}
printf ("\n 3qp Quadratic Wirelength=%f\n",QL);


}

////////////////////////////////////////////////////////////////////////
// other member functions here
////////////////////////////////////////////////////////////////////////

int circuit::matrixSolver (float** gate_xy, float** A_matrix, float* BX_matrix, int gcount, int xy )
       {

	int i, j,k;
    float *gate_x0;
    float *gate_r0;
    float *gate_p0;
	float *Ap;

	float alpha=0;
	float beta=0;
	float rsold=0;
	float rsnew=0;
	float temp=0;

	gate_x0 = new float[gcount+1];
	gate_r0 = new float[gcount+1];
	gate_p0 = new float[gcount+1];
	Ap = new float[gcount+1];

    for (i=1;i<=gcount;i++)
	{
	gate_x0[i]=0;
	gate_r0[i]=0;
	gate_p0[i]=0;
	Ap[i] = 0;
	}

	//r=b-A*x
	for(i=1;i<=gcount;i++)
	{
	for (j=1;j<=gcount;j++)
		{
    gate_r0[i]=gate_r0[i] + A_matrix[i][j]*gate_x0[i];	
		}
	gate_r0[i]=BX_matrix[i]-gate_r0[i];
	}

	//p=r
	for(i=1;i<=gcount;i++)
	{
	gate_p0[i]=gate_r0[i];
	}

	//rsold=r'*r
	for(i=1;i<=gcount;i++)
	{
	rsold=rsold + gate_r0[i]*gate_r0[i];
	}

for(i=1;i<=1000000;i++)
	{
 		rsnew=0;
		temp=0;
	     //   Ap=A*p;
		for(j=1;j<=gcount;j++)
			{
			Ap[j]=0;
			for(k=1;k<=gcount;k++)
				{
				Ap[j]=Ap[j] + A_matrix[j][k]*gate_p0[k];
				}
			}

        for(j=1;j<=gcount;j++)
			{
			temp=temp+gate_p0[j]*Ap[j];
			}

//if(temp==0)
//break;
		//alpha=rsold/(p'*Ap)
		alpha=rsold/temp;

		//x=x+alpha*p
        for(j=1;j<=gcount;j++)
			{
			gate_x0[j]=gate_x0[j] + alpha*gate_p0[j];
			}

        
        //r=r-alpha*Ap;

		for(j=1;j<=gcount;j++)
			{
				gate_r0[j]=gate_r0[j] - alpha*Ap[j];
			}

        //rsnew=r'*r;
		
		for(j=1;j<=gcount;j++)
			{
			rsnew= rsnew + gate_r0[j]*gate_r0[j];
			}

        if (sqrt(rsnew)< exp(-10))
              break;

        //p=r+rsnew/rsold*p;

	for (j=1;j<=gcount;j++)
		{
    	gate_p0[j]=gate_r0[j] + (rsnew/rsold)*gate_p0[j];	
		}

        rsold=rsnew;
	}

	if (xy==0)
  	{
    	for (i=1; i<=gcount; i++)
    	{
		gate_xy[i][0]=gate_x0[i];
		}
  	}
	else
	{
    	for (i=1; i<=gcount; i++)
    	{
		gate_xy[i][1]=gate_x0[i];
		}

	}

}


int circuit::verticalCut (float** gate_xy, float** full_gate_xy, float** A_matrix, int** fullnetlist, int** pads, float* weight, int** net, int gcount, int* num_regions, int* num_cuts, int* sort_gates, float xmin, float ymin, float xmax, float ymax )

	{
	int i,j,k,x,y,z,u,v,w,l,m;
	int *loc_gates;
	float **A_matrix_left;
	float *BX_matrix_left;
	float *BY_matrix_left;
	float **A_matrix_right;
	float *BX_matrix_right;
	float *BY_matrix_right;
    float **temp_gate_xy;
    float xmin_l, ymin_l, xmax_l, ymax_l;
    float xmin_r, ymin_r, xmax_r, ymax_r;

	xmin_l=xmin;
	ymin_l=ymin;
	xmax_l=xmin+((xmax-xmin)/2);
	ymax_l=ymax;

	xmin_r=xmin+((xmax-xmin)/2);
	ymin_r=ymin;
	xmax_r=xmax;
	ymax_r=ymax;

printf("\n xmin_l=%f ymin_l=%f xmax_l=%f ymax_l=%f\n",xmin_l, ymin_l,xmax_l, ymax_l);
printf("\n xmin_r=%f ymin_r=%f xmax_r=%f ymax_r=%f\n",xmin_r, ymin_r,xmax_r, ymax_r);

//Assignment

    temp_gate_xy = new float * [fullnetlist[0][0]+1];
	temp_gate_xy[0] = new float[2];
    loc_gates = new int[fullnetlist[0][0]+1];

   for (i=1;i<=fullnetlist[0][0];i++)
	{
	temp_gate_xy[i] = new float[2];
	temp_gate_xy[i][0] = gate_xy[i][0];
	temp_gate_xy[i][1] = gate_xy[i][1];
	loc_gates[i]=1;
	}

printf("\nSorting gates....\n");

int temp;
    for (i=1;i<gcount;i++)
		{
			for (j=i+1;j<=gcount;j++)
				{
				if(gate_xy[sort_gates[i]][0]>gate_xy[sort_gates[j]][0])
					{
					temp=sort_gates[i];
					sort_gates[i]=sort_gates[j];
					sort_gates[j]=temp;
					}
				if(gate_xy[sort_gates[i]][0]==gate_xy[sort_gates[j]][0])
					{
					if(gate_xy[sort_gates[i]][1]>gate_xy[sort_gates[j]][1])
						{
							temp=sort_gates[i];
							sort_gates[i]=sort_gates[j];
							sort_gates[j]=temp;
						}
					}
				}
		}

float chk ;
int lcount;
int rcount;
chk = gcount%2;
if(chk!=0)
  {
	lcount=(gcount-1)/2;
	rcount=gcount-lcount;
  }
 else
	{
	lcount=(gcount)/2;
	rcount=gcount-lcount;
	}
    printf("\nLeft partion count=%d Right Partition count=%d\n",lcount,rcount);	

//Left-Side containment step

//Location of gates

    for (i=1;i<=lcount;i++)
	{
    loc_gates[sort_gates[i]]=0; 
	}

//Form A matrix for left side

    // allocate space for Left Matrix data structure
    A_matrix_left = new float * [lcount+1];
    BX_matrix_left = new float[lcount+1];
    BY_matrix_left = new float[lcount+1];
	
    for (i=0; i<=lcount; i++)
    {
	A_matrix_left[i] = new float[lcount+1];
	}

    for (i=1; i<=lcount; i++)
    {
		BX_matrix_left[i] = 0;
		BY_matrix_left[i] = 0;
		for (j=1; j<=lcount; j++)
		  {
				if(i!=j)
					{
        			A_matrix_left[i][j] = A_matrix[sort_gates[i]][sort_gates[j]];
					}
				if(i==j)
					{
				 	A_matrix_left[i][j] = 0;   	
					}
		  }
	}

float xcoor,ycoor, xcoor_new, ycoor_new;

    for (i=1; i<=lcount; i++)
    {
		for (j=1;j<=lcount;j++)
			{
			if(i!=j)
				{
				A_matrix_left[i][i] = A_matrix_left[i][i] - A_matrix_left[i][j];
				}
			}

		v=sort_gates[i];
		for(k=1;k<=fullnetlist[v][0];k++)
					{
					if(fullnetlist[v][k]!=0)
						{
						for(x=1;x<=net[fullnetlist[v][k]][0];x++)
							{
							if(net[fullnetlist[v][k]][x]<=fullnetlist[0][0])
								{
									if(loc_gates[net[fullnetlist[v][k]][x]]==1)
									{
									xcoor=gate_xy[net[fullnetlist[v][k]][x]][0];
									ycoor=gate_xy[net[fullnetlist[v][k]][x]][1];

									if((xcoor>xmax_l)&&(ycoor>ymax_l))
										{
										xcoor_new=xmax_l;
										ycoor_new=ymax_l;
										}
									if((xcoor>xmax_l)&&(ycoor<ymin_l))
										{
										xcoor_new=xmax_l;
										ycoor_new=ymin_l;
										}
									if((xcoor<xmin_l)&&(ycoor>ymax_l))
										{
										xcoor_new=xmin_l;
										ycoor_new=ymax_l;
										}
									if((xcoor<xmin_l)&&(ycoor<ymin_l))
										{
										xcoor_new=xmin_l;
										ycoor_new=ymin_l;
										}
									if((xcoor>=xmax_l)&&(ycoor<=ymax_l)&&(ycoor>=ymin_l))
										{
										xcoor_new=xmax_l;
										ycoor_new=ycoor;
										}
									if((xcoor<=xmin_l)&&(ycoor<=ymax_l)&&(ycoor>=ymin_l))
										{
										xcoor_new=xmin_l;
										ycoor_new=ycoor;
										}
									if((xcoor>=xmin_l)&&(xcoor<=xmax_l)&&(ycoor<=ymin_l))
										{
										xcoor_new=xcoor;
										ycoor_new=ymin_l;
										}
									if((xcoor>=xmin_l)&&(xcoor<=xmax_l)&&(ycoor>=ymax_l))
										{
										xcoor_new=xcoor;
										ycoor_new=ymax_l;
										}
									if((xcoor>xmin_l)&&(xcoor<xmax_l)&&(ycoor>ymin_l)&&(ycoor<ymax_l))
										{
										xcoor_new=xmax_l;
										ycoor_new=ycoor;
										}
								  	A_matrix_left[i][i]=A_matrix_left[i][i]+weight[fullnetlist[v][k]];
									BX_matrix_left[i]=BX_matrix_left[i]+(weight[fullnetlist[v][k]]*xcoor_new);
									BY_matrix_left[i]=BY_matrix_left[i]+(weight[fullnetlist[v][k]]*ycoor_new);
									}
								}
							if(net[fullnetlist[v][k]][x]>fullnetlist[0][0])
								{
									xcoor=float(pads[net[fullnetlist[v][k]][x]-fullnetlist[0][0]][1]);
									ycoor=float(pads[net[fullnetlist[v][k]][x]-fullnetlist[0][0]][2]);

									if((xcoor>xmax_l)&&(ycoor>ymax_l))
										{
										xcoor_new=xmax_l;
										ycoor_new=ymax_l;
										}
									if((xcoor>xmax_l)&&(ycoor<ymin_l))
										{
										xcoor_new=xmax_l;
										ycoor_new=ymin_l;
										}
									if((xcoor<xmin_l)&&(ycoor>ymax_l))
										{
										xcoor_new=xmin_l;
										ycoor_new=ymax_l;
										}
									if((xcoor<xmin_l)&&(ycoor<ymin_l))
										{
										xcoor_new=xmin_l;
										ycoor_new=ymin_l;
										}
									if((xcoor>=xmax_l)&&(ycoor<=ymax_l)&&(ycoor>=ymin_l))
										{
										xcoor_new=xmax_l;
										ycoor_new=ycoor;
										}
									if((xcoor<=xmin_l)&&(ycoor<=ymax_l)&&(ycoor>=ymin_l))
										{
										xcoor_new=xmin_l;
										ycoor_new=ycoor;
										}
									if((xcoor>=xmin_l)&&(xcoor<=xmax_l)&&(ycoor<=ymin_l))
										{
										xcoor_new=xcoor;
										ycoor_new=ymin_l;
										}
									if((xcoor>=xmin_l)&&(xcoor<=xmax_l)&&(ycoor>=ymax_l))
										{
										xcoor_new=xcoor;
										ycoor_new=ymax_l;
										}
									if((xcoor>xmin_l)&&(xcoor<xmax_l)&&(ycoor>ymin_l)&&(ycoor<ymax_l))
										{
										xcoor_new=xmax_l;
										ycoor_new=ycoor;
										}

									BX_matrix_left[i]=BX_matrix_left[i]+(weight[fullnetlist[v][k]]*(xcoor_new));
									BY_matrix_left[i]=BY_matrix_left[i]+(weight[fullnetlist[v][k]]*(ycoor_new));
									A_matrix_left[i][i]=A_matrix_left[i][i]+weight[fullnetlist[v][k]];
								}
							
							}
						}
					}
	}
//Form A matrix for right side

    for (i=1;i<=lcount;i++)
	{
    loc_gates[sort_gates[i]]=1; 
	}
    for (i=1;i<=rcount;i++)
	{
    loc_gates[sort_gates[lcount+i]]=0; 
	}

    // allocate space for Right Matrix data structure
    A_matrix_right = new float * [rcount+1];
    BX_matrix_right = new float[rcount+1];
    BY_matrix_right = new float[rcount+1];
	
    for (i=0; i<=rcount; i++)
    {
	A_matrix_right[i] = new float[rcount+1];
	}

    for (i=1; i<=rcount; i++)
    {
		BX_matrix_right[i] = 0;
		BY_matrix_right[i] = 0;
		for (j=1; j<=rcount; j++)
		  {
				if(i!=j)
					{
        			A_matrix_right[i][j] = A_matrix[sort_gates[lcount+i]][sort_gates[lcount+j]];
					}
				if(i==j)
					{
				 	A_matrix_right[i][j] = 0;   	
					}

		  }
	}

    for (i=1; i<=rcount; i++)
    {
		for (j=1;j<=rcount;j++)
			{
				if(i!=j)
				{
				A_matrix_right[i][i] = A_matrix_right[i][i] - A_matrix_right[i][j];
				}
			}

		v=sort_gates[lcount+i];
		for(k=1;k<=fullnetlist[v][0];k++)
					{
					if(fullnetlist[v][k]!=0)
						{
						for(x=1;x<=net[fullnetlist[v][k]][0];x++)
							{
							if(net[fullnetlist[v][k]][x]<=fullnetlist[0][0])
								{
									if(loc_gates[net[fullnetlist[v][k]][x]]==1)
									{
									xcoor=gate_xy[net[fullnetlist[v][k]][x]][0];
									ycoor=gate_xy[net[fullnetlist[v][k]][x]][1];

									if((xcoor>xmax_r)&&(ycoor>ymax_r))
										{
										xcoor_new=xmax_r;
										ycoor_new=ymax_r;
										}
									if((xcoor>xmax_r)&&(ycoor<ymin_r))
										{
										xcoor_new=xmax_r;
										ycoor_new=ymin_r;
										}
									if((xcoor<xmin_r)&&(ycoor>ymax_r))
										{
										xcoor_new=xmin_r;
										ycoor_new=ymax_r;
										}
									if((xcoor<xmin_r)&&(ycoor<ymin_r))
										{
										xcoor_new=xmin_r;
										ycoor_new=ymin_r;
										}
									if((xcoor>=xmax_r)&&(ycoor<=ymax_r)&&(ycoor>=ymin_r))
										{
										xcoor_new=xmax_r;
										ycoor_new=ycoor;
										}
									if((xcoor<=xmin_r)&&(ycoor<=ymax_r)&&(ycoor>=ymin_r))
										{
										xcoor_new=xmin_r;
										ycoor_new=ycoor;
										}
									if((xcoor>=xmin_r)&&(xcoor<=xmax_r)&&(ycoor<=ymin_r))
										{
										xcoor_new=xcoor;
										ycoor_new=ymin_r;
										}
									if((xcoor>=xmin_r)&&(xcoor<=xmax_r)&&(ycoor>=ymax_r))
										{
										xcoor_new=xcoor;
										ycoor_new=ymax_r;
										}
									if((xcoor>xmin_r)&&(xcoor<xmax_r)&&(ycoor>ymin_r)&&(ycoor<ymax_r))
										{
										xcoor_new=xmin_r;
										ycoor_new=ycoor;
										}
								  	A_matrix_right[i][i]=A_matrix_right[i][i]+weight[fullnetlist[v][k]];
									BX_matrix_right[i]=BX_matrix_right[i]+(weight[fullnetlist[v][k]]*(xcoor_new));
									BY_matrix_right[i]=BY_matrix_right[i]+(weight[fullnetlist[v][k]]*(ycoor_new));
									}
								}

							if(net[fullnetlist[v][k]][x]>fullnetlist[0][0])
								{
									xcoor=float(pads[net[fullnetlist[v][k]][x]-fullnetlist[0][0]][1]);
									ycoor=float(pads[net[fullnetlist[v][k]][x]-fullnetlist[0][0]][2]);

									if((xcoor>xmax_r)&&(ycoor>ymax_r))
										{
										xcoor_new=xmax_r;
										ycoor_new=ymax_r;
										}
									if((xcoor>xmax_r)&&(ycoor<ymin_r))
										{
										xcoor_new=xmax_r;
										ycoor_new=ymin_r;
										}
									if((xcoor<xmin_r)&&(ycoor>ymax_r))
										{
										xcoor_new=xmin_r;
										ycoor_new=ymax_r;
										}
									if((xcoor<xmin_r)&&(ycoor<ymin_r))
										{
										xcoor_new=xmin_r;
										ycoor_new=ymin_r;
										}
									if((xcoor>=xmax_r)&&(ycoor<=ymax_r)&&(ycoor>=ymin_r))
										{
										xcoor_new=xmax_r;
										ycoor_new=ycoor;
										}
									if((xcoor<=xmin_r)&&(ycoor<=ymax_r)&&(ycoor>=ymin_r))
										{
										xcoor_new=xmin_r;
										ycoor_new=ycoor;
										}
									if((xcoor>=xmin_r)&&(xcoor<=xmax_r)&&(ycoor<=ymin_r))
										{
										xcoor_new=xcoor;
										ycoor_new=ymin_r;
										}
									if((xcoor>=xmin_r)&&(xcoor<=xmax_r)&&(ycoor>=ymax_r))
										{
										xcoor_new=xcoor;
										ycoor_new=ymax_r;
										}
									if((xcoor>xmin_r)&&(xcoor<xmax_r)&&(ycoor>ymin_r)&&(ycoor<ymax_r))
										{
										xcoor_new=xmax_r;
										ycoor_new=ycoor;
										}

									BX_matrix_right[i]=BX_matrix_right[i]+(weight[fullnetlist[v][k]]*(xcoor_new));
									BY_matrix_right[i]=BY_matrix_right[i]+(weight[fullnetlist[v][k]]*(ycoor_new));
									A_matrix_right[i][i]=A_matrix_right[i][i]+weight[fullnetlist[v][k]];
								}
							}
						}
					}
	}


	matrixSolver (temp_gate_xy, A_matrix_left, BX_matrix_left, lcount, 0 );
    for (i=1;i<=lcount;i++)
	{
	gate_xy[sort_gates[i]][0] =  temp_gate_xy[i][0];
	}

	matrixSolver (temp_gate_xy, A_matrix_left, BY_matrix_left, lcount, 1 );

    for (i=1;i<=lcount;i++)
	{
	gate_xy[sort_gates[i]][1] =  temp_gate_xy[i][1];
	}

	matrixSolver (temp_gate_xy, A_matrix_right, BX_matrix_right, rcount, 0 );
   for (i=1;i<=rcount;i++)
	{
	gate_xy[sort_gates[lcount+i]][0] = temp_gate_xy[i][0] ;
	}

	matrixSolver (temp_gate_xy, A_matrix_right, BY_matrix_right, rcount, 1 );

    for (i=1;i<=rcount;i++)
	{
	gate_xy[sort_gates[lcount+i]][1] = temp_gate_xy[i][1] ;
	}

	int *sort_l_gates;

    sort_l_gates = new int[lcount+1];
    for (i=1; i<=lcount; i++)
    {
	sort_l_gates[i]=sort_gates[i];
	}

	int *sort_r_gates;

    sort_r_gates = new int[rcount+1];
    for (i=1; i<=rcount; i++)
    {
	sort_r_gates[i]=sort_gates[lcount+i];
	}
num_regions[0]++;
printf("\n No. of region in ver cut=%d \n",num_regions[0]);
horizontalCut (gate_xy, full_gate_xy, A_matrix, fullnetlist, pads, weight, net, lcount, num_regions, num_cuts, sort_l_gates, xmin_l, ymin_l, xmax_l, ymax_l );
horizontalCut (gate_xy, full_gate_xy, A_matrix, fullnetlist, pads, weight, net, rcount, num_regions, num_cuts, sort_r_gates, xmin_r, ymin_r, xmax_r, ymax_r );
num_regions[0]--;
}



int circuit::horizontalCut (float** gate_xy, float** full_gate_xy, float** A_matrix, int** fullnetlist, int** pads, float* weight, int** net, int gcount, int* num_regions,int* num_cuts, int* sort_gates, float xmin, float ymin, float xmax, float ymax )

	{
	int i,j,k,x,y,z,u,v,w,l,m;
	int *loc_gates;
	float **A_matrix_down;
	float *BX_matrix_down;
	float *BY_matrix_down;
	float **A_matrix_up;
	float *BX_matrix_up;
	float *BY_matrix_up;
    float **temp_gate_xy;
    float xmin_d, ymin_d, xmax_d, ymax_d;
    float xmin_u, ymin_u, xmax_u, ymax_u;
	float xcoor,ycoor, xcoor_new, ycoor_new;

	xmin_d=xmin;
	ymin_d=ymin;
	xmax_d=xmax;
	ymax_d=ymin+((ymax-ymin)/2);

	xmin_u=xmin;
	ymin_u=ymin+((ymax-ymin)/2);
	xmax_u=xmax;
	ymax_u=ymax;
printf("\n xmin_d=%f ymin_d=%f xmax_d=%f ymax_d=%f\n",xmin_d, ymin_d,xmax_d, ymax_d);
printf("\n xmin_u=%f ymin_u=%f xmax_u=%f ymax_u=%f\n",xmin_u, ymin_u,xmax_u, ymax_u);

//Assignment

    temp_gate_xy = new float * [fullnetlist[0][0]+1];
	temp_gate_xy[0] = new float[2];
    loc_gates = new int[fullnetlist[0][0]+1];

   for (i=1;i<=fullnetlist[0][0];i++)
	{
	temp_gate_xy[i] = new float[2];
	temp_gate_xy[i][0] = gate_xy[i][0];
	temp_gate_xy[i][1] = gate_xy[i][1];
	loc_gates[i]=1;
	}


printf("\nSorting gates....\n");

int temp;
    for (i=1;i<gcount;i++)
		{
			for (j=i+1;j<=gcount;j++)
				{
				if(gate_xy[sort_gates[i]][1]>gate_xy[sort_gates[j]][1])
					{
					temp=sort_gates[i];
					sort_gates[i]=sort_gates[j];
					sort_gates[j]=temp;
					}
				if(gate_xy[sort_gates[i]][1]==gate_xy[sort_gates[j]][1])
					{
					if(gate_xy[sort_gates[i]][0]>gate_xy[sort_gates[j]][0])
						{
							temp=sort_gates[i];
							sort_gates[i]=sort_gates[j];
							sort_gates[j]=temp;
						}
					}
				}
		}

float chk ;
int dcount;
int ucount;
chk = gcount%2;
if(chk!=0)
  {
	dcount=(gcount-1)/2;
	ucount=gcount-dcount;
  }
 else
	{
	dcount=(gcount)/2;
	ucount=gcount-dcount;
	}
    printf("\nDown partion count=%d Up Partition count=%d\n",dcount,ucount);	

//Left-Side containment step

//Location of gates
  //  loc_gates = new int[fullnetlist[0][0]+1];
    for (i=1;i<=dcount;i++)
	{
    loc_gates[sort_gates[i]]=0; 
	}
//Form A matrix for left side

    // allocate space for Left Matrix data structure
    A_matrix_down = new float * [dcount+1];
    BX_matrix_down = new float[dcount+1];
    BY_matrix_down = new float[dcount+1];
	
    for (i=0; i<=dcount; i++)
    {
	A_matrix_down[i] = new float[dcount+1];
	}

    for (i=1; i<=dcount; i++)
    {
		BX_matrix_down[i] = 0;
		BY_matrix_down[i] = 0;
		for (j=1; j<=dcount; j++)
		  {
				if(i!=j)
					{
        			A_matrix_down[i][j] = A_matrix[sort_gates[i]][sort_gates[j]];
					}
				if(i==j)
					{
				 	A_matrix_down[i][j] = 0;   	
					}
		  }
	}

    for (i=1; i<=dcount; i++)
    {
		for (j=1;j<=dcount;j++)
			{
			if(i!=j)
				{
				A_matrix_down[i][i] = A_matrix_down[i][i] - A_matrix_down[i][j];
				}
			}

		v=sort_gates[i];
		for(k=1;k<=fullnetlist[v][0];k++)
					{
					if(fullnetlist[v][k]!=0)
						{
						for(x=1;x<=net[fullnetlist[v][k]][0];x++)
							{
							if(net[fullnetlist[v][k]][x]<=fullnetlist[0][0])
								{
									if(loc_gates[net[fullnetlist[v][k]][x]]==1)
									{
									xcoor=gate_xy[net[fullnetlist[v][k]][x]][0];
									ycoor=gate_xy[net[fullnetlist[v][k]][x]][1];

									if((xcoor>xmax_d)&&(ycoor>ymax_d))
										{
										xcoor_new=xmax_d;
										ycoor_new=ymax_d;
										}
									if((xcoor>xmax_d)&&(ycoor<ymin_d))
										{
										xcoor_new=xmax_d;
										ycoor_new=ymin_d;
										}
									if((xcoor<xmin_d)&&(ycoor>ymax_d))
										{
										xcoor_new=xmin_d;
										ycoor_new=ymax_d;
										}
									if((xcoor<xmin_d)&&(ycoor<ymin_d))
										{
										xcoor_new=xmin_d;
										ycoor_new=ymin_d;
										}
									if((xcoor>=xmax_d)&&(ycoor<=ymax_d)&&(ycoor>=ymin_d))
										{
										xcoor_new=xmax_d;
										ycoor_new=ycoor;
										}
									if((xcoor<=xmin_d)&&(ycoor<=ymax_d)&&(ycoor>=ymin_d))
										{
										xcoor_new=xmin_d;
										ycoor_new=ycoor;
										}
									if((xcoor>=xmin_d)&&(xcoor<=xmax_d)&&(ycoor<=ymin_d))
										{
										xcoor_new=xcoor;
										ycoor_new=ymin_d;
										}
									if((xcoor>=xmin_d)&&(xcoor<=xmax_d)&&(ycoor>=ymax_d))
										{
										xcoor_new=xcoor;
										ycoor_new=ymax_d;
										}
									if((xcoor>xmin_d)&&(xcoor<xmax_d)&&(ycoor>ymin_d)&&(ycoor<ymax_d))
										{
										xcoor_new=xcoor;
										ycoor_new=ymax_d;
										}
								  	A_matrix_down[i][i]=A_matrix_down[i][i]+weight[fullnetlist[v][k]];
									BY_matrix_down[i]=BY_matrix_down[i]+(weight[fullnetlist[v][k]]*ycoor_new);
									BX_matrix_down[i]=BX_matrix_down[i]+(weight[fullnetlist[v][k]]*xcoor_new);
									}
								}
							if(net[fullnetlist[v][k]][x]>fullnetlist[0][0])
								{
									xcoor=float(pads[net[fullnetlist[v][k]][x]-fullnetlist[0][0]][1]);
									ycoor=float(pads[net[fullnetlist[v][k]][x]-fullnetlist[0][0]][2]);
									if((xcoor>xmax_d)&&(ycoor>ymax_d))
										{
										xcoor_new=xmax_d;
										ycoor_new=ymax_d;
										}
									if((xcoor>xmax_d)&&(ycoor<ymin_d))
										{
										xcoor_new=xmax_d;
										ycoor_new=ymin_d;
										}
									if((xcoor<xmin_d)&&(ycoor>ymax_d))
										{
										xcoor_new=xmin_d;
										ycoor_new=ymax_d;
										}
									if((xcoor<xmin_d)&&(ycoor<ymin_d))
										{
										xcoor_new=xmin_d;
										ycoor_new=ymin_d;
										}
									if((xcoor>=xmax_d)&&(ycoor<=ymax_d)&&(ycoor>=ymin_d))
										{
										xcoor_new=xmax_d;
										ycoor_new=ycoor;
										}
									if((xcoor<=xmin_d)&&(ycoor<=ymax_d)&&(ycoor>=ymin_d))
										{
										xcoor_new=xmin_d;
										ycoor_new=ycoor;
										}
									if((xcoor>=xmin_d)&&(xcoor<=xmax_d)&&(ycoor<=ymin_d))
										{
										xcoor_new=xcoor;
										ycoor_new=ymin_d;
										}
									if((xcoor>=xmin_d)&&(xcoor<=xmax_d)&&(ycoor>=ymax_d))
										{
										xcoor_new=xcoor;
										ycoor_new=ymax_d;
										}
									if((xcoor>xmin_d)&&(xcoor<xmax_d)&&(ycoor>ymin_d)&&(ycoor<ymax_d))
										{
										xcoor_new=xcoor;
										ycoor_new=ymax_d;
										}

									BY_matrix_down[i]=BY_matrix_down[i]+(weight[fullnetlist[v][k]]*ycoor_new);
									BX_matrix_down[i]=BX_matrix_down[i]+(weight[fullnetlist[v][k]]*xcoor_new);
									A_matrix_down[i][i]=A_matrix_down[i][i]+weight[fullnetlist[v][k]];
								}
							
							}
						}
					}
	}
//Form A matrix for right side

    for (i=1;i<=dcount;i++)
	{
    loc_gates[sort_gates[i]]=1; 
	}

    for (i=1;i<=ucount;i++)
	{
    loc_gates[sort_gates[dcount+i]]=0; 
	}

    // allocate space for Right Matrix data structure
    A_matrix_up = new float * [ucount+1];
    BX_matrix_up = new float[ucount+1];
    BY_matrix_up = new float[ucount+1];
	
    for (i=0; i<=ucount; i++)
    {
	A_matrix_up[i] = new float[ucount+1];
	}

    for (i=1; i<=ucount; i++)
    {
		BX_matrix_up[i] = 0;
		BY_matrix_up[i] = 0;
		for (j=1; j<=ucount; j++)
		  {
				if(i!=j)
					{
        			A_matrix_up[i][j] = A_matrix[sort_gates[dcount+i]][sort_gates[dcount+j]];
					}
				if(i==j)
					{
				 	A_matrix_up[i][j] = 0;   	
					}

		  }
	}

    for (i=1; i<=ucount; i++)
    {
		for (j=1;j<=ucount;j++)
			{
				if(i!=j)
				{
				A_matrix_up[i][i] = A_matrix_up[i][i] - A_matrix_up[i][j];
				}
			}

		v=sort_gates[dcount+i];
		for(k=1;k<=fullnetlist[v][0];k++)
					{
					if(fullnetlist[v][k]!=0)
						{
						for(x=1;x<=net[fullnetlist[v][k]][0];x++)
							{
							if(net[fullnetlist[v][k]][x]<=fullnetlist[0][0])
								{
									if(loc_gates[net[fullnetlist[v][k]][x]]==1)
									{
									xcoor=gate_xy[net[fullnetlist[v][k]][x]][0];
									ycoor=gate_xy[net[fullnetlist[v][k]][x]][1];

									if((xcoor>xmax_u)&&(ycoor>ymax_u))
										{
										xcoor_new=xmax_u;
										ycoor_new=ymax_u;
										}
									if((xcoor>xmax_u)&&(ycoor<ymin_u))
										{
										xcoor_new=xmax_u;
										ycoor_new=ymin_u;
										}
									if((xcoor<xmin_u)&&(ycoor>ymax_u))
										{
										xcoor_new=xmin_u;
										ycoor_new=ymax_u;
										}
									if((xcoor<xmin_u)&&(ycoor<ymin_u))
										{
										xcoor_new=xmin_u;
										ycoor_new=ymin_u;
										}
									if((xcoor>=xmax_u)&&(ycoor<=ymax_u)&&(ycoor>=ymin_u))
										{
										xcoor_new=xmax_u;
										ycoor_new=ycoor;
										}
									if((xcoor<=xmin_u)&&(ycoor<=ymax_u)&&(ycoor>=ymin_u))
										{
										xcoor_new=xmin_u;
										ycoor_new=ycoor;
										}
									if((xcoor>=xmin_u)&&(xcoor<=xmax_u)&&(ycoor<=ymin_u))
										{
										xcoor_new=xcoor;
										ycoor_new=ymin_u;
										}
									if((xcoor>=xmin_u)&&(xcoor<=xmax_u)&&(ycoor>=ymax_u))
										{
										xcoor_new=xcoor;
										ycoor_new=ymax_u;
										}
									if((xcoor>xmin_u)&&(xcoor<xmax_u)&&(ycoor>ymin_u)&&(ycoor<ymax_u))
										{
										xcoor_new=xcoor;
										ycoor_new=ymin_u;
										}

								  	A_matrix_up[i][i]=A_matrix_up[i][i]+weight[fullnetlist[v][k]];
									BY_matrix_up[i]=BY_matrix_up[i]+(weight[fullnetlist[v][k]]*ycoor_new);
									BX_matrix_up[i]=BX_matrix_up[i]+(weight[fullnetlist[v][k]]*xcoor_new);
									}
								}

							if(net[fullnetlist[v][k]][x]>fullnetlist[0][0])
								{
									xcoor=float(pads[net[fullnetlist[v][k]][x]-fullnetlist[0][0]][1]);
									ycoor=float(pads[net[fullnetlist[v][k]][x]-fullnetlist[0][0]][2]);

									if((xcoor>xmax_u)&&(ycoor>ymax_u))
										{
										xcoor_new=xmax_u;
										ycoor_new=ymax_u;
										}
									if((xcoor>xmax_u)&&(ycoor<ymin_u))
										{
										xcoor_new=xmax_u;
										ycoor_new=ymin_u;
										}
									if((xcoor<xmin_u)&&(ycoor>ymax_u))
										{
										xcoor_new=xmin_u;
										ycoor_new=ymax_u;
										}
									if((xcoor<xmin_u)&&(ycoor<ymin_u))
										{
										xcoor_new=xmin_u;
										ycoor_new=ymin_u;
										}
									if((xcoor>=xmax_u)&&(ycoor<=ymax_u)&&(ycoor>=ymin_u))
										{
										xcoor_new=xmax_u;
										ycoor_new=ycoor;
										}
									if((xcoor<=xmin_u)&&(ycoor<=ymax_u)&&(ycoor>=ymin_u))
										{
										xcoor_new=xmin_u;
										ycoor_new=ycoor;
										}
									if((xcoor>=xmin_u)&&(xcoor<=xmax_u)&&(ycoor<=ymin_u))
										{
										xcoor_new=xcoor;
										ycoor_new=ymin_u;
										}
									if((xcoor>=xmin_u)&&(xcoor<=xmax_u)&&(ycoor>=ymax_u))
										{
										xcoor_new=xcoor;
										ycoor_new=ymax_u;
										}
									if((xcoor>xmin_u)&&(xcoor<xmax_u)&&(ycoor>ymin_u)&&(ycoor<ymax_u))
										{
										xcoor_new=xcoor;
										ycoor_new=ymin_u;
										}

									BY_matrix_up[i]=BY_matrix_up[i]+(weight[fullnetlist[v][k]]*ycoor_new);
									BX_matrix_up[i]=BX_matrix_up[i]+(weight[fullnetlist[v][k]]*xcoor_new);
									A_matrix_up[i][i]=A_matrix_up[i][i]+weight[fullnetlist[v][k]];
								}

							}
						}
					}
	}

	matrixSolver (temp_gate_xy, A_matrix_down, BX_matrix_down, dcount, 0 );
    for (i=1;i<=dcount;i++)
	{
	gate_xy[sort_gates[i]][0] =  temp_gate_xy[i][0];
	}

	matrixSolver (temp_gate_xy, A_matrix_down, BY_matrix_down, dcount, 1 );

    for (i=1;i<=dcount;i++)
	{
	gate_xy[sort_gates[i]][1] =  temp_gate_xy[i][1];
	}

	matrixSolver (temp_gate_xy, A_matrix_up, BX_matrix_up, ucount, 0 );
   for (i=1;i<=ucount;i++)
	{
	gate_xy[sort_gates[dcount+i]][0] = temp_gate_xy[i][0] ;
	}

	matrixSolver (temp_gate_xy, A_matrix_up, BY_matrix_up, ucount, 1 );

    for (i=1;i<=ucount;i++)
	{
	gate_xy[sort_gates[dcount+i]][1] = temp_gate_xy[i][1] ;
	}

	int *sort_u_gates;

    sort_u_gates = new int[ucount+1];
    for (i=1; i<=ucount; i++)
    {
	sort_u_gates[i]=sort_gates[i];
	}

	int *sort_d_gates;

    sort_d_gates = new int[dcount+1];
    for (i=1; i<=dcount; i++)
    {
	sort_d_gates[i]=sort_gates[dcount+i];
	}
num_regions[0]++;

printf("\n No. of regions in hor cut =%d \n",num_regions[0]);


if(num_regions[0]==6)
  {
num_regions[0]=num_regions[0]-1;
return 0;
  }

verticalCut(gate_xy, full_gate_xy, A_matrix, fullnetlist, pads, weight, net, ucount, num_regions,num_cuts, sort_u_gates, xmin_d, ymin_d, xmax_d, ymax_d );
verticalCut(gate_xy, full_gate_xy, A_matrix, fullnetlist, pads, weight, net, dcount, num_regions,num_cuts, sort_d_gates, xmin_u, ymin_u, xmax_u, ymax_u );
num_regions[0]=num_regions[0]-1;
	}


//void strcpy (char *s, char *t)
//{
//    while (*s++ = *t++)
//	;
//}


//char * strcat ( char * destination, const char * source )
//{
//char *d = destination;
//while (*d) ++d;
//while ((*d++ = *source++) != '\0') ;
//return (destination);
//}

////////////////////////////////////////////////////////////////////////
// main starts here
////////////////////////////////////////////////////////////////////////
main(int argc, char *argv[])
{
    char cktName[20];
    circuit *ckt;
    int vecWidth, poNum;

    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s <ckt>\n", argv[0]);
        exit(-1);
    }

    strcpy(cktName, argv[1]);

    // read in the circuit and build data structure
    ckt = new circuit(cktName);

    // NOW ADD YOUR FUNCTION HERE
}
