/*
 * heat_mpi.cpp
 * make
 * mpirun -n 4 ./heat_mpi
 */

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "mpi.h"

# define n 200
# define nodeedge 200
# define nblock n/nodeedge
# define nproc nblock*nblock

int main ( int argc, char **argv );
void doblack ( double w, double M[][nodeedge+2] );
void dored ( double w, double M[][nodeedge+2] );
void exchange ( double M[][nodeedge+2], int comm[], int rank );
void iterate ( double w, double M[][nodeedge+2], double result[][n], int rank, int comm[] );
void setcomm ( int rank, int comm[] );
void setex ( double ex[], double M[][nodeedge+2], int which );
void initialize_matrix ( double M[][nodeedge+2] );
void unpack ( double M[][nodeedge+2], int where, double in[] );

int main ( int argc, char **argv )
{
  int comm[4];
  FILE *fp;
  int i;
  int j;
  double M[nodeedge+2][nodeedge+2];
  int ntasks;
  int rank;
  double result[n][n];
  double w;
  double wtime;

  MPI_Init ( &argc, &argv );

  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  MPI_Comm_size ( MPI_COMM_WORLD, &ntasks );

  wtime = MPI_Wtime ( );

  if ( rank == 0 ) 
  {
    printf ( "\n" );
    printf ( "LAPLACE_MPI:\n" );
    printf ( "  C/MPI version\n" );
    printf ( "  Solve the Laplace equation using MPI.\n" );
  }

  if ( ntasks != nproc )
  {
    if ( rank == 0 ) 
    {
      printf ( "\n" );
      printf ( "Fatal error!\n" );
      printf ( "  MP_PROCS should be set to %i!\n", nproc );
    }
    MPI_Finalize ( );
    exit ( 1 );
  }

  if ( rank == 0 ) 
  {
    printf ( "\n" );
    printf ( "  MPI has been set up.\n" );
  }

  if ( rank == 0 ) 
  {
    printf ( "  Initialize the matrix M.\n" );
  }
  initialize_matrix ( M );

  if ( rank == 0 ) 
  {
    printf ( "  Set the list of neighbors.\n" );
  }
  setcomm ( rank, comm );

  if ( rank == 0 ) 
  {
    printf ( "  Begin the iteration.\n" );
  }
  w = 1.2;
  iterate ( w, M, result, rank, comm );

  wtime = MPI_Wtime ( ) - wtime;

  printf ( "  Task %i took %6.3f seconds\n", rank, wtime );

  if ( rank == 0 )
  {
    fp = fopen ( "laplace_solution.txt", "w" );

    for ( i = 0; i < n; i++ ) 
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf ( fp, "%f \n", result[i][j] );
      }  
    }
    fclose ( fp );
    printf ( "  Solution written to \"laplace_solution.txt\".\n" );
  }

  MPI_Finalize ( );

  if ( rank == 0 )
  {
    printf ( "\n" );
    printf ( "LAPLACE_MPI:\n" );
    printf ( "  Normal end of execution.\n" );
  }
  return 0;
}

void doblack ( double w, double M[][nodeedge+2] )
{
  int i;
  int j;

  for ( i = 1; i <= nodeedge / 2; i++ )
  {
    for ( j = nodeedge / 2 + 1; j <= nodeedge; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }

  for ( i = nodeedge / 2 + 1; i <= nodeedge; i++ )
  {
    for ( j = 1; j <= nodeedge / 2; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
  return;
}

void dored ( double w, double M[][nodeedge+2] ) 
{
  int i;
  int j;

  for ( i = 1; i <= nodeedge / 2; i++ )
  {
    for ( j = 1; j <= nodeedge / 2; j++ ) 
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }

  for ( i = nodeedge / 2 + 1; i <= nodeedge; i++ )
  {
    for ( j = nodeedge / 2 + 1; j <= nodeedge; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
  return;
}

void exchange ( double M[][nodeedge+2], int comm[], int rank )
{
  double ex0[nodeedge];
  double ex1[nodeedge];
  double ex2[nodeedge];
  double ex3[nodeedge];
  int i;
  double in0[nodeedge];
  double in1[nodeedge];
  double in2[nodeedge];
  double in3[nodeedge];
  int partner;
  MPI_Request requests[8];
  MPI_Status status[8];
  int tag;

  for ( i = 0; i < 8; i++ ) 
  {
    requests[i] = MPI_REQUEST_NULL; 
  }

  if ( comm[0] == 1 )
  {
    partner = rank - nblock;
    tag = 0;
    MPI_Irecv ( &in0, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD, 
      &requests[0] );
  }

  if ( comm[1] == 1 )
  {
    partner = rank + 1;
    tag = 1;
    MPI_Irecv ( &in1, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[1] );
  }

  if ( comm[2] == 1 )
  {
    partner = rank + nblock;
    tag = 2;
    MPI_Irecv ( &in2, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[2] );
  }

  if ( comm[3] == 1 )
  {
    partner = rank - 1;
    tag = 3;
    MPI_Irecv ( &in3, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[3] );
  }

  if ( comm[0] == 1 )
  {
    partner = rank - nblock;
    tag = 2;
    setex ( ex0, M, 0 );
    MPI_Isend ( &ex0, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[4] );
  }

  if (comm[1] == 1 )
  {
    partner = rank + 1;
    tag = 3;
    setex ( ex1, M, 1 );
    MPI_Isend ( &ex1, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[5] );
  }

  if ( comm[2] == 1 )
  {
    partner = rank + nblock;
    tag = 0;
    setex ( ex2, M, 2 );
    MPI_Isend ( &ex2, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[6] );
  }

  if ( comm[3] == 1 )
  {
    partner = rank - 1;
    tag = 1;
    setex ( ex3, M, 3 );
    MPI_Isend ( &ex3, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[7] );
  }

  MPI_Waitall ( 8, requests, status );

  if ( comm[0] == 1 ) 
  {
    unpack ( M, 0, in0 );
  }
  if ( comm[1] == 1 ) 
  {
    unpack ( M, 1, in1 );
  }
  if ( comm[2] == 1 ) 
  {
    unpack ( M, 2, in2 );
  }
  if ( comm[3] == 1 ) 
  {
    unpack ( M, 3, in3 );
  }

  return;
}

void initialize_matrix ( double M[][nodeedge+2] )
{
  double avg;
  double bv[4];
  int i;
  int j;

  bv[0] = 100.0;
  bv[1] = 0.0;
  bv[2] = 0.0;
  bv[3] = 0.0;

  for ( i = 1; i <= nodeedge; i++ )
  { 
    M[0][i] =          bv[0];
    M[i][nodeedge+1] = bv[1];
    M[nodeedge+1][i] = bv[2];
    M[i][0] =          bv[3];
  }

  avg = ( bv[0] + bv[1] + bv[2] + bv[3] ) / 4.0;

  for ( i = 1; i <= nodeedge; i++ )
  {
    for ( j = 1; j <= nodeedge; j++ )
    {
      M[i][j] = avg;
    }
  }

  return;
}

void iterate ( double w, double M[][nodeedge+2], double result[][n], int rank, 
  int comm[] )
{
  int count;
  double diff;
  int done;
  double ediff;
  int i;
  double in;
  int index;
  int it;
  int j;
  int k;
  int l;
  double MM[n*n];
  double mold[nodeedge+2][nodeedge+2];
  double send[nodeedge][nodeedge];

  it = 0;
  done = 0;
  for ( i = 1; i <= nodeedge; i++ )
  {
    for ( j = 1; j <= nodeedge; j++ )
    {
      mold[i][j] = M[i][j];
    }
  }

  while ( done == 0 )
  {
    it++;
    exchange ( M, comm, rank );
    dored ( w, M );
    exchange ( M, comm, rank );
    doblack ( w, M );

    if ( 10000 < it )
    {
      done = 1;
    }

    if ( ( it % 20 == 0 ) && ( done != 1 ) )
    { 
      diff = 0.0;
      for ( i = 1; i <= nodeedge; i++ )
      {
        for ( j = 1; j <= nodeedge; j++ )
        {
          ediff = M[i][j] - mold[i][j];
          if ( ediff < 0.0 ) 
          {
            ediff = - ediff;
          }
          diff = diff + ediff;
          mold[i][j] = M[i][j];
        }
      }
      diff = diff / ( ( double ) ( nodeedge * nodeedge ) );

      MPI_Allreduce ( &diff, &in, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      if ( in < ( double ) nproc * 0.001 ) 
      {
        done = 1;
      }
    }
  }

  for ( i = 0; i < nodeedge; i++ )
  {
    for ( j = 0; j < nodeedge; j++ )
    {
      send[i][j] = M[i+1][j+1];
    }
  }

  count = nodeedge * nodeedge;

  MPI_Gather ( &send, count, MPI_DOUBLE, &MM, count, MPI_DOUBLE, 0, 
    MPI_COMM_WORLD );

  printf ( "  ITERATE gathered updated results to process 0.\n" );

  if ( rank == 0 ) 
  {
    printf ( "did %i iterations\n", it );

    index = 0;
    for ( k = 0; k < nblock; k++ )
    {
      for ( l = 0; l < nblock; l++ )
      {
        for ( i = k * nodeedge; i < ( k + 1 ) * nodeedge; i++ )
        {
          for ( j = l * nodeedge; j < ( l + 1 ) * nodeedge; j++ )
          {
            result[i][j] = MM[index];
            index++;
          }
        }
      }
    }
  }
  return;
}

void setcomm ( int rank, int comm[] )
{
  int i;

  for ( i = 0; i < 4; i++ ) 
  {
    comm[i] = 1;
  }

  if ( rank < nblock )
  {
    comm[0] = 0;    
  }

  if ( ( rank + 1 ) % nblock == 0 )
  {
    comm[1] = 0;
  }

  if ( rank > (nblock*(nblock-1)-1) )
  {
    comm[2] = 0;
  }

  if ( ( rank % nblock ) == 0 )
  {
    comm[3] = 0;
  }

  return;
}

void setex ( double ex[], double M[][nodeedge+2], int which )             
{
  int i;

  switch ( which ) 
  {
    case 0:
    {
      for ( i = 1; i <= nodeedge; i++) 
      {
        ex[i-1] = M[1][i];
      }
      break;
    }
    case 1:
    {
      for ( i = 1; i <= nodeedge; i++)
      {
        ex[i-1] = M[i][nodeedge];
      }
      break;
    }
    case 2:
    {
      for ( i = 1; i <= nodeedge; i++)
      {
        ex[i-1] = M[nodeedge][i];
      }
      break;
    }
    case 3:
    {
      for ( i = 1; i <= nodeedge; i++)
      {
        ex[i-1] = M[i][1];
      }
      break;
    }
  }
  return;
}

void unpack ( double M[][nodeedge+2], int where, double in[] )
{
  int i;

  if ( where == 0 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[0][i+1] = in[i]; 
    }
  }
  else if ( where == 1 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[i+1][nodeedge+1] = in[i];
    }
  }
  else if ( where == 2 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[nodeedge+1][i+1] = in[i];
    }
  }
  else if ( where == 3 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[i+1][0] = in[i];
    }
  }

  return;
}