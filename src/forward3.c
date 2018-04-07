#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>

#include <stdlib.h>
#include <string.h> /* pour le memcpy */
#include <time.h>   /* chronometrage */
#include <mpi.h>


double hFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return HPHY(t, i, j);
  return HPHY(t - 1, i, j) +
    alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

double uFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return UPHY(t, i, j);
  return UPHY(t - 1, i, j) +
    alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

double vFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return VPHY(t, i, j);
  return VPHY(t - 1, i, j) +
    alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

double hPhy_forward(int t, int i, int j) {
  double c, d;
  
  c = 0.;
  if (i > 0)
    c = UPHY(t - 1, i - 1, j);

  d = 0.;
  if (j < g_size_y - 1)
    d = VPHY(t - 1, i, j + 1);

  return HFIL(t - 1, i, j) -
    dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +
		 (d - VPHY(t - 1, i, j)) / dy);
}

double uPhy_forward(int t, int i, int j) {
  double b, e, f, g;
  
  if (i == g_size_x - 1)
    return 0.;

  b = 0.;
  if (i < g_size_x - 1)
    b = HPHY(t - 1, i + 1, j);

  e = 0.;
  if (j < g_size_y - 1)
    e = VPHY(t - 1, i, j + 1);

  f = 0.;
  if (i < g_size_x - 1)
    f = VPHY(t - 1, i + 1, j);

  g = 0.;
  if (i < g_size_x - 1 && j < g_size_y - 1)
    g = VPHY(t - 1, i + 1, j + 1);

  return UFIL(t - 1, i, j) +
    dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
	  (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
	  (dissip * UFIL(t - 1, i, j)));
}

double vPhy_forward(int t, int i, int j) {
  double c, d, e, f;

  if (j == 0)
    return 0.;

  c = 0.;
  if (j > 0)
    c = HPHY(t - 1, i, j - 1);

  d = 0.;
  if (i > 0 && j > 0)
    d = UPHY(t - 1, i -1, j -1);

  e = 0.;
  if (i > 0)
    e = UPHY(t - 1, i - 1, j);

  f = 0.;
  if (j > 0)
    f = UPHY(t - 1, i, j - 1);

  return VFIL(t - 1, i, j) +
    dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
	  (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
	  (dissip * VFIL(t - 1, i, j)));
}


void forward() {
printf("FORWARD\n");


    FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  
  if (file_export) {
    file = create_file();
    export_step(file, t);
  }

printf("FORWARD1\n");


  
  int p, Np;
  MPI_Status status;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&Np);
  //MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_rank(MPI_COMM_WORLD,&p);
  int tag=0;
  int tag1=tag;




printf("FORWARD2\n");



 /* debut du chronometrage */
  double debut = MPI_Wtime();
  double fin;


  /* Envoi des dimensions de la grille à tous les proc */
MPI_Bcast(&g_size_x,1,MPI_INT,0,MPI_COMM_WORLD);
MPI_Bcast(&g_size_y,1,MPI_INT,0,MPI_COMM_WORLD);



printf("FORWARD3\n");

int size_x, size_y;
size_x=g_size_x/Np;
//size_x = ( p==0 || p==Np-1 )?(g_size_x/Np+1):(g_size_x/Np+2);
size_y = g_size_y;






//Dans le cas t=0:






  for (t = 1; t < nb_steps; t++) {

  	//MPI_Scatter(hFil, g_size_y*g_size_x/Np, MPI_DOUBLE, hFil+size_y*(p!=0), size_y*g_size_x/Np, MPI_DOUBLE,0, MPI_COMM_WORLD);

    if (t == 1) {
      svdt = dt;
      dt = 0;
    }
    if (t == 2){
      dt = svdt / 2.;
    }


    for (int j = 0; j < g_size_y; j++) {
      for (int i = 0; i < g_size_x; i++) {


    if (p>0){
		printf("FORWARD31\n");

		MPI_Send(&HFIL(t,i,j)+2*size_y, 2*size_y, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
    MPI_Send(&VFIL(t,i,j)+2*size_y, 2*size_y, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
    MPI_Send(&UFIL(t,i,j)+2*size_y, 2*size_y, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
    MPI_Send(&UPHY(t,i,j)+2*size_y, 2*size_y, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
    MPI_Send(&VPHY(t,i,j)+2*size_y, 2*size_y, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
    MPI_Send(&HPHY(t,i,j)+2*size_y, 2*size_y, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);




		printf("FORWARD32\n");




    printf("FORWARD37\n");


    MPI_Recv(&HFIL(t,i,j), 2*size_y, MPI_DOUBLE, p-1, tag1, MPI_COMM_WORLD, &status);
    MPI_Recv(&VFIL(t,i,j), 2*size_y, MPI_DOUBLE, p-1, tag1, MPI_COMM_WORLD, &status);
    MPI_Recv(&UFIL(t,i,j), 2*size_y, MPI_DOUBLE, p-1, tag1, MPI_COMM_WORLD, &status);
    MPI_Recv(&UPHY(t,i,j), 2*size_y, MPI_DOUBLE, p-1, tag1, MPI_COMM_WORLD, &status);
    MPI_Recv(&VPHY(t,i,j), 2*size_y, MPI_DOUBLE, p-1, tag1, MPI_COMM_WORLD, &status);
    MPI_Recv(&HPHY(t,i,j), 2*size_y, MPI_DOUBLE, p-1, tag1, MPI_COMM_WORLD, &status);


    printf("FORWARD38\n");







	 }
    if (p<Np-1){

	    printf("FORWARD33\n");

    	MPI_Send(&HFIL(t,i,j)+(size_x-2)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag1, MPI_COMM_WORLD);
      MPI_Send(&VFIL(t,i,j)+(size_x-2)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag1, MPI_COMM_WORLD);
      MPI_Send(&UFIL(t,i,j)+(size_x-2)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag1, MPI_COMM_WORLD);
      MPI_Send(&UPHY(t,i,j)+(size_x-2)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag1, MPI_COMM_WORLD);
      MPI_Send(&VPHY(t,i,j)+(size_x-2)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag1, MPI_COMM_WORLD);
      MPI_Send(&HPHY(t,i,j)+(size_x-2)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag1, MPI_COMM_WORLD);
	
    
		printf("FORWARD35\n");

	   MPI_Recv(&HFIL(t,i,j)+(size_x-1)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
     MPI_Recv(&VFIL(t,i,j)+(size_x-1)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
     MPI_Recv(&UFIL(t,i,j)+(size_x-1)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
     MPI_Recv(&UPHY(t,i,j)+(size_x-1)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
     MPI_Recv(&VPHY(t,i,j)+(size_x-1)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
     MPI_Recv(&HPHY(t,i,j)+(size_x-1)*2*size_y, 2*size_y, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);


		printf("FORWARD36\n");
}
  


	HFIL(t, i, j) = hFil_forward(t, i, j);
	VFIL(t, i, j) = vFil_forward(t, i, j);
  UFIL(t, i, j) = uFil_forward(t, i, j);

  HPHY(t, i, j) = hPhy_forward(t, i, j);
  UPHY(t, i, j) = uPhy_forward(t, i, j);
  VPHY(t, i, j) = vPhy_forward(t, i, j);


      }
    } 


}
      MPI_Gather(hFil+size_y*(p!=0), size_y*g_size_x/Np, MPI_DOUBLE, hFil, size_y*g_size_x/Np, MPI_DOUBLE, 0, MPI_COMM_WORLD);


 

  /* fin du chronometrage */
  fin = MPI_Wtime();
  fprintf( stderr, "Temps total de calcul : %g sec\n rang = %d \n",  fin - debut, p);
  /* Fin de la parallélisation */  
  MPI_Finalize();
  /*//int local_h=(rank==0||rank==(Np-1))?(size_x+1):(size_x+2);
  //double *hfil_loc=(unsigned char*)malloc(size_y*local_h*sizeof(unsigned char));;



  
  //Envoyer a chacun son travail
  
  if(rank > 0){
	 //MPISendrecv(lastline,1,MPI_UNSIGNED_CHAR, rank,tag,rank+1,1, MPI_UNSIGNED_CHAR, rank+1,tag,MPICOMMWORLD,&status) ;   
    //MPI_Send(local_data+w,w,MPI_UNSIGNED_CHAR,rank-1,tag,MPI_COMM_WORLD);
    //MPI_Recv(local_data,w,MPI_UNSIGNED_CHAR,rank-1,tag,MPI_COMM_WORLD, &status);

		MPISendrecv(&HPHY, size_y*size_x, MPI_UNSIGNED_CHAR,rank+1,tag,size_x*size_y,MPI_UNSIGNED_CHAR,rank-1,tag, MPI_COMM_WORLD,&status);
  }



  if (rank < Np-1){
//MPI_Send(local_data+(local_h-2)*w,w, MPI_UNSIGNED_CHAR,rank+1,tag,MPI_COMM_WORLD);
//MPI_Recv(local_data+w*(local_h-1), w,MPI_UNSIGNED_CHAR,rank+1,tag,MPI_COMM_WORLD,&status);
		MPISendrecv(&HFIL, size_y*size_x, MPI_UNSIGNED_CHAR,rank+1,tag,size_x*size_y,MPI_UNSIGNED_CHAR,rank-1,tag, MPI_COMM_WORLD,&status);
  }
*/


	/******************************************/

    // if (file_export) {
    //   export_step(file, t);
    // }
    
    if (t == 2) {
      dt = svdt;
    }
  

  if (file_export) {
    finalize_export(file);
  }

}
