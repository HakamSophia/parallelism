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

double hPhy_forward(int t, int i, int j, int size_x, int size_y) {
  double c, d;
  
  c = 0.;
  if (i > 0)
    c = UPHY(t - 1, i - 1, j);

  d = 0.;
  if (j < size_y - 1)
    d = VPHY(t - 1, i, j + 1);

  return HFIL(t - 1, i, j) -
    dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +
		 (d - VPHY(t - 1, i, j)) / dy);
}

double uPhy_forward(int t, int i, int j, int size_x, int size_y) {
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

double vPhy_forward(int t, int i, int j, int size_x, int size_y) {
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

   /******************************************/
//int convolution(unsigned char tab[],int nbl,int nbc) {


// int convolution(int nbl,int nbc) {

//   int i,j;
//   double svdt = 0.;
//   int t=0;


//   for (t = 1; t < nb_steps; t++) {
//     if (t == 1) {
//       svdt = dt;
//       dt = 0;
//     }
//     if (t == 2){
//       dt = svdt / 2.;
//     }
//   }


//   unsigned char *tmp_UFIL, *tmp_VFIL, *tmp_HFIL, *tmp_UPHY, *tmp_VPHY, *tmp_HPHY, *tmp;
  
//   /* Allocation memoire du tampon intermediaire : */
//   tmp_UFIL = (unsigned char*) malloc(sizeof(unsigned char)*nbc*nbl);
//   tmp_VFIL = (unsigned char*) malloc(sizeof(unsigned char)*nbc*nbl);
//   tmp_HFIL = (unsigned char*) malloc(sizeof(unsigned char)*nbc*nbl);
//   tmp_UPHY = (unsigned char*) malloc(sizeof(unsigned char)*nbc*nbl);
//   tmp_VPHY = (unsigned char*) malloc(sizeof(unsigned char)*nbc*nbl);
//   tmp_HPHY = (unsigned char*) malloc(sizeof(unsigned char)*nbc*nbl);

//   if (tmp_UFIL == NULL || tmp_HPHY == NULL || tmp_VPHY == NULL || tmp_VFIL == NULL || tmp_UPHY == NULL || tmp_HFIL==NULL ) {
//     printf("Erreur dans l'allocation de tmp dans convolution \n");
//     return 1;
//   }

//   /* on laisse tomber les bords */
//   // // for(i=1 ; i<nbl-1 ; i++){
//   //   for(j=1 ; j<nbc-1 ; j++){
//   //     tmp_UFIL[i*nbc+j] = uFil_forward(t, i, j)
// 		// 	    tab[(i+1)*nbc+j-1],tab[(i+1)*nbc+j],tab[(i+1)*nbc+j+1],
// 		// 	    tab[(i  )*nbc+j-1],tab[(i)*nbc+j],tab[(i)*nbc+j+1],
// 		// 	    tab[(i-1)*nbc+j-1],tab[(i-1)*nbc+j],tab[(i-1)*nbc+j+1]);
//   //     tmp_VFIL[i*nbc+j] = vFil_forward(t, i,  j)
// 		// 	    tab[(i+1)*nbc+j-1],tab[(i+1)*nbc+j],tab[(i+1)*nbc+j+1],
// 		// 	    tab[(i  )*nbc+j-1],tab[(i)*nbc+j],tab[(i)*nbc+j+1],
// 		// 	    tab[(i-1)*nbc+j-1],tab[(i-1)*nbc+j],tab[(i-1)*nbc+j+1]);
//   //     tmp_HFIL[i*nbc+j] = hFil_forward( t,  i,  j)
// 		// 	    tab[(i+1)*nbc+j-1],tab[(i+1)*nbc+j],tab[(i+1)*nbc+j+1],
// 		// 	    tab[(i  )*nbc+j-1],tab[(i)*nbc+j],tab[(i)*nbc+j+1],
// 		// 	    tab[(i-1)*nbc+j-1],tab[(i-1)*nbc+j],tab[(i-1)*nbc+j+1]);
//   //     tmp_UPHY[i*nbc+j] = uPhy_forward(t, i, j, nbl, nbc)
// 		// 	    tab[(i+1)*nbc+j-1],tab[(i+1)*nbc+j],tab[(i+1)*nbc+j+1],
// 		// 	    tab[(i  )*nbc+j-1],tab[(i)*nbc+j],tab[(i)*nbc+j+1],
// 		// 	    tab[(i-1)*nbc+j-1],tab[(i-1)*nbc+j],tab[(i-1)*nbc+j+1]);
//   //     tmp_VPHY[i*nbc+j] = vPhy_forward( t,  i,  j, nbl, nbc)
// 		// 	    tab[(i+1)*nbc+j-1],tab[(i+1)*nbc+j],tab[(i+1)*nbc+j+1],
// 		// 	    tab[(i  )*nbc+j-1],tab[(i)*nbc+j],tab[(i)*nbc+j+1],
// 		// 	    tab[(i-1)*nbc+j-1],tab[(i-1)*nbc+j],tab[(i-1)*nbc+j+1]);
//   //     tmp_HPHY[i*nbc+j] = hPhy_forward(t, i, j,  nbl,  nbc)
// 		// 	    tab[(i+1)*nbc+j-1],tab[(i+1)*nbc+j],tab[(i+1)*nbc+j+1],
// 		// 	    tab[(i  )*nbc+j-1],tab[(i)*nbc+j],tab[(i)*nbc+j+1],
// 		// 	    tab[(i-1)*nbc+j-1],tab[(i-1)*nbc+j],tab[(i-1)*nbc+j+1]);
//   //   } /* for j */
//   // } /* for i */
  
//   /* Recopie de l'image apres traitement dans l'image initiale,
//    * On remarquera que la premiere, la derniere ligne, la premiere
//    * et la derniere colonne ne sont pas copiées (ce qui force a faire
//    * la copie ligne par ligne). */







//   // for( i=1; i<nbl-1; i++){
//   //   memcpy(tab+nbc*i+1, tmp_UFIL+nbc*i+1, (nbc-2)*sizeof(unsigned char));
//   //   memcpy(tab+nbc*i+1, tmp_VFIL+nbc*i+1, (nbc-2)*sizeof(unsigned char));
//   //   memcpy(tab+nbc*i+1, tmp_HFIL+nbc*i+1, (nbc-2)*sizeof(unsigned char));
//   //   memcpy(tab+nbc*i+1, tmp_UPHY+nbc*i+1, (nbc-2)*sizeof(unsigned char));
//   //   memcpy(tab+nbc*i+1, tmp_VPHY+nbc*i+1, (nbc-2)*sizeof(unsigned char));
//   //   memcpy(tab+nbc*i+1, tmp_HPHY+nbc*i+1, (nbc-2)*sizeof(unsigned char));
//   // } /* for i */
  
//   /* Liberation memoire du tampon intermediaire : */
//   free(tmp_UFIL); 
//   free(tmp_VFIL); 
//   free(tmp_HFIL); 
//   free(tmp_UPHY); 
//   free(tmp_VPHY); 
//   free(tmp_HPHY);  

// }
   /******************************************/

void forward() {
printf("FORWARD\n");



  //int rank;
  
  //int i;
  //int TAG_FIRST_ROW = 1;
  //int TAG_LAST_ROW = g_size_x;
printf("FORWARD1\n");





//printf("HELLOO 1");



    FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  
  if (file_export) {
    file = create_file();
    export_step(file, t);
  }

printf("FORWARD2\n");



int p, Np;


  MPI_Status status;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&Np);
  //MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_rank(MPI_COMM_WORLD,&p);




printf("FORWARD3\n");



 /* debut du chronometrage */
  double debut = MPI_Wtime();
  double fin;



  int size_x, size_y;
  // size_x = g_size_x/Np;
  // size_y = g_size_y;
size_x = ( p==0 || p==Np-1 )?(g_size_x/Np+1):(g_size_x/Np+2);
size_y = g_size_y;






  /* Envoi des dimensions de la grille à tous les proc */
MPI_Bcast(&g_size_x,1,MPI_INT,0,MPI_COMM_WORLD);
MPI_Bcast(&g_size_y,1,MPI_INT,0,MPI_COMM_WORLD);



printf("FORWARD4\n");




   /******************************************/




  // local_UFil = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
  // local_VFil = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
  // local_HFil = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
  // local_UPHY = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
  // local_VPHY = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
  // local_HPHY = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));

printf("FORWARD5\n");





/*
if (local_UFil == NULL && local_VFil && local_HFil && local_UPHY && local_VPHY && local_HPHY) {
    printf("Erreur dans l'allocation de tmp dans convolution \n");
    return 1;
  }
*/


// local_UFIL=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));
// local_VFIL=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));
// local_HFIL=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));
// local_UPHY=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));
// local_VPHY=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));
// local_HPHY=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));


  /* Envoi des parties de l'image à chaque proc */
  /* Chaque processeur dispose d'un local de taille size_y*size_x*/


//MPI_Scatter(hFil, g_size_y*g_size_x/Np, MPI_DOUBLE, local_HFil+size_y*(p!=0), size_y*g_size_x/Np, MPI_DOUBLE,0, MPI_COMM_WORLD);


MPI_Scatter(hFil, g_size_y*g_size_x/Np, MPI_DOUBLE, hFil+size_y*(p!=0), size_y*g_size_x/Np, MPI_DOUBLE,0, MPI_COMM_WORLD);




printf("FORWARD6\n");




  for (t = 1; t < nb_steps; t++) {
    if (t == 1) {
      svdt = dt;
      dt = 0;
    }
    if (t == 2){
      dt = svdt / 2.;
    }

int tag=0;
int tag1=10;

    /* Récupération et envoi des lignes à la frontière avec les proc voisins */
    if (p!=0){
		// MPI_Send(local_UFIL+size_y, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_FIRST_ROW, MPI_COMM_WORLD);
		// MPI_Send(local_VFIL+size_y, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_FIRST_ROW, MPI_COMM_WORLD);
		printf("FORWARD61\n");

		//MPI_Send(local_HFil+size_y, size_y, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);

		MPI_Send(&hFil, size_y, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);


		//MPI_Sendrecv(local_HFil+size_y*(size_x-2), size_y, MPI_UNSIGNED_CHAR,p-1,size_y,MPI_UNSIGNED_CHAR,p+1,tag, MPI_COMM_WORLD,&status);


		printf("FORWARD62\n");

		// MPI_Send(local_UPHY+size_y, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_FIRST_ROW, MPI_COMM_WORLD);
		// MPI_Send(local_VPHY+size_y, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_FIRST_ROW, MPI_COMM_WORLD);
		// MPI_Send(local_HPHY+size_y, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_FIRST_ROW, MPI_COMM_WORLD);
	 }
    if (p!=Np-1){
	   // MPI_Send(local_UFIL+size_y*(size_x-2), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_LAST_ROW, MPI_COMM_WORLD);
	   // MPI_Send(local_VFIL+size_y*(size_x-2), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_LAST_ROW, MPI_COMM_WORLD);
	    //MPI_Send(local_HFil+size_y*(size_x-2), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_LAST_ROW, MPI_COMM_WORLD);



	    		printf("FORWARD63\n");


    	//MPI_Send(local_HFil+size_y*(g_size_x/Np-2), size_y, MPI_DOUBLE, p+1, tag1, MPI_COMM_WORLD);

    	MPI_Send(hFil+size_y*(g_size_x/Np-2), size_y, MPI_DOUBLE, p+1, tag1, MPI_COMM_WORLD);


    			printf("FORWARD64\n");



		//MPI_Sendrecv(local_HFil+size_y*(size_x-2), size_y, MPI_UNSIGNED_CHAR,p+1,size_y,MPI_UNSIGNED_CHAR,p-1,tag, MPI_COMM_WORLD,&status);







	   // MPI_Send(local_UPHY+size_y*(size_x-2), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_LAST_ROW, MPI_COMM_WORLD);
	   // MPI_Send(local_VPHY+size_y*(size_x-2), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_LAST_ROW, MPI_COMM_WORLD);
	   // MPI_Send(local_HPHY+size_y*(size_x-2), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_LAST_ROW, MPI_COMM_WORLD);}
	}
		
    if (p!=Np-1){
		// MPI_Recv(local_UFIL+size_y*(size_x-1), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_FIRST_ROW, MPI_COMM_WORLD, &status);
		// MPI_Recv(local_VFIL+size_y*(size_x-1), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_FIRST_ROW, MPI_COMM_WORLD, &status);


	   //MPI_Recv(local_HFil+size_y*(size_x-1), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_FIRST_ROW, MPI_COMM_WORLD, &status);


		printf("FORWARD65\n");




	   //MPI_Recv(local_HFil+size_y*(size_x-1), size_y, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);

	   MPI_Recv(&hFil+size_y*(size_x-1), size_y, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);


		printf("FORWARD66\n");



	 //   MPI_Recv(local_UPHY+size_y*(size_x-1), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_FIRST_ROW, MPI_COMM_WORLD, &status);


	 //   MPI_Recv(local_VPHY+size_y*(size_x-1), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_FIRST_ROW, MPI_COMM_WORLD, &status);
	 //   MPI_Recv(local_HPHY+size_y*(size_x-1), size_y, MPI_UNSIGNED_CHAR, p+1, TAG_FIRST_ROW, MPI_COMM_WORLD, &status);
	 }
    if (p!=0){
		// MPI_Recv(local_UFIL, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_LAST_ROW, MPI_COMM_WORLD, &status);
		// MPI_Recv(local_VFIL, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_LAST_ROW, MPI_COMM_WORLD, &status);
		 //MPI_Recv(local_HFil, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_LAST_ROW, MPI_COMM_WORLD, &status);



    			printf("FORWARD67\n");






	   //MPI_Recv(local_HFil, size_y, MPI_DOUBLE, p-1, tag1, MPI_COMM_WORLD, &status);


	   //MPI_Recv(local_HFil, size_y, MPI_DOUBLE, p-1, tag1, MPI_COMM_WORLD, &status);

	   MPI_Recv(&hFil, size_y, MPI_DOUBLE, p-1, tag1, MPI_COMM_WORLD, &status);

		printf("FORWARD68\n");




		// MPI_Recv(local_UPHY, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_LAST_ROW, MPI_COMM_WORLD, &status);
		// MPI_Recv(local_VPHY, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_LAST_ROW, MPI_COMM_WORLD, &status);
		// MPI_Recv(local_HPHY, size_y, MPI_UNSIGNED_CHAR, p-1, TAG_LAST_ROW, MPI_COMM_WORLD, &status);
	 }


printf("FORWARD7\n");



  }

    for (int j = 0; j < size_y; j++) {
      for (int i = 0; i < size_x; i++) {
	HPHY(t, i, j) = hPhy_forward(t, i, j, size_x, size_y);
	UPHY(t, i, j) = uPhy_forward(t, i, j, size_x, size_y);
	VPHY(t, i, j) = vPhy_forward(t, i, j, size_x, size_y);
	HFIL(t, i, j) = hFil_forward(t, i, j);
	UFIL(t, i, j) = uFil_forward(t, i, j);
	VFIL(t, i, j) = vFil_forward(t, i, j);
      }
    }










	// convolution(local_UFIL, size_x, size_y);
	// convolution(local_VFIL, size_x, size_y);
	// convolution(local_HFIL, size_x, size_y);
	// convolution(local_UPHY, size_x, size_y);
	// convolution(local_VPHY, size_x, size_y);
	// convolution(local_HPHY, size_x, size_y);
//   convolution(size_x, size_y);
  // convolution(size_x, size_y);
   //convolution(size_x, size_y);
   //convolution(size_x, size_y);
   //convolution(size_x, size_y);
   //convolution(size_x, size_y);
	

  /* Rassembler les images locale dans l'image globale chez le proc root */
  //MPI_Gather(local_UFIL+w*(p!=0), size_y*g_size_x/Np, MPI_UNSIGNED_CHAR, file_export, size_y*g_size_x/Np, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  //MPI_Gather(local_VFIL+w*(p!=0), size_y*g_size_x/Np, MPI_UNSIGNED_CHAR, file_export, size_y*g_size_x/Np, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  



  //MPI_Gather(local_HFil+size_y*(p!=0), size_y*g_size_x/Np, MPI_DOUBLE, hFil, size_y*g_size_x/Np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(hFil+size_y*(p!=0), size_y*g_size_x/Np, MPI_DOUBLE, hFil, size_y*g_size_x/Np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  


  //MPI_Gather(local_UPHY+w*(p!=0), size_y*g_size_x/Np, MPI_UNSIGNED_CHAR, file_export, size_y*g_size_x/Np, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  //MPI_Gather(local_VPHY+w*(p!=0), size_y*g_size_x/Np, MPI_UNSIGNED_CHAR, file_export, size_y*g_size_x/Np, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  //MPI_Gather(local_HPHY+w*(p!=0), size_y*g_size_x/Np, MPI_UNSIGNED_CHAR, file_export, size_y*g_size_x/Np, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  // free(local_UFil);
  // free(local_VFil);
  // free(local_HFil);
  // free(local_UPHY);
  // free(local_VPHY);
  // free(local_HPHY);

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
