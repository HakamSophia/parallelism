#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
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

void forward(int argc, char *argv[]) {



  int Np; 
  int rank;
  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&Np);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int tag=0;




  FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  
  if (file_export) {
    file = create_file();
    export_step(file, t);
  }
    
//size_x et size_y définis localement 
int size_x, size_y;
size_x=g_size_x/Np;
size_y=g_size_y;

//On envoie la taille 

MPI_Bcast(&size_x,1,MPI_INT,0,MPI_COMM_WORLD);
MPI_Bcast(&size_y,1,MPI_INT,0,MPI_COMM_WORLD);


local_UFil=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));
local_VFil=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));
local_HFil=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));
local_UPHY=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));
local_VPHY=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));
local_HPHY=(unsigned char*)malloc(size_x*size_y*sizeof(unsigned char));





//MPI_Scatter(hFil, size_y*size_x, MPI_UNSIGNED_CHAR, hfil+size_y*(rank!=0), size_y*size_x, MPI_UNSIGNED_CHAR,0, MPI_COMM_WORLD);
//MPI_Scatter(uFil, size_y*size_x, MPI_UNSIGNED_CHAR, hfil+size_y*(rank!=0), size_y*size_x, MPI_UNSIGNED_CHAR,0, MPI_COMM_WORLD);
//MPI_Scatter(vFil, size_y*size_x, MPI_UNSIGNED_CHAR, hfil+size_y*(rank!=0), size_y*size_x, MPI_UNSIGNED_CHAR,0, MPI_COMM_WORLD);
//MPI_Scatter(hPhy, size_y*size_x, MPI_UNSIGNED_CHAR, hfil+size_y*(rank!=0), size_y*size_x, MPI_UNSIGNED_CHAR,0, MPI_COMM_WORLD);
//MPI_Scatter(vPhy, size_y*size_x, MPI_UNSIGNED_CHAR, hfil+size_y*(rank!=0), size_y*size_x, MPI_UNSIGNED_CHAR,0, MPI_COMM_WORLD);
//MPI_Scatter(uPhy, size_y*size_x, MPI_UNSIGNED_CHAR, hfil+size_y*(rank!=0), size_y*size_x, MPI_UNSIGNED_CHAR,0, MPI_COMM_WORLD);




  for (t = 1; t < nb_steps; t++) {
    if (t == 1) {
      svdt = dt;
      dt = 0;
    }
    if (t == 2){
      dt = svdt / 2.;
    }



  //int local_h=(rank==0||rank==(Np-1))?(size_x+1):(size_x+2);
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




    for (int j = 0; j < size_y; j++) {
      for (int i = 0; i < size_x; i++) {


  HPHY(t, i, j) = hPhy_forward(t, i, j);
  UPHY(t, i, j) = uPhy_forward(t, i, j);
  VPHY(t, i, j) = vPhy_forward(t, i, j);
  HFIL(t, i, j) = hFil_forward(t, i, j);
  UFIL(t, i, j) = uFil_forward(t, i, j);
  VFIL(t, i, j) = vFil_forward(t, i, j);

      }
    }




}



    for (nb=0; nb<nb_steps; nb++){
      if()
    }

    if (file_export) {
      export_step(file, t);
    }
    
    if (t == 2) {
      dt = svdt;
    }
  }

  if (file_export) {
    finalize_export(file);
  }





    if(rank==0){ 


    char nom_sortie[100] = "";
    sprintf(nom_sortie, "post-convolution_para_filtre%d_nbIter%d.ras", filtre, nbiter);
    sauve_rasterfile(nom_sortie, &r);
  }

  MPI_Finalize();

}