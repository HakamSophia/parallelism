#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>

//double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy, *local_HPHY, *local_VPHY, *local_UPHY, *local_UFil, *local_VFil, *local_HFil ;
double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;

int nb_steps;
int g_size_x, g_size_y;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export;
std::string export_path;

int main(int argc, char **argv) {

// 	int p;
// 	int Np;

// //  MPI_Status status;
// //  MPI_Init(NULL, NULL);
// //  MPI_Comm_size(MPI_COMM_WORLD,&Np);
// //  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
// //  MPI_Comm_rank(MPI_COMM_WORLD,&p);




 
  parse_args(argc, argv);
  printf("Command line options parsed\n");
  
  alloc();
  printf("Memory allocated\n");
  
  // a faire if rang =0 peut etre 
  
	  gauss_init(); // la premiere image au temps t=0
	  printf("State initialised\n");
  

 // printf("size_x=%d et size_y=%d\n",g_size_x,g_size_y);  
  //loc_alloc();
  printf("HELLOOOOOOOOOO 0\n");
  forward(); // va construire une image au temps d'aprés au temps t1 image t1, image t2 etc.....
  printf("State computed\n");
  


  dealloc();
  printf("Memory freed\n");
  
  return EXIT_SUCCESS;
}
