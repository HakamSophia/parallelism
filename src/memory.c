#include <stdlib.h>
#include <shalw.h>








void alloc(void) {
  hFil = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
  uFil = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
  vFil = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
  hPhy = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
  uPhy = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
  vPhy = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
}

// void loc_alloc(void){
//  /* Allocation memoire du tampon intermediaire : */
//   local_UFil = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
//   local_VFil = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
//   local_HFil = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
//   local_UPHY = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
//   local_VPHY = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
//   local_HPHY = (double *) calloc(2*g_size_x/Np*g_size_y, sizeof(double));
// }
void dealloc(void) {
  free(hFil);
  free(uFil);
  free(vFil);
  free(hPhy);
  free(uPhy);
  free(vPhy);
}
