#include <string>
extern double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy, *local_HPHY, *local_VPHY, *local_UPHY, *local_UFil, *local_VFil, *local_HFil;
//extern int size_x, size_y, nb_steps;
extern int nb_steps;
extern int g_size_y, g_size_x;
extern double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
extern bool file_export;
extern int Np;
extern std::string export_path;



#define HFIL(t, i, j) hFil[ (j) +			\
			    (i) * g_size_y +		\
			    ((t)%2) * g_size_x * g_size_y ]
#define UFIL(t, i, j) uFil[ (j) +			\
			    (i) * g_size_y +		\
			    ((t)%2) * g_size_x * g_size_y ]
#define VFIL(t, i, j) vFil[ (j) +			\
			    (i) * g_size_y +		\
			    ((t)%2) * g_size_x * g_size_y ]
#define HPHY(t, i, j) hPhy[ (j) +			\
			    (i) * g_size_y +		\
			    ((t)%2) * g_size_x * g_size_y ]
#define UPHY(t, i, j) uPhy[ (j) +			\
			    (i) * g_size_y +		\
			    ((t)%2) * g_size_x * g_size_y ]
#define VPHY(t, i, j) vPhy[ (j) +			\
			    (i) * g_size_y +		\
			    ((t)%2) * g_size_x * g_size_y ]
