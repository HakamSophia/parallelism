#include <stdio.h>
#include <shalw.h>

FILE *create_file(void) {
  FILE *f;
  char fname[256];

  sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), g_size_x, g_size_y, nb_steps);

  f = fopen(fname, "w+b");

  return f;
}

void export_step(FILE *f, int t) {
  fwrite((void *)&HFIL(t, 0, 0), sizeof(double), g_size_x * g_size_y, f);
} // DOIT ETRE FAIT QUE PAR LE 0 PEUT ETRE 

void finalize_export(FILE *f) {
  fclose(f);
}
