#include "bem2_emf_qte.h"

int main(int argc,char *argv[]) 
{
  DOMD md;
  FILE *fp1,*fp2;
  double complex ez,hr[2];
  double rang,dr,r[2],*ie,*ih;
  int max,i,j,type;

  dat_read(argv[1],&md);  // read data file outputed by d2qte_bv_solver
  print_data(&md);        // print data
  
  max=200;                // number of sampling 
  rang=1.5*md.wd.lambda0; // range of sampling
  dr=rang*2.0/(double)(max-1);
  type=1;                 // type=1 : 7 point GK
  
  ie=(double *)m_alloc2(max,sizeof(double),"exampl2.c,ie");
  ih=(double *)m_alloc2(max,sizeof(double),"exampl2.c,ih");
  
  if((fp1=fopen("Ie_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  } 
  fprintf(fp1,"%s\n","# x y electric_field_intensity");
  if((fp2=fopen("Ih_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  } 
  fprintf(fp2,"%s\n","# x y magnetic_field_intensity");
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    #pragma omp parallel for schedule(dynamic) firstprivate(r) private(ez,hr) // omp parallel
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      EH_t(&ez,hr,r,type,&md);
      ie[j]=creal(ez*conj(ez));
      ih[j]=creal(hr[0]*conj(hr[0]))+creal(hr[1]*conj(hr[1]));
    }
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      fprintf(fp1,"%g %g %15.14e\n",r[0],r[1],ie[j]);
      fprintf(fp2,"%g %g %15.14e\n",r[0],r[1],ih[j]);
    }
    fprintf(fp1,"\n");
    fprintf(fp2,"\n");
  }
  fclose(fp1);
  fclose(fp2);
  
  printf("Intensity plot is finished\n");
  
  free(ie);
  free(ih);
  mfree_domd(&md);
  return 0;
}
