#include "bem2_emf_qte.h"


void read_data(int argc,char **argv,DOMD *md)
{
  void filename_chk(int argc,char **argv);
  void read_medium_data(char *med_fn,DOMD *md);
  void read_mesh_data(char *msh_fn,DOMD *md);
  void read_periodicity_data(char *prd_fn,DOMD *md);

  filename_chk(argc,argv); 
  read_infd(argv[1],&(md->wd));
  read_medium_data(argv[2],md);
  read_mesh_data(argv[3],md);
  read_periodicity_data(argv[4],md);
  
  // rotation and translation data
  if(argc==9){
    md->ra=atof(argv[6]);
    md->tx=atof(argv[7]);
    md->ty=atof(argv[8]);
  }
  else {
    md->ra=0.0;
    md->tx=0.0;
    md->ty=0.0;    
  }
}

void print_data(DOMD *md)
{
  void print_medium_data(DOMD *md);
  void print_mesh_data(DOMD *md);
  void print_periodicity_data(DOMD *md);

  print_infd(&(md->wd));
  print_medium_data(md);
  print_mesh_data(md);
  print_periodicity_data(md);
  printf("\n");
  
  if(md->ra!=0.0 || md->tx!=0.0 || md->ty!=0.0){
    printf("-- rotation and translation settings --\n");
    printf("rotation angle                    [rad]: %8.7g\n",md->ra);
    printf("x-component of translation vector      : %8.7g\n",md->tx);
    printf("y-component of translation vector      : %8.7g\n",md->ty);
    printf("\n");
  }  
}

void print_data_MKSA(DOMD *md)
{
  void print_medium_data(DOMD *md);
  void print_mesh_data(DOMD *md);
  void print_periodicity_data_MKSA(DOMD *md);
  
  print_infd_MKSA(&(md->wd));
  print_medium_data(md);
  print_mesh_data(md);
  print_periodicity_data_MKSA(md);
  printf("\n");
  
  if(md->ra!=0.0 || md->tx!=0.0 || md->ty!=0.0){
    printf("-- rotation and translation settings --\n");
    printf("rotation angle                    [rad]: %8.7g\n",md->ra);
    printf("x-component of translation vector   [m]: %8.7g\n",OSUtoMKSA_length(md->tx));
    printf("y-component of translation vector   [m]: %8.7g\n",OSUtoMKSA_length(md->ty));
    printf("\n");
  }  
}

void initialize_domd(DOMD *md)
{
  void init_elem_const(BOUD *bd);
  void init_sub_domain(DOMD *md);
  void malloc_sub_domain(DOMD *md);
  void init_boundary_data(DOMD *md);
  void rotation_translation_obj(DOMD *md);
  
  // medium 
  md->n[0]=md->wd.ne;
  // rotation and translation
  if(md->ra!=0.0 || md->tx!=0.0 || md->ty!=0.0) rotation_translation_obj(md);
  // element constant 
  init_elem_const(&(md->bd));
  // sub domain
  malloc_sub_domain(md);  
  init_sub_domain(md);
  // boundary data
  init_boundary_data(md);
  // d2_qpgf_d1 data
  md->qd.k=md->wd.ne*md->wd.k0;
  md->qd.kx=md->qd.k*md->wd.kxh;
  if(md->PD==1){
    if(md->qd.k==md->qd.kx){
      printf("The x-component of wave vector must not be equal to the wave number.\n");
      printf("There is a contradiction between the wave source and the object. Exit...\n"); 
      exit(1);
    }
    if(md->qd.kx==0.0 || md->qd.d==0.0){
      printf("The quasi-periodic boundary condition is not satisfied. Exit...\n");
      exit(1);
    }
  }  
}

void mfree_domd(DOMD *md)
{
  void mfree_sub_domain(DOMD *md);
  void mfree_elem(BOUD *bd);
  void mfree_node(BOUD *bd);
  
  mfree_sub_domain(md);
  mfree_elem(&(md->bd));
  mfree_node(&(md->bd));
  free(md->n);
}

int domain_id(double *rt,DOMD *md)
{
  int domain_crt(int id,double *rt,DOMD *md);
  
  int id,re;

  id=0;
  re=domain_crt(id,rt,md);
  if(re==0) return id;
  else {
    for(id=1;id<=md->MN;id++){
      re=domain_crt(id,rt,md);
      if(re==1) return id;
    }
    printf("domain_id(), Failed to determine the domain id\n");
    printf("rt=(%15.14e, %15.14e)\n",rt[0],rt[1]);
    printf("Exit\n");
    exit(1);
  }
}

void dat_read (char *fname,DOMD *md)
{
  void malloc_node(BOUD *bd);
  void malloc_elem(BOUD *bd);
  void malloc_sub_domain(DOMD *md);
  
  FILE *fp;
  int i,j;

  if((fp=fopen(fname,"rb"))==NULL){     printf("Failed to open the %s file.\n",fname);    exit(1); }
  // fname
  if(fread(md->med_fn,sizeof(char),128,fp)!=128){
    printf("d2qte_setup.c, dat_read(), failed to read the med_fn. exit...\n");
    exit(1);
  }
  if(fread(md->msh_fn,sizeof(char),128,fp)!=128){
    printf("d2qte_setup.c, dat_read(), failed to read the msh_fn. exit...\n");
    exit(1);
  }
  if(fread(md->prd_fn,sizeof(char),128,fp)!=128){
    printf("d2qte_setup.c, dat_read(), failed to read the prd_fn. exit...\n");
    exit(1);
  }
  // rotation and translation data
  if(fread(&(md->ra),sizeof(double),1,fp)!=1){
    printf("d2qte_setup.c, dat_read(), failed to read the ra. exit...\n");
    exit(1);
  }
  if(fread(&(md->tx),sizeof(double),1,fp)!=1){
    printf("d2qte_setup.c, dat_read(), failed to read the tx. exit...\n");
    exit(1);
  }
  if(fread(&(md->ty),sizeof(double),1,fp)!=1){
    printf("d2qte_setup.c, dat_read(), failed to read the med_ty. exit...\n");
    exit(1);
  }
  // material def
  if(fread(&(md->MN),sizeof(int),1,fp)!=1){
    printf("d2qte_setup.c, dat_read(), failed to read the MN. exit...\n");
    exit(1);
  }
  md->n=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"dat_read(),md->n"); // malloc 
  if(fread(md->n,sizeof(double complex),md->MN+1,fp)!=md->MN+1){
    printf("d2qte_setup.c, dat_read(), failed to read the n. exit...\n");
    exit(1);
  }
  if(fread(&(md->PD),sizeof(int),1,fp)!=1){
    printf("d2qte_setup.c, dat_read(), failed to read the PD. exit...\n");
    exit(1);
  }
  // incident field def
  if(fread(&(md->wd),sizeof(INFD),1,fp)!=1){
    printf("d2qte_setup.c, dat_read(), failed to read the wd. exit...\n");
    exit(1);
  }
  // hm_qgf def
  if(fread(&(md->qd),sizeof(QPDT),1,fp)!=1){
    printf("d2qte_setup.c, dat_read(), failed to read the qd. exit...\n");
    exit(1);
  }
  // boundary def
  if(fread(&(md->bd.Nn),sizeof(int),1,fp)!=1){
    printf("d2qte_setup.c, dat_read(), failed to read the Nn. exit...\n");
    exit(1);
  }
  if(fread(&(md->bd.Ne),sizeof(int),1,fp)!=1){
    printf("d2qte_setup.c, dat_read(), failed to read the Ne. exit...\n");
    exit(1);
  }
  malloc_node(&(md->bd));  malloc_elem(&(md->bd)); // malloc
  if(fread(md->bd.et,sizeof(double),7,fp)!=7){
    printf("d2qte_setup.c, dat_read(), failed to read the et. exit...\n");
    exit(1);
  }
  if(fread(md->bd.wg,sizeof(double),7,fp)!=7){
    printf("d2qte_setup.c, dat_read(), failed to read the wg. exit...\n");
    exit(1);
  }  
  if(fread(md->bd.wk,sizeof(double),7,fp)!=7){
    printf("d2qte_setup.c, dat_read(), failed to read the wk. exit...\n");
    exit(1);
  }
  if(fread(md->bd.M,sizeof(double),21,fp)!=21){
    printf("d2qte_setup.c, dat_read(), failed to read the M. exit...\n");
    exit(1);
  }
  if(fread(md->bd.xn,sizeof(double),md->bd.Nn+1,fp)!=md->bd.Nn+1){
    printf("d2qte_setup.c, dat_read(), failed to read the xn. exit...\n");
    exit(1);
  }  
  if(fread(md->bd.yn,sizeof(double),md->bd.Nn+1,fp)!=md->bd.Nn+1){
    printf("d2qte_setup.c, dat_read(), failed to read the yn. exit...\n");
    exit(1);
  }
  if(fread(md->bd.md, sizeof(int),md->bd.Ne+1,fp)!=md->bd.Ne+1){
    printf("d2qte_setup.c, dat_read(), failed to read the md. exit...\n");
    exit(1);
  }  
  if(fread(md->bd.sd, sizeof(int),md->bd.Ne+1,fp)!=md->bd.Ne+1){
    printf("d2qte_setup.c, dat_read(), failed to read the sd. exit...\n");
    exit(1);
  }
  if(fread(md->bd.egd,sizeof(int),md->bd.Ne+1,fp)!=md->bd.Ne+1){
    printf("d2qte_setup.c, dat_read(), failed to read the egd. exit...\n");
    exit(1);
  }
  if(fread(md->bd.xh,sizeof(double),GLH,fp)!=GLH){
    printf("d2qte_setup.c, dat_read(), failed to read the xh. exit...\n");
    exit(1);
  }  
  if(fread(md->bd.wh,sizeof(double),GLH,fp)!=GLH){
    printf("d2qte_setup.c, dat_read(), failed to read the wh. exit...\n");
    exit(1);
  }
  for(i=0;i<=md->bd.Ne;i++){
    if(fread(md->bd.ed[i],sizeof(int),3,fp)!=3){
      printf("d2qte_setup.c, dat_read(), failed to read the ed[i]. exit...\n");
      exit(1);
    }
  }
  for(i=0;i<=md->bd.Ne;i++){
    if(fread(md->bd.x[i],sizeof(double),7,fp)!=7){
      printf("d2qte_setup.c, dat_read(), failed to read the x[i]. exit...\n");
      exit(1);
    }
  }
  for(i=0;i<=md->bd.Ne;i++){
    if(fread(md->bd.y[i],sizeof(double),7,fp)!=7){
      printf("d2qte_setup.c, dat_read(), failed to read the y[i]. exit...\n");
      exit(1);
    }
  }
  for(i=0;i<=md->bd.Ne;i++){
    if(fread(md->bd.dx[i],sizeof(double),7,fp)!=7){
      printf("d2qte_setup.c, dat_read(), failed to read the dx[i]. exit...\n");
      exit(1);
    }
  }
  for(i=0;i<=md->bd.Ne;i++){
    if(fread(md->bd.dy[i],sizeof(double),7,fp)!=7){
      printf("d2qte_setup.c, dat_read(), failed to read the dy[i]. exit...\n");
      exit(1);
    }
  }
  for(i=0;i<=md->bd.Ne;i++){
    if(fread(md->bd.ui[i],sizeof(double complex),3,fp)!=3){
      printf("d2qte_setup.c, dat_read(), failed to read the ui[i]. exit...\n");
      exit(1);
    }
  }
  for(i=0;i<=md->bd.Ne;i++){
    if(fread(md->bd.duidn[i],sizeof(double complex),3,fp)!=3){
      printf("d2qte_setup.c, dat_read(), failed to read the duidn[i]. exit...\n");
      exit(1);
    }
  }
  // sub domain def
  malloc_sub_domain(md);
  for(i=0;i<=md->MN;i++){
    if(fread(md->bd.sb[i].eid,sizeof(int),md->bd.sb[i].Ne+1,fp)!=md->bd.sb[i].Ne+1){
      printf("d2qte_setup.c, dat_read(), failed to read the sb[i]. exit...\n");
      exit(1);
    }
    for(j=0;j<=md->bd.sb[i].Ne;j++){
      if(fread(md->bd.sb[i].u[j],sizeof(double complex),3,fp)!=3){
        printf("d2qte_setup.c, dat_read(), failed to read the u[j]. exit...\n");
        exit(1);
      }
      if(fread(md->bd.sb[i].v[j],sizeof(double complex),3,fp)!=3){
        printf("d2qte_setup.c, dat_read(), failed to read the v[j]. exit...\n");
        exit(1);
      }
      if(fread(md->bd.sb[i].w[j],sizeof(double complex),3,fp)!=3){
        printf("d2qte_setup.c, dat_read(), failed to read the w[j]. exit...\n");
        exit(1);
      }
      if(fread(md->bd.sb[i].dudn[j],sizeof(double complex),3,fp)!=3){
        printf("d2qte_setup.c, dat_read(), failed to read the dudn[j]. exit...\n");
        exit(1);
      }
      if(fread(md->bd.sb[i].dvdn[j],sizeof(double complex),3,fp)!=3){
        printf("d2qte_setup.c, dat_read(), failed to read the dvdn[j]. exit...\n");
        exit(1);
      }
      if(fread(md->bd.sb[i].dwdn[j],sizeof(double complex),3,fp)!=3){
        printf("d2qte_setup.c, dat_read(), failed to read the dwdn[j]. exit...\n");
        exit(1);
      }
    }
  }
  fclose(fp);
}

void dat_write(char *fname,DOMD *md)
{
  FILE *fp;
  int i,j;

  if((fp=fopen(fname,"wb"))==NULL){    printf("Failed to open the %s file.\n",fname);    exit(1);  }
  // fname
  fwrite(md->med_fn,sizeof(char),128,fp);
  fwrite(md->msh_fn,sizeof(char),128,fp);
  fwrite(md->prd_fn,sizeof(char),128,fp);
  // rotation and translation data
  fwrite(&(md->ra),sizeof(double),1,fp);
  fwrite(&(md->tx),sizeof(double),1,fp);
  fwrite(&(md->ty),sizeof(double),1,fp);
  // material def
  fwrite(&(md->MN),sizeof(int),1,fp);
  fwrite(md->n,sizeof(double _Complex),md->MN+1,fp);
  fwrite(&(md->PD),sizeof(int),1,fp);
  // incident field def
  fwrite(&(md->wd),sizeof(INFD),1,fp);
  // hm_qgf def
  fwrite(&(md->qd),sizeof(QPDT),1,fp);
  // boundary def
  fwrite(&(md->bd.Nn),sizeof(int),1,fp);
  fwrite(&(md->bd.Ne),sizeof(int),1,fp);
  fwrite(md->bd.et,sizeof(double),7,fp);
  fwrite(md->bd.wg,sizeof(double),7,fp);  fwrite(md->bd.wk,sizeof(double),7,fp);
  fwrite(md->bd.M,sizeof(double),21,fp);
  fwrite(md->bd.xn,sizeof(double),md->bd.Nn+1,fp);  fwrite(md->bd.yn,sizeof(double),md->bd.Nn+1,fp);
  fwrite(md->bd.md, sizeof(int),md->bd.Ne+1,fp);  fwrite(md->bd.sd, sizeof(int),md->bd.Ne+1,fp);
  fwrite(md->bd.egd,sizeof(int),md->bd.Ne+1,fp);
  fwrite(md->bd.xh,sizeof(double),GLH,fp);  fwrite(md->bd.wh,sizeof(double),GLH,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.ed[i],sizeof(int),3,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.x[i],sizeof(double),7,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.y[i],sizeof(double),7,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.dx[i],sizeof(double),7,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.dy[i],sizeof(double),7,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.ui[i],sizeof(double _Complex),3,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.duidn[i],sizeof(double _Complex),3,fp);
  // sub domain def
  for(i=0;i<=md->MN;i++){
    fwrite(md->bd.sb[i].eid,sizeof(int),md->bd.sb[i].Ne+1,fp);
    for(j=0;j<=md->bd.sb[i].Ne;j++){
      fwrite(md->bd.sb[i].u[j],sizeof(double _Complex),3,fp);
      fwrite(md->bd.sb[i].v[j],sizeof(double _Complex),3,fp);
      fwrite(md->bd.sb[i].w[j],sizeof(double _Complex),3,fp);
      fwrite(md->bd.sb[i].dudn[j],sizeof(double _Complex),3,fp);
      fwrite(md->bd.sb[i].dvdn[j],sizeof(double _Complex),3,fp);
      fwrite(md->bd.sb[i].dwdn[j],sizeof(double _Complex),3,fp);
    }
  }
  fclose(fp);
}

void output_node_particles(char *fname,DOMD *md)
{
  FILE *fp;
  int s1,i,j;
  char *sd,fo[256]="",tf[200]="";

  s1=strlen(fname);
  if(s1>200){
    printf("d2qte_setup.c, output_node_particles(), file name is too long. exit...\n");
    exit(1);
  }
  sprintf(fo,"%s",fname);
  sd=strrchr(fo,'.');
  if(sd!=NULL){
    strncpy(tf,fname,s1-strlen(sd));
    sprintf(fo,"%s.particles",tf);
  }
  
  if((fp=fopen(fo,"wt"))==NULL){    printf("Can not open the %s file.\n",fo);    exit(1);  }
  fprintf(fp,"# x y object_id\n");
  
  for(i=1;i<=md->bd.Ne;i++){
    for(j=0;j<3;j++){
      fprintf(fp,"%15.14e %15.14e %d\n",md->bd.x[i][j],md->bd.y[i][j],0);
    }
  }

  fclose(fp);
}

//////////////////////////////////////////////////////////////////////
void filename_chk(int argc,char **argv)
{
  if(argc!=6 && argc!=9){     
    printf("This program needs command line arguments as follows.\n");
    printf("%s incident_field_datafile_name medium_datafile_name mesh_datafile_name periodicity_datafile_name output_datafile_name [angle trans_x trans_y](optional)\n",argv[0]);
    exit(0);
  }
}

void read_medium_data(char *med_fn,DOMD *md)
{
  FILE *fp;
  double td,td2;
  char buf[256]="";
  int i,ti;

  if((fp=fopen(med_fn,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",med_fn);    exit(1);  }
  strcpy(md->med_fn,med_fn);
  if(fgets(buf,256,fp)==NULL){
    printf("d2qte_setup.c, read_medium_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("d2qte_setup.c, read_medium_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%d\n",&ti)!=1){
    printf("d2qte_setup.c, read_medium_data(), failed to read the MN. exit...\n");
    exit(1);
  }
  md->MN=ti;
  md->n=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"read_medium_data(), md->n");
  if(fgets(buf,256,fp)==NULL){
    printf("d2qte_setup.c, read_medium_data(), failed to read the line. exit...\n");
    exit(1);
  }
  for(i=1;i<=md->MN;i++){
    if(fscanf(fp,"%lf",&td)!=1){
      printf("d2qte_setup.c, read_medium_data(), failed to read the rean(n). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&td2)!=1){
      printf("d2qte_setup.c, read_medium_data(), failed to read the imag(n). exit...\n");
      exit(1);
    }
    md->n[i]=td+td2*I;
  }
  fclose(fp);
}

void print_medium_data(DOMD *md)
{
  int i;
  printf("-- medium data --\n");
  printf("medium data file name                  : %s\n",md->med_fn);
  for(i=1;i<=md->MN;i++){
    printf("medium (domain) id %2d refractive index :%8.7g + %8.7gI\n",i,creal(md->n[i]),cimag(md->n[i]));
  }
}

void read_mesh_data(char *msh_fn,DOMD *md)
{
    void malloc_node(BOUD *bd);
  void malloc_elem(BOUD *bd);
  
  FILE *fp;
  char buf[256]="";
  double td;
  int ti,i,j,ti2;

  if((fp=fopen(msh_fn,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",msh_fn);    exit(1);  }
  strcpy(md->msh_fn,msh_fn);
  if(fgets(buf,256,fp)==NULL){
    printf("d2qte_setup.c, read_mesh_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&td)!=1){
    printf("d2qte_setup.c, read_mesh_data(), failed to read the file version. exit...\n");
    exit(1);
  }
  // check file version
  if(td<MSHVER){
    printf("This program supports mesh file version %g later. Reading file version is %g. Exit...\n",MSHVER,td);
    fclose(fp);
    exit(1);
  }
  if(fscanf(fp,"%d",&ti)!=1){
    printf("d2qte_setup.c, read_mesh_data(), failed to read the data format. exit...\n");
    exit(1);
  }
  //check data format
  if(ti!=MSHASCI){
    printf("This program supports 'ASCII' data format mesh file. Exit...\n");
    fclose(fp);
    exit(1);
  }
  if(fscanf(fp,"%d\n",&ti)!=1){
    printf("d2qte_setup.c, read_mesh_data(), failed to read the precision. exit...\n");
    exit(1);
  }
  //check data precision
  if(ti!=MSHPREC){
    printf("This program supports double precision mesh data. Exit...\n");
    fclose(fp);
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("d2qte_setup.c, read_mesh_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("d2qte_setup.c, read_mesh_data(), failed to read the line. exit...\n");
    exit(1);
  }

  if(fscanf(fp,"%d",&ti)!=1){
    printf("d2qte_setup.c, read_mesh_data(), failed to read the Nn. exit...\n");
    exit(1);
  }
  md->bd.Nn=ti;
  malloc_node(&(md->bd));
  for(i=1;i<=md->bd.Nn;i++){
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the id. exit...\n");
      exit(1);
    }
    if(ti!=i) printf("bad id %d\n",ti);
    if(fscanf(fp,"%lf",&td)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the node x. exit...\n");
      exit(1);
    }
    md->bd.xn[i]=td;
    if(fscanf(fp,"%lf",&td)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the node y. exit...\n");
      exit(1);
    }
    md->bd.yn[i]=td;
    if(fscanf(fp,"%lf\n",&td)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the node z. exit...\n");
      exit(1);
    }
  }
  if(fgets(buf,256,fp)==NULL){
    printf("d2qte_setup.c, read_mesh_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("d2qte_setup.c, read_mesh_data(), failed to read the line. exit...\n");
    exit(1);
  }

  if(fscanf(fp,"%d",&ti)!=1){
    printf("d2qte_setup.c, read_mesh_data(), failed to read the Ne. exit...\n");
    exit(1);
  }
  md->bd.Ne=ti/2;
  malloc_elem(&(md->bd));
  for(i=1;i<=md->bd.Ne;i++){
    // main
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the main id. exit...\n");
      exit(1);
    }
    if(ti!=i*2-1)    printf("bad id %d\n",ti);
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the element type. exit...\n");
      exit(1);
    }
    if(ti!=ELEMTYPE) printf("bad element type. element type must be %d\n",ELEMTYPE);
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the ti. exit...\n");
      exit(1);
    }
    for(j=0;j<ti;j++){
      if(fscanf(fp,"%d",&ti2)!=1){
        printf("d2qte_setup.c, read_mesh_data(), failed to read the ti2. exit...\n");
        exit(1);
      }
      if(j==0){
        if(ti2==OPENDID) ti2=0;
        if(md->MN>=ti2) md->bd.md[i]=ti2;
        else {
          printf("domain id %d is not defined for medium data. check domain and medium data. exit..\n",ti2);
          exit(1);
        }
      }
      else if(j==1){
        md->bd.egd[i]=ti2;
      }
    }
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the ed0. exit...\n");
      exit(1);
    }
    md->bd.ed[i][0]=ti;
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the ed1. exit...\n");
      exit(1);
    }
    md->bd.ed[i][1]=ti;
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the line. exit...\n");
      exit(1);
    }
    md->bd.ed[i][2]=ti;
    // shared
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the id. exit...\n");
      exit(1);
    }
    if(ti!=i*2)    printf("bad id :%d\n",ti);
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the element type. exit...\n");
      exit(1);
    }
    if(ti!=ELEMTYPE) printf("bad element type. element type must be %d\n",ELEMTYPE);
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the ti. exit...\n");
      exit(1);
    }
    for(j=0;j<ti;j++){
      if(fscanf(fp,"%d",&ti2)!=1){
        printf("d2qte_setup.c, read_mesh_data(), failed to read the ti2. exit...\n");
        exit(1);
      }
      if(j==0){
        if(ti2==OPENDID) ti2=0;
        if(md->MN>=ti2) md->bd.sd[i]=ti2;
        else {
          printf("domain id %d is not defined for medium data. check domain and medium data. exit..\n",ti2);
          exit(1);
        }
      }
    }
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the ed0. exit...\n");
      exit(1);
    }
    if(md->bd.ed[i][1]!=ti) printf("node id miss matched. check element id %d\n",i*2);
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the ed1. exit...\n");
      exit(1);
    }
    if(md->bd.ed[i][0]!=ti) printf("node id miss matched. check element id %d\n",i*2);
    if(fscanf(fp,"%d",&ti)!=1){
      printf("d2qte_setup.c, read_mesh_data(), failed to read the ed2. exit...\n");
      exit(1);
    }
    if(md->bd.ed[i][2]!=ti) printf("node id miss matched. check element id %d\n",i*2);
    
    // exchange open region to main domain
    if(md->bd.sd[i]==0){ // open region is sub domain
      md->bd.sd[i]=md->bd.md[i];
      md->bd.md[i]=0;
      ti=md->bd.ed[i][0];
      md->bd.ed[i][0]=md->bd.ed[i][1];
      md->bd.ed[i][1]=ti;
    }
  }
  fclose(fp);
}

void print_mesh_data(DOMD *md)
{
  printf("-- mesh data --\n");
  printf("mesh data file name        : %s\n",md->msh_fn);
  printf("numnber of nodes           : %8d\n",md->bd.Nn);
  printf("number of defined elements : %8d\n",md->bd.Ne*2);
}

void rotation_translation_obj(DOMD *md)
{
  double ct,st,r[2],M[4];
  int s;

  // rotation matrix
  st=sin(md->ra);
  ct=cos(md->ra);
  M[0]= ct;  M[1]=-st;
  M[2]= st;  M[3]= ct;;

  for(s=1;s<=md->bd.Nn;s++){
    r[0]=M[2*0+0]*md->bd.xn[s]+M[2*0+1]*md->bd.yn[s]+md->tx;
    r[1]=M[2*1+0]*md->bd.xn[s]+M[2*1+1]*md->bd.yn[s]+md->ty;
    md->bd.xn[s]=r[0];
    md->bd.yn[s]=r[1];
  }
}

void read_periodicity_data(char *prd_fn,DOMD *md)
{
  FILE *fp;
  double td;
  int ti;
  char buf[256]="";

  if((fp=fopen(prd_fn,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",prd_fn);    exit(1);  }
  strcpy(md->prd_fn,prd_fn);

  if(fgets(buf,256,fp)==NULL){
    printf("d2qte_setup.c, read_periodicity_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("d2qte_setup.c, read_periodicity_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%d %lf",&ti,&td)!=2){
    printf("d2qte_setup.c, read_periodicity_data(), failed to read the PD and d. exit...\n");
    exit(1);
  } 
  md->PD=ti;
  md->qd.d=td;
  
  fclose(fp);
}

void print_periodicity_data(DOMD *md)
{
  printf("-- periodicity data --\n");
  printf("periodicity data file  : %s\n",md->prd_fn);
  printf("periodicity setting    : %d ",md->PD);
  if(md->PD==0) printf("( no periodicity )\n");
  else if(md->PD==1){
    printf("( quasi-periodic boundary condition )\n");
    printf("lattice constant       : %8g\n",md->qd.d);
  }
  else {
    printf("periodicity setting error! periodicity setting %d is not available. Exit...\n",md->PD);
    exit(1);
  }
}

void print_periodicity_data_MKSA(DOMD *md)
{
  printf("-- periodicity data --\n");
  printf("periodicity data file  : %s\n",md->prd_fn);
  printf("periodicity setting    : %d ",md->PD);
  if(md->PD==0) printf("( no periodicity )\n");
  else if(md->PD==1){
    printf("( quasi-periodic boundary condition )\n");
    printf("lattice constant    [m]: %8g\n",OSUtoMKSA_length(md->qd.d));
  }
  else {
    printf("periodicity setting error! periodicity setting %d is not available. Exit...\n",md->PD);
    exit(1);
  }
}

void malloc_node(BOUD *bd)
{
  int N=bd->Nn;
  bd->xn=(double *)malloc(sizeof(double)*(N+1));
  bd->yn=(double *)malloc(sizeof(double)*(N+1));
}

void mfree_node(BOUD *bd)
{
  free(bd->xn);  free(bd->yn);
  bd->Nn=0;
}

void malloc_elem(BOUD *bd)
{
  int i,Ne=bd->Ne;
  bd->ed=(int **)m_alloc2(Ne+1,sizeof(int *),"malloc_elem(),bd->ed");
  bd-> x=(double **)m_alloc2(Ne+1,sizeof(double *),"malloc_elem(),bd->x");
  bd-> y=(double **)m_alloc2(Ne+1,sizeof(double *),"malloc_elem(),bd->y");
  bd->dx=(double **)m_alloc2(Ne+1,sizeof(double *),"malloc_elem(),bd->dx");
  bd->dy=(double **)m_alloc2(Ne+1,sizeof(double *),"malloc_elem(),bd->dy");
  bd->   ui=(double complex **)m_alloc2(Ne+1,sizeof(double complex *),"malloc_elem(),bd->ui");
  bd->duidn=(double complex **)m_alloc2(Ne+1,sizeof(double complex *),"malloc_elem(),bd->duidn");
  for(i=0;i<=Ne;i++){
    bd->ed[i]=(int *)m_alloc2(3,sizeof(int ),"malloc_elem(),bd->ed[i]");
    bd-> x[i]=(double *)m_alloc2(7,sizeof(double),"malloc_elem(),bd->x[i]");
    bd-> y[i]=(double *)m_alloc2(7,sizeof(double),"malloc_elem(),bd->y[i]");
    bd->dx[i]=(double *)m_alloc2(7,sizeof(double),"malloc_elem(),bd->dx[i]");
    bd->dy[i]=(double *)m_alloc2(7,sizeof(double),"malloc_elem(),bd->dy[i]");
    bd->   ui[i]=(double complex *)m_alloc2(3,sizeof(double complex),"malloc_elem(),bd->ui[i]");
    bd->duidn[i]=(double complex *)m_alloc2(3,sizeof(double complex),"malloc_elem(),bd->duidn[i]");
  }

  bd->md=(int *)m_alloc2(Ne+1,sizeof(int),"malloc_elem(),bd->md");
  bd->sd=(int *)m_alloc2(Ne+1,sizeof(int),"malloc_elem(),bd->sd");
  bd->egd=(int *)m_alloc2(Ne+1,sizeof(int),"malloc_elem(),bd->egd");
  
  // gauss nodes
  bd->xh=(double *)m_alloc2(GLH,sizeof(double),"malloc_elem(),bd->xh");
  bd->wh=(double *)m_alloc2(GLH,sizeof(double),"malloc_elem(),bd->wh");
}

void mfree_elem(BOUD *bd)
{
  int i,Ne=bd->Ne;
  for(i=0;i<=Ne;i++){
    free(bd->ed[i]);
    free(bd->x[i]);     free(bd->y[i]);
    free(bd->dx[i]);    free(bd->dy[i]);
    free(bd->ui[i]);    free(bd->duidn[i]);
  }
  free(bd->ed);
  free(bd->x);   free(bd->y);
  free(bd->dx);  free(bd->dy);
  free(bd->ui);  free(bd->duidn);

  free(bd->md);
  free(bd->sd);
  free(bd->egd);
  bd->Ne=0;
  
  free(bd->xh);
  free(bd->wh);
}

void malloc_sub_domain(DOMD *md)
{
  int *Nc,i,j;
  Nc=(int *)m_alloc2(md->MN+1,sizeof(int),"malloc_sub_domain(),Nc");
  for(i=0;i<=md->MN;i++) Nc[i]=0;

  for(i=1;i<=md->bd.Ne;i++){
    Nc[md->bd.md[i]]++;
    Nc[md->bd.sd[i]]++;
  }

  md->bd.sb=(SUBD *)m_alloc2(md->MN+2,sizeof(SUBD),"malloc_sub_domain(),md->bd.sb");
  for(i=0;i<=md->MN;i++){
    md->bd.sb[i].Ne=Nc[i];
    md->bd.sb[i].eid=(int *)m_alloc2(Nc[i]+1,sizeof(int),"malloc_sub_domain(),md->bd.sb[i],eid");
    md->bd.sb[i].u   =(double complex **)m_alloc2(Nc[i]+1,sizeof(double complex *),"malloc_sub_domain(),md->bd.sb[i],u");
    md->bd.sb[i].v   =(double complex **)m_alloc2(Nc[i]+1,sizeof(double complex *),"malloc_sub_domain(),md->bd.sb[i].v");
    md->bd.sb[i].w   =(double complex **)m_alloc2(Nc[i]+1,sizeof(double complex *),"malloc_sub_domain(),md->bd.sb[i].w");
    md->bd.sb[i].dudn=(double complex **)m_alloc2(Nc[i]+1,sizeof(double complex *),"malloc_sub_domain(),md->bd.sb[i].dudn");
    md->bd.sb[i].dvdn=(double complex **)m_alloc2(Nc[i]+1,sizeof(double complex *),"malloc_sub_domain(),md->bd.sb[i].dvdn");
    md->bd.sb[i].dwdn=(double complex **)m_alloc2(Nc[i]+1,sizeof(double complex *),"malloc_sub_domain(),md->bd.sb[i].dwdn");
    for(j=0;j<=Nc[i];j++){
      md->bd.sb[i].u[j]=(double complex *)m_alloc2(3,sizeof(double complex),"md->bd.sb[i].u[j]");
      md->bd.sb[i].v[j]=(double complex *)m_alloc2(3,sizeof(double complex),"md->bd.sb[i].v[j]");
      md->bd.sb[i].w[j]=(double complex *)m_alloc2(3,sizeof(double complex),"md->bd.sb[i].w[j]");
      md->bd.sb[i].dudn[j]=(double complex *)m_alloc2(3,sizeof(double complex),"md->bd.sb[i].dudn[j]");
      md->bd.sb[i].dvdn[j]=(double complex *)m_alloc2(3,sizeof(double complex),"md->bd.sb[i].dvdn[j]");
      md->bd.sb[i].dwdn[j]=(double complex *)m_alloc2(3,sizeof(double complex),"md->bd.sb[i].dwdn[j]");
    }
  }

  free(Nc);
}

void mfree_sub_domain(DOMD *md)
{
  int i,j;

  for(i=0;i<md->MN+1;i++){
    for(j=0;j<=md->bd.sb[i].Ne;j++){
      free(md->bd.sb[i].u[j]);      free(md->bd.sb[i].dudn[j]);
      free(md->bd.sb[i].v[j]);      free(md->bd.sb[i].dvdn[j]);
      free(md->bd.sb[i].w[j]);      free(md->bd.sb[i].dwdn[j]);
    }
    free(md->bd.sb[i].eid);
    free(md->bd.sb[i].u);    free(md->bd.sb[i].dudn);
    free(md->bd.sb[i].v);    free(md->bd.sb[i].dvdn);
    free(md->bd.sb[i].w);    free(md->bd.sb[i].dwdn);
  }

  free(md->bd.sb);
}

void init_elem_const(BOUD *bd)
{
  int i,j;

  bd->et[0]=-A_GL;  bd->et[1]= A_GL;  bd->et[2]= 0.0;
  bd->et[3]=-B_GK;  bd->et[4]= B_GK;
  bd->et[5]=-G_GK;  bd->et[6]= G_GK;
  bd->wg[0]=WA_GL;  bd->wg[1]=WA_GL;  bd->wg[2]=W0_GL;
  bd->wg[3]=0.0;    bd->wg[4]=0.0;
  bd->wg[5]=0.0;    bd->wg[6]=0.0;
  bd->wk[0]=WA_GK;  bd->wk[1]=WA_GK;  bd->wk[2]=W0_GK;
  bd->wk[3]=WB_GK;  bd->wk[4]=WB_GK;
  bd->wk[5]=WG_GK;  bd->wk[6]=WG_GK;
  for(i=0;i<3;i++)    for(j=0;j<7;j++) bd->M[i][j]=Mn(i,bd->et[j]);
  
  // gauleg 
  gauleg(-1.0,1.0,bd->xh,bd->wh,GLH);
}

void init_sub_domain(DOMD *md)
{
  int d,i,c;

  for(d=0;d<=md->MN;d++){
    c=1;
    for(i=1;i<=md->bd.Ne;i++){
      if(md->bd.md[i]==d){
        md->bd.sb[d].eid[c]=i;
        c++;
      }
      else if(md->bd.sd[i]==d){
        md->bd.sb[d].eid[c]=-i;
        c++;
      }
    }
  }
}

void init_boundary_data(DOMD *md)
{
  double complex ui,gui[2];
  double rn[2][3],r[2],dr[2],i_J;
  int i,j,id;

  for(i=1;i<=md->bd.Ne;i++){
    for(j=0;j<3;j++){
      id=md->bd.ed[i][j];
      rn[0][j]=md->bd.xn[id];
      rn[1][j]=md->bd.yn[id];
    }
    for(j=0;j<7;j++){
      rn_eta ( r,rn,md->bd.et[j]);
      drn_eta(dr,rn,md->bd.et[j]);
      md->bd.x [i][j]= r[0];    md->bd.y [i][j]= r[1];
      md->bd.dx[i][j]=dr[0];    md->bd.dy[i][j]=dr[1];
    }
  }
  // incident field
  for(i=0;i<=md->bd.Ne;i++){
    for(j=0;j<3;j++){
      md->bd.ui[i][j]=0.0;
      md->bd.duidn[i][j]=0.0;
    }
  }
  for(i=1;i<=md->bd.Ne;i++){
    if(md->bd.md[i]==0){
      for(j=0;j<3;j++){
        r [0]=md->bd.x [i][j];     r [1]=md->bd.y [i][j];
        dr[0]=md->bd.dx[i][j];     dr[1]=md->bd.dy[i][j];
        i_J=1.0/sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
        infd_grad(&ui,gui,r,&(md->wd));
        md->bd.ui[i][j]=ui;
        md->bd.duidn[i][j]=i_J*(gui[0]*dr[1]-gui[1]*dr[0]);
      }
    }
  }
}

int domain_crt(int id,double *rt,DOMD *md)
{
  double rs0[2],rs1[2],rs2[2],S02,C02,S21,C21,angle;
  int s,eid,aeid;

  angle=0.0;
  for(s=1;s<=md->bd.sb[id].Ne;s++){
    eid=md->bd.sb[id].eid[s];
    aeid=abs(eid);
    rs0[0]=md->bd.xn[md->bd.ed[aeid][0]]-rt[0];    rs0[1]=md->bd.yn[md->bd.ed[aeid][0]]-rt[1];
    rs1[0]=md->bd.xn[md->bd.ed[aeid][1]]-rt[0];    rs1[1]=md->bd.yn[md->bd.ed[aeid][1]]-rt[1];
    rs2[0]=md->bd.xn[md->bd.ed[aeid][2]]-rt[0];    rs2[1]=md->bd.yn[md->bd.ed[aeid][2]]-rt[1];

    C02=rs0[0]*rs2[0]+rs0[1]*rs2[1];
    S02=rs0[0]*rs2[1]-rs2[0]*rs0[1];

    C21=rs2[0]*rs1[0]+rs2[1]*rs1[1];
    S21=rs2[0]*rs1[1]-rs1[0]*rs2[1];

    if(eid>0)  angle+=atan2(S02,C02)+atan2(S21,C21);
    else       angle-=atan2(S02,C02)+atan2(S21,C21);
  }
  if(fabs(angle)<MEPS*20.0) return 0;  // external
  else return 1; // internal include boundary
}
