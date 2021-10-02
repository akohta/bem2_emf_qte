#include "bem2_emf_qte.h"

void force_FN(double *Fr,double *Nz,double *rc,int type,DOMD *md)
{
  void force_FN_GL(double *Fr,double *Nz,double *rc,DOMD *md);
  void force_FN_GK(double *Fr,double *Nz,double *rc,DOMD *md);
  
  if(type==0) force_FN_GL(Fr,Nz,rc,md);
  else        force_FN_GK(Fr,Nz,rc,md);
}

///////////////////////////////////////////////////////////////////////
void force_FN_GL(double *Fr,double *Nz,double *rc,DOMD *md)
{
  double complex Ez,Hr[2];
  double r[2],dx,dy,Txx,Txy,Tyx,Tyy,aez2,ahx2,ahy2,dFx,dFy,epsd,tE,tH;
  int s,i,eid;

  for(i=0;i<2;i++) Fr[i]=0.0;
  *Nz=0.0;

  epsd=md->n[0]*md->n[0];
  for(s=1;s<=md->bd.sb[0].Ne;s++){
    eid=md->bd.sb[0].eid[s];
    for(i=0;i<3;i++){
      r[0]=md->bd.x [eid][i];       r[1]=md->bd.y [eid][i];
      dx  =md->bd.dx[eid][i];       dy  =md->bd.dy[eid][i];
      infd_EH(&Ez,Hr,r,&(md->wd));
      Ez   +=md->bd.sb[0].u[s][i];
      Hr[0]+=md->bd.sb[0].v[s][i];
      Hr[1]+=md->bd.sb[0].w[s][i];
      aez2=creal(Ez*conj(Ez));
      ahx2=creal(Hr[0]*conj(Hr[0]));
      ahy2=creal(Hr[1]*conj(Hr[1]));
      tE=0.25*epsd*aez2;
      tH=0.25*(ahx2-ahy2);
      Txx=-tE+tH;
      Txy=0.50*creal(Hr[0]*conj(Hr[1]));
      Tyx=Txy;
      Tyy=-tE-tH;
      dFx=Txx*dy-Txy*dx;
      dFy=Tyx*dy-Tyy*dx;
      Fr[0]+=md->bd.wg[i]*dFx;
      Fr[1]+=md->bd.wg[i]*dFy;
      *Nz+=md->bd.wg[i]*((r[0]-rc[0])*dFy-(r[1]-rc[1])*dFx);
    }
  }

  Fr[0]*=-1.0;
  Fr[1]*=-1.0;
  *Nz*=-1.0;
}

void force_FN_GK(double *Fr,double *Nz,double *rc,DOMD *md)
{
  double complex Ez,Hr[2],Ezs,Hrs[2];
  double r[2],dx,dy,Txx,Txy,Tyx,Tyy,aez2,ahx2,ahy2,dFx,dFy,epsd,tE,tH,fx,fy,nz;
  int s,i,eid;

  fx=0.0;
  fy=0.0;
  nz=0.0;

  epsd=md->n[0]*md->n[0];
  #pragma omp parallel for schedule(dynamic) reduction(+:fx,fy,nz) private(eid,i,r,dx,dy,Ez,Hr,Ezs,Hrs,aez2,ahx2,ahy2,tE,tH,Txx,Txy,Tyx,Tyy,dFx,dFy)
  for(s=1;s<=md->bd.sb[0].Ne;s++){
    eid=md->bd.sb[0].eid[s];
    for(i=0;i<7;i++){
      r[0]=md->bd.x [eid][i];       r[1]=md->bd.y [eid][i];
      dx  =md->bd.dx[eid][i];       dy  =md->bd.dy[eid][i];
      infd_EH(&Ez,Hr,r,&(md->wd));
      if(i<3){
        Ez   +=md->bd.sb[0].u[s][i];
        Hr[0]+=md->bd.sb[0].v[s][i];
        Hr[1]+=md->bd.sb[0].w[s][i];
      }
      else {
        EH_bv(&Ezs,Hrs,0,md->bd.et[i],s,md);
        Ez   +=Ezs;
        Hr[0]+=Hrs[0];
        Hr[1]+=Hrs[1];
      }
      aez2=creal(Ez*conj(Ez));
      ahx2=creal(Hr[0]*conj(Hr[0]));
      ahy2=creal(Hr[1]*conj(Hr[1]));
      tE=0.25*epsd*aez2;
      tH=0.25*(ahx2-ahy2);
      Txx=-tE+tH;
      Txy=0.50*creal(Hr[0]*conj(Hr[1]));
      Tyx=Txy;
      Tyy=-tE-tH;
      dFx=Txx*dy-Txy*dx;
      dFy=Tyx*dy-Tyy*dx;
      fx+=md->bd.wk[i]*dFx;
      fy+=md->bd.wk[i]*dFy;
      nz+=md->bd.wk[i]*((r[0]-rc[0])*dFy-(r[1]-rc[1])*dFx);
    }
  }

  Fr[0]=-1.0*fx;
  Fr[1]=-1.0*fy;
  *Nz=-1.0*nz;
}
