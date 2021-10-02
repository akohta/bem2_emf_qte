#include "pb_elem.h"

void coefficient_qpgf_GL(double _Complex *CC,double *rt,int s,DOMD *md)
{
  double _Complex qgf,dqgf[2];
  double r[2],dx,dy,aJ,i_R,tD;
  int i,erc;

  CC[6]=0.0;
  for(i=0;i<3;i++){
    r[0]=md->bd.x[s][i]-rt[0];
    r[1]=md->bd.y[s][i]-rt[1];
    dx=md->bd.dx[s][i];
    dy=md->bd.dy[s][i];
    aJ=sqrt(dx*dx+dy*dy);
    erc=d2hm_qpgf_d1_ew(&qgf,dqgf,r,MEPS,&(md->qd));
    i_R=1.0/sqrt(r[0]*r[0]+r[1]*r[1]);
    tD=( r[0]*dy-r[1]*dx)*i_R;

    if(erc>=0){
      CC[  i]=md->bd.wg[i]*qgf*aJ;
      CC[3+i]=md->bd.wg[i]*(dqgf[0]*dy-dqgf[1]*dx);
    }
    else {
      printf("d2hm_qpgf_d1_ew() error in coefficient_qpgf_GL! erc=%d. Exit...\n",erc);
      exit(1);
    }
    CC[6]+=md->bd.wg[i]*tD*i_R;
  }
  CC[6]*=1.0/(2.0*M_PI);
}

void coefficient_qpgf_GK(double complex *CC,double *rt,int s,DOMD *md)
{
  double _Complex qgf,dqgf[2],tc0,tc1;
  double dx,dy,r[2],i_R,aJ,tD;
  int i,j,erc;

  for(i=0;i<7;i++) CC[i]=0.0;
  for(i=0;i<7;i++){
    r[0]=md->bd.x[s][i]-rt[0];
    r[1]=md->bd.y[s][i]-rt[1];
    dx=md->bd.dx[s][i];
    dy=md->bd.dy[s][i];
    aJ=sqrt(dx*dx+dy*dy);
    erc=d2hm_qpgf_d1_ew(&qgf,dqgf,r,MEPS,&(md->qd));
    i_R=1.0/sqrt(r[0]*r[0]+r[1]*r[1]);
    tD=( r[0]*dy-r[1]*dx)*i_R;

    if(erc>=0){
      tc0=qgf*aJ;
      tc1=dqgf[0]*dy-dqgf[1]*dx;
      for(j=0;j<3;j++){
        CC[  j]+=md->bd.wk[i]*md->bd.M[j][i]*tc0;      CC[3+j]+=md->bd.wk[i]*md->bd.M[j][i]*tc1; 
      }
    }
    else {
      printf("d2hm_qpgf_d1_ew() error in coefficient_qpgf_GK! erc=%d. Exit...\n",erc);
      exit(1);
    }
    CC[6]+=md->bd.wk[i]*tD*i_R;
  }
  CC[6]*=1.0/(2.0*M_PI);
}

void coefficient_qpgf_HP(double complex *CC,double *rt,int s,DOMD *md)
{
  double complex qgf,dqgf[2],tc0,tc1;
  double r[2],dr[2],rs[2][3],i_R,aJ,tD;
  int i,j,erc;

  for(i=0;i<3;i++){
    rs[0][i]=md->bd.x[s][i];
    rs[1][i]=md->bd.y[s][i];
  }

  for(i=0;i<7;i++) CC[i]=0.0;
  for(i=0;i<GLH;i++){
    rs_eta(r,rs,md->bd.xh[i]);
    drs_eta(dr,rs,md->bd.xh[i]);
    r[0]-=rt[0];
    r[1]-=rt[1];
    aJ=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
    erc=d2hm_qpgf_d1_ew(&qgf,dqgf,r,MEPS,&(md->qd));
    i_R=1.0/sqrt(r[0]*r[0]+r[1]*r[1]);
    tD=( r[0]*dr[1]-r[1]*dr[0])*i_R;

    if(erc>=0){
      tc0=qgf*aJ;
      tc1=dqgf[0]*dr[1]-dqgf[1]*dr[0];
      for(j=0;j<3;j++){
        CC[  j]+=md->bd.wh[i]*Mn(j,md->bd.xh[i])*tc0;      CC[3+j]+=md->bd.wh[i]*Mn(j,md->bd.xh[i])*tc1; 
      }
    }
    else {
      printf("d2hm_qpgf_d1_ew() error in coefficient_qpgf_HP()! erc=%d. Exit...\n",erc);
      exit(1);
    }
    CC[6]+=md->bd.wh[i]*tD*i_R;
  }
  CC[6]*=1.0/(2.0*M_PI);
}

void coefficient_qpgf_DE(double complex *CC,double *rt,int s,DOMD *md)
{
  double complex sqG(double eta,void *tmp);
  double complex sqH(double eta,void *tmp);
  double sF(double eta,void *tmp);

  TDATA td;
  double err0;
  int i;

  td.k=md->qd.k;
  td.rt[0]=rt[0];  td.rt[1]=rt[1];
  for(i=0;i<3;i++){
    td.rs[0][i]=md->bd.x[s][i];
    td.rs[1][i]=md->bd.y[s][i];
  }
  td.qd=&(md->qd);
  
  td.type=0;
  CC[0]=deintz(sqG,-1.0,1.0,&td,DEPS,&err0);    if(err0<0.0){ printf("DE integration error coefficient_qpgf_DE(), sqG0! Exit...\n");  exit(1);  }
  td.type=1;
  CC[1]=deintz(sqG,-1.0,1.0,&td,DEPS,&err0);    if(err0<0.0){ printf("DE integration error coefficient_qpgf_DE(), sqG1! Exit...\n");  exit(1);  }
  td.type=2;
  CC[2]=deintz(sqG,-1.0,1.0,&td,DEPS,&err0);    if(err0<0.0){ printf("DE integration error coefficient_qpgf_DE(), sqG2! Exit...\n");  exit(1);  }
  
  td.type=0;
  CC[3]=deintz(sqH,-1.0,1.0,&td,DEPS,&err0);    if(err0<0.0){ printf("DE integration error coefficient_qpgf_DE(), sqH0! Exit...\n");  exit(1);  }
  td.type=1;
  CC[4]=deintz(sqH,-1.0,1.0,&td,DEPS,&err0);    if(err0<0.0){ printf("DE integration error coefficient_qpgf_DE(), sqH1! Exit...\n");  exit(1);  }
  td.type=2;
  CC[5]=deintz(sqH,-1.0,1.0,&td,DEPS,&err0);    if(err0<0.0){ printf("DE integration error coefficient_qpgf_DE(), sqH2! Exit...\n");  exit(1);  }
  
  if( ((2.0*fabs(td.rs[0][0]-rt[0])/fabs(td.rs[0][0]+rt[0])<MEPS) && (2.0*fabs(td.rs[0][1]-rt[0])/fabs(td.rs[0][1]+rt[0])<MEPS) 
       && (2.0*fabs(td.rs[0][2]-rt[0])/fabs(td.rs[0][2]+rt[0])<MEPS)) || 
      ((2.0*fabs(td.rs[1][0]-rt[1])/fabs(td.rs[1][0]+rt[1])<MEPS) && (2.0*fabs(td.rs[1][1]-rt[1])/fabs(td.rs[1][1]+rt[1])<MEPS)
       && (2.0*fabs(td.rs[1][2]-rt[1])/fabs(td.rs[1][2]+rt[1])<MEPS))){
    CC[6]=0.0;
  }
  else{
    CC[6]=1.0/(2.0*M_PI)*deintd(sF,-1.0,1.0,&td,DEPS,&err0);     if(err0<0.0){ printf("DE integration error coefficient_qpgf_DE(),sF! Exit...\n");  exit(1);  }
  }
}

double complex sqG(double eta,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex qgf,dqgf[2];
  double r[2],rb[2],dr[2],aJ;
  int erc;

  rs_eta(rb,td->rs,eta);
  drs_eta(dr,td->rs,eta);
  aJ=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
  r[0]=rb[0]-td->rt[0];
  r[1]=rb[1]-td->rt[1];
  erc=d2hm_qpgf_d1_ew(&qgf,dqgf,r,MEPS,td->qd);
  if(erc>=0) return qgf*Mn(td->type,eta)*aJ;
  else {
    printf("d2hm_qpgf_d1_ew() error in sqG! erc=%d. Exit...\n",erc);
    exit(1);
  }
}

double complex sqH(double eta,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex qgf,dqgf[2];
  double r[2],rb[2],dr[2];
  int erc;

  rs_eta(rb,td->rs,eta);
  drs_eta(dr,td->rs,eta);
  r[0]=rb[0]-td->rt[0];
  r[1]=rb[1]-td->rt[1];
  erc=d2hm_qpgf_d1_ew(&qgf,dqgf,r,MEPS,td->qd);
  if(erc>=0) return (dqgf[0]*dr[1]-dqgf[1]*dr[0])*Mn(td->type,eta);
  else {
    printf("d2hm_qpgf_d1_ew() error in sqH! erc=%d.Exit...\n",erc);
    exit(1);
  }
}

int coefficient_qpgf_NV(double complex *CGK,double *rt,int s,DOMD *md)
{
  double complex CGL[7],qgf,dqgf[2],tc0,tc1;
  double dx,dy,r[2],i_R,aJ,tD,CD;
  int i,j,erc;

  for(i=0;i<7;i++){
    CGK[i]=0.0;
    CGL[i]=0.0;
  }
  for(i=0;i<7;i++){
    r[0]=md->bd.x[s][i]-rt[0];
    r[1]=md->bd.y[s][i]-rt[1];
    dx=md->bd.dx[s][i];
    dy=md->bd.dy[s][i];
    aJ=sqrt(dx*dx+dy*dy);
    erc=d2hm_qpgf_d1_ew(&qgf,dqgf,r,MEPS,&(md->qd));
    i_R=1.0/sqrt(r[0]*r[0]+r[1]*r[1]);
    tD=( r[0]*dy-r[1]*dx)*i_R;

    if(erc>=0){
      tc0=qgf*aJ;
      tc1=dqgf[0]*dy-dqgf[1]*dx;
      for(j=0;j<3;j++){
        CGK[  j]+=md->bd.wk[i]*md->bd.M[j][i]*tc0;      CGK[3+j]+=md->bd.wk[i]*md->bd.M[j][i]*tc1;
        CGL[  j]+=md->bd.wg[i]*md->bd.M[j][i]*tc0;      CGL[3+j]+=md->bd.wg[i]*md->bd.M[j][i]*tc1;
      }
    }
    else {
      printf("d2hm_qpgf_d1_ew() error in coefficient_qpgf_NV! erc=%d. Exit...\n",erc);
      printf("r=(%15.14e,%15.14e)\n",r[0],r[1]);
      exit(1);
    }
    CGK[6]+=md->bd.wk[i]*tD*i_R;
    CGL[6]+=md->bd.wg[i]*tD*i_R;
  }
  CGK[6]*=1.0/(2.0*M_PI);
  CGL[6]*=1.0/(2.0*M_PI);
    
  // check
  j=0;
  for(i=0;i<7;i++){
    CD=2.0*cabs(CGK[i]-CGL[i])/cabs(CGK[i]+CGL[i]);
    if(CD>IEPS) j++;
  }
  if(j!=0){
    coefficient_qpgf_HP(CGK,rt,s,md); 
    return 1;
  }
  else {
    return 0;
  }
  
  return 0;
}

void coefficient_qpgf_bd_eta(double complex *CC,double eta_t,int s,DOMD *md)
{
  double complex sqG_bd(double xi,void *tmp);
  double complex sqH_bd(double xi,void *tmp);
  double sF_bd(double xi,void *tmp);

  TDATA td;
  double err0,err1,ts,tr;
  int i;
  
  td.k=md->qd.k;
  td.eta_t=eta_t;
  for(i=0;i<3;i++){
    td.rs[0][i]=md->bd.x[s][i];
    td.rs[1][i]=md->bd.y[s][i];
  }
  td.qd=&(md->qd);

  ts=(td.rs[0][1]-td.rs[0][2])*td.rs[1][0]+(td.rs[0][2]-td.rs[0][0])*td.rs[1][1]+(td.rs[0][0]-td.rs[0][1])*td.rs[1][2];
  tr=pow(td.rs[0][1]-td.rs[0][0],2)+pow(td.rs[1][1]-td.rs[1][0],2);
  if( fabs(ts)/tr < MEPS){ // Cs=0.0;
    td.Cs=0.0;
  }
  else td.Cs=0.5*i_A_GL*i_A_GL*i_A_GL*ts;

  td.type=0;  
  CC[0]=-deintz(sqG_bd,0.0,-1.0-eta_t,&td,DEPS,&err0) 
    +    deintz(sqG_bd,0.0, 1.0-eta_t,&td,DEPS,&err1);    if(err0<0.0||err1<0.0){ printf("DE integration error coefficient_qpgf_bd_eta(), sqG0_bd! Exit...\n");  exit(1);  }
  td.type=1;
  CC[1]=-deintz(sqG_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
    +    deintz(sqG_bd,0.0, 1.0-eta_t,&td,DEPS,&err1);    if(err0<0.0||err1<0.0){ printf("DE integration error coefficient_qpgf_bd_eta(), sqG1_bd! Exit...\n");  exit(1);  }
  td.type=2;
  CC[2]=-deintz(sqG_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
    +    deintz(sqG_bd,0.0, 1.0-eta_t,&td,DEPS,&err1);    if(err0<0.0||err1<0.0){ printf("DE integration error coefficient_qpgf_bd_eta(), sqG2_bd! Exit...\n");  exit(1);  }
    
  td.type=0;
  CC[3]=-deintz(sqH_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
    +    deintz(sqH_bd,0.0, 1.0-eta_t,&td,DEPS,&err1); if(err0<0.0||err1<0.0){ printf("DE integration error coefficient_qpgf_bd_eta(), sqH0_bd! Exit...\n");    exit(1);  }
  td.type=1;
  CC[4]=-deintz(sqH_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
    +    deintz(sqH_bd,0.0, 1.0-eta_t,&td,DEPS,&err1); if(err0<0.0||err1<0.0){ printf("DE integration error coefficient_qpgf_bd_eta(), sqH1_bd! Exit...\n");    exit(1);  }
  td.type=2;
  CC[5]=-deintz(sqH_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
    +    deintz(sqH_bd,0.0, 1.0-eta_t,&td,DEPS,&err1); if(err0<0.0||err1<0.0){ printf("DE integration error coefficient_qpgf_bd_eta(), sqH2_bd! Exit...\n");    exit(1);  }
  
  if(td.Cs==0.0) CC[6]=0.0;
  else CC[6]=td.Cs/(2.0*M_PI)*(-deintd(sF_bd,0.0,-1.0-eta_t,&td,DEPS,&err0)
			       +deintd(sF_bd,0.0, 1.0-eta_t,&td,DEPS,&err1));
  if(err0<0.0||err1<0.0){ printf("DE integration error coefficient_qpgf_bd_eta(), sF_bd! Exit...\n");    exit(1);  }
}

double complex sqG_bd(double xi,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex qgf,dqgf[2];
  double r[2],dr[2],dR[2],aJ;
  int erc;

  drs_eta(dr,td->rs,xi+td->eta_t);
  aJ=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
  
  drs_eta(dR,td->rs,0.5*xi+td->eta_t);
  r[0]=xi*dR[0];
  r[1]=xi*dR[1];
  erc=d2hm_qpgf_d1_ew_cs2(&qgf,dqgf,r,MEPS,td->qd);
  if(erc>=0)    return qgf*Mn(td->type,xi+td->eta_t)*aJ;
  else {
    printf("d2hm_qpgf_d1_ew_cs2() error in sqG_bd! erc=%d. Exit...\n",erc);
    printf("r[0]=%15.14e,r[1]=%15.14e\n",r[0],r[1]);
    exit(1);
  }
}

double complex sqH_bd(double xi,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex qgf,dqgf[2];
  double r[2],dr[2],dR[2],i_aJR,tD,tE;
  int erc;

  drs_eta(dr,td->rs,xi+td->eta_t);

  drs_eta(dR,td->rs,0.5*xi+td->eta_t);
  r[0]=xi*dR[0];
  r[1]=xi*dR[1];
  erc=d2hm_qpgf_d1_ew_cs2(&qgf,dqgf,r,MEPS,td->qd);
  if(erc>=0){
    i_aJR=1.0/sqrt(dR[0]*dR[0]+dR[1]*dR[1]);
    tD=fabs(xi)*i_aJR*td->Cs;
    tE=-xi/fabs(xi)*i_aJR*(dR[0]*dr[0]+dR[1]*dr[1]);
    return (dqgf[0]*tD+dqgf[1]*tE)*Mn(td->type,xi+td->eta_t);
  }
  else {
    printf("d2hm_qpgf_d1_ew_cs2() error in sqH_bd! erc=%d Exit...\n",erc);
    exit(1);
  }
}
