// 9/17/2014: Hexc & Olesya - Checking on the code (plus cleaning up)
//
#include "BSIntegrator1.hh"
#include "BSEquation.hh"
#include "BSEquationInG4.hh"
#include "MAGCOSUnits.hh"

//constructor destructor

BSIntegrator1* BSIntegrator1::instance = 0;
////////////////////////////////////////////////////////////////////////////////
//
BSIntegrator1::BSIntegrator1()
{
  theEquation= 0;
  eps=1.0e-5;
  hmin=1.0e-50;
  hmax=1.0e50;
  kmax=-1000;
  kount=0;
  maxstep=1000;
  
  /*xp_p=0;
    yp_p=0;*/
  /* x_p=0;
     d_p=0;*/
  
  
  xmax=1.0e50;
  
  verbose =0;
  
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
BSIntegrator1::BSIntegrator1(BSEquation* anEquation)
{
  eps=1.0e-5;
  hmin=1.0e-50;
  hmax=1.0e50;
  kmax=1000;
  kount=0;
  maxstep=1000;
  theEquation=anEquation;
  
  /*xp_p=0;
    yp_p=0;*/
  
  /*x_p=0;
    d_p=0;*/
  
  
  xmax=1.0e50;
  
  verbose =1;
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
BSIntegrator1::~BSIntegrator1()
{
  if (theEquation) delete theEquation;
  /*if (xp_p) delete xp_p;
    if (yp_p) delete yp_p;*/
  
}

////////////////////////////////////////////////////////////////////////////////
//
BSIntegrator1* BSIntegrator1::getInstance()
{
  if (instance == 0) instance = new BSIntegrator1;
  return instance;
} 

////////////////////////////////////////////////////////////////////////////////
//
//integration
bool BSIntegrator1::do_integration(G4double ystart[], const G4double x1, const G4double h1, 
				   G4int &nok, G4int &nbad)
{ 
  //check if equation to be integrated exist
  //-----------------------
  if  (! theEquation){
    std::cout<<" No equation to be integrated "<<std::endl;
    return false;
  }
  
  // integration 
  const G4double TINY=1.0e-30;
  G4int i,nstp;
  G4double xsav,x,hnext,hdid,h;
  xsav=0.;
  G4int nvar= theEquation->GetNvar();
  G4double yscal[10];
  G4double y[10];
  G4double dydx[10];
  
  /*if (xp_p) delete xp_p;
    if (yp_p) 
    { for (G4int i=0; i<last_nvar ; i++) delete yp_p[i];
    delete yp_p;
    }*/
  xp_p.clear();
  for (unsigned int  i=0; i<yp_p.size() ; i++) yp_p[i].clear(); 
  yp_p.clear();   
  
  last_nvar=nvar;
  
  if (kmax > 0){
    /*xp_p=new G4double[kmax];
      yp_p=new G4double*[nvar];
      for (G4int i=0; i<nvar; i++)  yp_p[i]= new G4double[kmax];*/
    for (G4int i=0; i<nvar; i++)  yp_p.push_back(std::vector<double>());
  }
  
  
  x=x1;
  h=std::min(h1,xmax-x);
  nbad = kount = nok=0;
  
  for (i=0;i<nvar;i++) y[i]=ystart[i];  
  
  // register last positions in  equation	
  theEquation->SetLastPosition(G4ThreeVector(y[0],y[1],y[2]));
  theEquation->SetLastX(x);
  
  
  if (kmax > 0) xsav=x-dxsav*2.0;
  for (nstp=0; nstp<maxstep;nstp++){ 
    y[6]=x;
    theEquation->derivs(x,y,dydx);
    h=std::min(h,xmax-x);
    h=std::min(h,hmax);
    if (verbose >1) verbose = 1;
    if (verbose >0 && nstp >maxstep-5) {
      std::cout <<"h "<<h/Re<<std::endl;
      std::cout <<"xlast "<<xlast/Re<<" x "<<x/re<<std::endl; 
      G4ThreeVector pos=G4ThreeVector(y[0],y[1],y[2])/Re;	
      std::cout <<"pos "<<pos<<'\t'<<pos.mag()<<std::endl;
      verbose = 2;
      
    }
    for (i=0;i<nvar;i++){
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
      if (dydx[i]==0 && y[i]==0)  yscal[i]=h;
    }	 
    
    if (kmax > 0 && kount < kmax-2 && fabs(x-xsav) > fabs(dxsav)){ 
      for (i=0;i<nvar;i++) yp_p[i].push_back(y[i]);
      xp_p.push_back(x);
      kount++;
      xsav=x;
    }
    //double re=6371.2*1000.*1000.;
    std::cout.precision(20); 	    
    //std::cout <<"xlast "<<xlast/re<<" x "<<x/re<<endl; 	    
    if (theEquation->InterruptCondition(x, y))
      {for (i=0;i<nvar;i++) yp_p[i].push_back(y[i]);
	xp_p.push_back(x);
	kount++;
	std::cout << "Interuption condition has been fullfilled"<<std::endl;
	return false;}
    
    xlast=x;
    for (i=0;i<nvar;i++) ylast[i]=y[i];
    if (IsTheEquationOfMotionForG4Tracking){
      dynamic_cast <BSEquationInG4*> (theEquation)
	->SetDidCrossBoundaryDuringBsstep(false);
    }			  
    bsstep(y,dydx,x,h,yscal,hdid,hnext);
    if (verbose >1){
      std::cout<<"After bsstep"<<std::endl;
      G4ThreeVector pos = G4ThreeVector(y[0],y[1],y[2])/Re;
      std::cout<<"Pos "<<pos<<'\t'<<pos.mag()<<std::endl;
      std::cout<<"x "<< x/Re<<std::endl;
      std::cout<<"hdid "<< hdid/Re<<std::endl;
      std::cout<<"hnext "<< hnext/Re<<std::endl;
      std::cout<<"y[6] "<< y[6]/Re<<std::endl;
      
    }
    if (IsTheEquationOfMotionForG4Tracking){
      if (!(dynamic_cast<BSEquationInG4*> (theEquation)
	    ->GetDidCrossBoundaryDuringBsstep())){
	y[6]=x;
	if (verbose >1) std::cout<<"y[6] set to x "<< y[6]/Re<<std::endl;
      }
      else if ((x-y[6])/x <eps*2.){    
	y[6]=x;
	dynamic_cast <BSEquationInG4*> (theEquation)
	  ->SetDidCrossBoundaryDuringBsstep(false);
	if (verbose >1) std::cout<<"y[6] set to x "<< y[6]/Re<<std::endl;	
      }
    }  
    
    //std::cout<<x<<std::endl;
    int ncond =theEquation->StopCondition(x, y);
    //std::cout<<"ncond "<<ncond<<std::endl;
    if  (ncond > 0){ 
      if (hdid == h) ++nok; else ++nbad;
      if (ncond == 1){
	for (i=0;i<nvar;i++) ystart[i]=y[i];
	if (kmax > 0) {
	  for (i=0;i<nvar;i++) yp_p[i].push_back(y[i]);
	  xp_p.push_back(x);
	  kount++;
	}
	//G4cout<<"nstep "<<nstp<<std::endl;
	
	return true;
      }
      h=hnext;
      // new y recompute ylast in equation
      theEquation->SetLastPosition(G4ThreeVector(y[0],y[1],y[2]));
      theEquation->SetLastX(x);
      y[6]=xlast;
    }
    
    else {
      /*G4cout<<"nstep for step_before "<<nstp<<std::endl;
	G4cout<<"xlast "<<xlast/6371.2/km<<std::endl;
	G4cout<<"ylast[6] "<<ylast[6]/6371.2/km<<std::endl;
	G4cout<<"ylast[0] "<<ylast[0]/6371.2/km<<std::endl;
	G4cout<<"ylast[1] "<<ylast[1]/6371.2/km<<std::endl;
	G4cout<<"ylast[2] "<<ylast[2]/6371.2/km<<std::endl;
	G4cout<<"x "<<x/6371.2/km<<std::endl;
	G4cout<<"y[6] "<<y[6]/6371.2/km<<std::endl;
	G4cout<<"y[0] "<<y[0]/6371.2/km<<std::endl;
	G4cout<<"y[1] "<<y[1]/6371.2/km<<std::endl;
	G4cout<<"y[2] "<<y[2]/6371.2/km<<std::endl;
	G4cout<<"h "<<h/6371.2/km<<std::endl;*/
      h= theEquation->StepBeforeStop(xlast, ylast,x,y);
      x=xlast;
      for (i=0;i<nvar;i++) y[i]=ylast[i];
      y[6]=xlast;
    } 
    
  }
  std::cout << "Too many steps in integration"<<std::endl;
  return false;
  
}
////////////////////////////////////////////////////////////////////////////////
//

void BSIntegrator1::bsstep(G4double y[], G4double dydx[], G4double &xx, 
			   const G4double htry, const G4double yscal[], 
			   G4double &hdid, G4double &hnext)
{
  const G4int KMAXX=8, IMAXX=(KMAXX+1);
  const G4double SAFE1=0.25, SAFE2=0.7, REDMAX=1.0e-5, REDMIN=0.7;
  const G4double TINY=1.0e-30, SCALMX=0.1;
  //	static const G4int nseq_d[IMAXX]={2,4,6,8,10,12,14,16,18};
  static G4int first=1,kmax,kopt;
  static G4double epsold = -1.0,xnew;
  static G4double a[IMAXX];
  static G4double alf[KMAXX][KMAXX];
  bool exitflag=false;
  G4int i,iq,k,kk,kkm,reduct;
  // added by Laurent to initialise kkm
  kkm =0;
  G4double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
  // added by Laurent to initialise scale
  scale =1.;	
  G4int nseq[IMAXX]={2,4,6,8,10,12,14,16,18};
  G4double err[KMAXX];
  
  int nv=theEquation->GetNvar();
  if (IsTheEquationOfMotionForG4Tracking){
    dynamic_cast <BSEquationInG4*> (theEquation)
      ->SetDidCrossBoundaryDuringMmid(false);
  }			  
  G4double yerr[10],ysav[10],yseq[10];
  /*x_p=new G4double[KMAXX];
    d_p=new G4double*[nv];*/
  if (x_p.size()<=0) {
    x_p.insert(x_p.end(),KMAXX,0.);
    d_p.insert(d_p.end(),nv, std::vector< double >());
    for (unsigned i =0;i<unsigned(nv);i++)
      d_p[i].insert(d_p[i].end(),KMAXX,0.);
  }
  
  
  //for (i=0; i<nv; i++) d_p[i] = new G4double[KMAXX];
  if (eps != epsold) {
    hnext = xnew = -1.0e29;
    eps1=SAFE1*eps;
    a[0]=nseq[0]+1;
    for (k=0;k<KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
    for (iq=1;iq<KMAXX;iq++) {
      for (k=0;k<iq;k++)
	alf[k][iq]=pow(eps1,(a[k+1]-a[iq+1])/
		       ((a[iq+1]-a[0]+1.0)*(2*k+3)));
    }
    epsold=eps;
    for (kopt=1;kopt<KMAXX-1;kopt++)
      if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
    kmax=kopt;
  }
  h=htry;
  for (i=0;i<nv;i++) ysav[i]=y[i];
  if (xx != xnew || h != hnext) {
    first=1;
    kopt=kmax;
  }
  reduct=0;
  for (;;) {
    for (k=0;k<=kmax;k++) {
      xnew=xx+h;
      if (xnew == xx) 
	{G4cout<<"first integrator step size underflow in bsstep"<<std::endl;
	}
      mmid(ysav,dydx,xx,h,nseq[k],yseq);
      //xest=std::std::sqr(h/nseq[k]);
      xest=(h/nseq[k])*(h/nseq[k]);
      pzextr(k,xest,yseq,y,yerr);
      if (k != 0) {
	errmax=TINY;
	for (i=0;i<nv;i++) errmax=std::max(errmax,fabs(yerr[i]/yscal[i]));
	errmax /= eps;
	kkm=k-1;
	err[kkm]=std::pow(errmax/SAFE1,1.0/(2*kkm+3));
      }
      if (k != 0 && (k >= kopt-1 || first)) {
	if (errmax < 1.0) {
	  exitflag=true;
	  break;
	}
	if (k == kmax || k == kopt+1) {
	  red=SAFE2/err[kkm];
	  break;
	}
	else if (k == kopt && alf[kopt-1][kopt] < err[kkm]) {
	  red=1.0/err[kkm];
	  break;
	}
	else if (kopt == kmax && alf[kkm][kmax-1] < err[kkm]) {
	  red=alf[kkm][kmax-1]*SAFE2/err[kkm];
	  break;
	}
	else if (alf[kkm][kopt] < err[kkm]) {
	  red=alf[kkm][kopt-1]/err[kkm];
	  break;
	}
      }
    }
    if (exitflag) break;
    red=std::min(red,REDMIN);
    red=std::max(red,REDMAX);
    h *= red;
    reduct=1;
  }
  xx=xnew;
  hdid=h;
  first=0;
  wrkmin=1.0e35;
  for (kk=0;kk<=kkm;kk++) {
    fact=std::max(err[kk],SCALMX);
    work=fact*a[kk+1];
    if (work < wrkmin) {
      scale=fact;
      wrkmin=work;
      kopt=kk+1;
    }
  }
  hnext=h/scale;
  if (kopt >= k && kopt != kmax && !reduct) {
    fact=std::max(scale/alf[kopt-1][kopt],SCALMX);
    if (a[kopt+1]*fact <= wrkmin) {
      hnext=h/fact;
      kopt++;
    }
  }
  
  /*	delete x_p;
	for (i=0; i<nv ; i++) delete d_p[i];
        delete d_p;*/
}

////////////////////////////////////////////////////////////////////////////////
//
void BSIntegrator1::mmid(const G4double y[], const G4double dydx[], 
			 const G4double xs, const G4double htot,
			 const G4int nstep, G4double yout[])
{
  G4int i,n;
  G4double x,swap,h2,h;
  
  G4int nvar=theEquation->GetNvar();
  G4double* ym = new G4double[nvar];
  G4double* yn =new G4double[nvar];
  h=htot/nstep; 
  //G4cout<<htot<<" htot "<<h<<" h "<<endl;
  if (IsTheEquationOfMotionForG4Tracking){
    dynamic_cast <BSEquationInG4*> (theEquation)
      ->SetDidCrossBoundaryDuringMmid(false);
  }			  
  // SetLastSubPosition in equation
  theEquation->SetLastSubX(xs);
  theEquation->SetLastSubPosition
    (G4ThreeVector(y[0],y[1],y[2]));
  
  for (i=0;i<nvar;i++) {
    ym[i]=y[i];
    yn[i]=y[i]+h*dydx[i];
  }
  x=xs+h;
  theEquation->derivs(x,yn,yout);
  // G4cout<<"ok "<<endl;
  
  
  // SetLastSubPosition in equation
  theEquation->SetLastSubX(x);
  theEquation->SetLastSubPosition
    (G4ThreeVector(yn[0],yn[1],yn[2]));
  
  
  h2=2.0*h;
  for (n=1;n<nstep;n++) {
    for (i=0;i<nvar;i++) {
      swap=ym[i]+h2*yout[i];
      ym[i]=yn[i];
      yn[i]=swap;
    }
    x += h;
    theEquation->derivs(x,yn,yout); 
    theEquation->SetLastSubX(x);
    theEquation->SetLastSubPosition
      (G4ThreeVector(yn[0],yn[1],yn[2]));
  }
  for (i=0;i<nvar;i++)
    yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
  
  delete ym;
  delete yn;
  // G4cout<<"ok "<<endl;
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void BSIntegrator1::pzextr(const G4int iest, const G4double xest, 
                          const G4double yest[], G4double yz[], G4double dy[])
{
  G4int j,k1;
  G4double q,f2,f1,delta;
  
  G4int nv=theEquation->GetNvar();
  std::vector< double > c;
  // = new G4double[nv];
  /*G4double* x=x_p;
    G4double** d=d_p;*/
  x_p[iest]=xest;
  for (j=0;j<nv;j++) dy[j]=yz[j]=yest[j];
  if (iest == 0) {
    for (j=0;j<nv;j++) d_p[j][0]=yest[j];
  } 
  else {
    for (j=0;j<nv;j++) c.push_back(yest[j]);
    for (k1=0;k1<iest;k1++) {
      delta=1.0/(x_p[iest-k1-1]-xest);
      f1=xest*delta;
      f2=x_p[iest-k1-1]*delta;
      for (j=0;j<nv;j++) {
	q=d_p[j][k1];
	d_p[j][k1]=dy[j];
	delta=c[j]-q;
	dy[j]=f1*delta;
	c[j]=f2*delta;
	yz[j] += dy[j];
      }
    }
    for (j=0;j<nv;j++) d_p[j][iest]=dy[j];
  }
  
}
