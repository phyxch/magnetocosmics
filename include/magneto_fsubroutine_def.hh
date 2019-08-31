// define equivalent of fortran common in c++
struct cctimereference
 {double year_ref;
  double day_ref;
  int nday;
 }; 
struct ccprevious
 {double tprev;
  };
struct ccgeopack
 {float st0,ct0,sl0,cl0,ctcl,stcl,ctsl,stsl,sfi,cfi,sps;
  float cps,shi,chi,hi,psi,xmut,a11,a21,a31,a12,a22,a32,a13,a23,a33,ds3;
  int k,iy;
  float cgst,sgst;
  float ba[6];
 };
struct ccgeopack_dbl
 {double st0,ct0,sl0,cl0,ctcl,stcl,ctsl,stsl,sfi,cfi,sps;
  double cps,shi,chi,hi,psi,xmut,a11,a21,a31,a12,a22,a32,a13,a23,a33,ds3;
  int k,iy;
  double cgst,sgst;
  double ba[6];
 };
struct ccgeocoef_dbl
 {double g10,g11,h11,g_10,g20,g21,g22,h21,h22;
 }; 
struct ccwhere_in_magnetopause2001
 {int location;
 }; 
 
  
struct ccbdip
 {double b0dip; 
  };   
struct ccigrf
 {float gg[105];
  float hh[105];
  float rec[105];
  int nm;};   
extern "C"
{extern void f_chang_coord__(int* , char* ,float* ,
        float* ,  float* ,  float* ,  float* , float* , float* );
 extern void igrf_to_cc__(float*,float*,float*,
                          float*,float*,float*);
 extern void  t89c__(int*,float*, float*,float*,float*,float*,
                                            float*,float*,float* );
 extern void  t89c_boberg__(int*,float*, float*,float*,float*,float*,
                                            float*,float*,float* );					    
 extern void  mcos_t96_01__(int*,float*, float*,float*,float*,float*,
                                            float*,float*,float* );			  	
 extern void  t01_01__(int*,float*, float*,float*,float*,float*,
                                            float*,float*,float* );
#ifdef  USE_TSY04
 extern void  t04_s__(int*,float*, float*,float*,float*,float*,
                                            float*,float*,float* );
#endif 
#ifdef  USE_PALEO
 extern void  calculate_cals7k_coef__(char*, double* );
 extern void  compute_cals7k_field_for_magcos__(double*,double*,double*,
                                            double*,double*,double*);
#endif 

 extern void mcos_magnetopause_2001__
            (double*, double*,double*,double*,double*,double*);
 
 extern struct cctimereference timereference_;
 extern struct ccprevious previous_; 
 extern struct ccgeopack geopack_;
 extern struct ccgeopack_dbl geopack_dbl__;
 extern struct ccbdip bdip_;
 extern struct ccgeocoef_dbl geocoef_dbl__;
 extern struct ccwhere_in_magnetopause2001 where_in_magnetopause2001__;	
 extern struct ccigrf igrfcc_;	  

 //extern void  mcos_t01_01__(int*,float*, float*,float*,float*,float*,
 //                                           float*,float*,float* );
 /*extern void igrf1_(float*,int*,float*,float*,float*,
                                        float*,float*,float*);
 extern void  read_storm__();
 extern void  julday2date_(double* ,int* ,int* ,int*,int*,int*,int*);
 extern void  recalc_(int*,int*,int*,int*,int*);
 extern void  recalc_dbl__(int*,int*,int*,int*,int*);
 extern void  t96_01__(int*,float*, float*,float*,float*,float*,

 extern void  t01_01_dbl__
             (int*,double*, double*,double*,double*,double*,
                                    double*,double*,double*);	
 extern void  t01_01_apot__
             (int*,double*, double*,double*,double*,double*,
                                    double*,double*,double*,
				     double*,double*,double*);
 extern void geogsm_(float*,float*,float*,float*,float*,float*,int*);
 extern void geogsm_dbl__(double*,double*,double*,double*,double*,double*,int*); 
 extern void smgsm_dbl__(double*,double*,double*,double*,double*,double*,int*);
 extern void igrf_(double*,int*,double*,double*,double*,
                                        double*,double*,double*);*/					
 
 }
