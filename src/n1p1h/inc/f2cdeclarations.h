#ifndef _f2cdeclarations_
#define _f2cdeclarations_

extern"C" {

  void constantsinitialization_(double *ma, int *irpa); 

  void initialization_(double *xmagev,int *irpa,double *rm,int *anti,double *type,int *contrlepton);
  
  void rpaparameters_(double *a0,double *a1,double *a2,double *a3,double *a4,double *a5,double *a6,double *a7,double *a8,double *a9,int *i);
 
  double  sigthnacho_(double *eingev,double *tmugev,double *cost);
  
  double  sigthbruno_(double *eingev,double *tmugev,double *cost,int *NRC,int *NPC,double *R,double *dP,double *vc);
  
  double sigthfree_(double *eingev,double *tmugev,double *cost);
  
  void comptmomentum_(double *vPnu);

  void selectlepton_(int *ilepton); 

  void setleptonmass_(int *ilepton,double *mass); 
    
  void vcinterpolate_(double *R,double *vc);

  double targetminmon_(double *q0,double *dq,double *dmt,double *dmp);

  void convolucion_(double *dncxp,double *dnca0p,double *dncxn,double *dnca0n, int *klave, int *lc, double *dxp, double *da0p,double *dxn,double *da0m);  
  
  void xro_(double *val); 
  
  void dsg20r_(double *val,double *rmax, int *mnr, double dr[], int *nr); 
  
  void drg20r_(double *val, double *r, int *vali, double df0[], double *f1); 
  
  double densq_(double *rp);
  
  extern struct{
    int nxro;
  } nxro_; 
  
  extern struct{
    int ncxro;
  } ncxro_; 
  
  extern struct {
    double drop; 
    double dron; 
  } vfpi_; 
  

  extern struct {
    double dxp; 
    double da0p; 
    double dro0p; 
    double dxn; 
    double da0n; 
    double dro0n; 
  } cteropn_; 
    
  extern struct {
    double dncxp; 
    double dnca0p; 
    double dronc0p; 
  } cteropnc_;

  extern struct {
    int klave; 
  } roferos_; 
  
  extern struct {
    double dzz; 
    double daa; 
 } nuc_; 
  
  extern struct {
    double dpi; 
    double hbarc; 
    double gf0; 
    double dmnu; 
    double dma; 
  } datos_; 
  
  extern struct {
    double dmneutrino; 
    double dmlepton; 
    double dmi; 
    double dmf; 
    double coscabibbo; 
    double dmuon; 
    double dmelectron;
    double dmtau;
    double xuma; 
  } datos2_; 
  
  extern struct {
    double vcd[2000];
    double dr[2000];
    int nr; 
  } vcfto_; 
  
  extern struct {
    double xma; 
  } mvalue_; 
  
  extern struct {
    int ieta; 
  } testieta_; 
  
  extern struct {
    double rmax; 
  } rinit_; 
  
  
  extern struct {
    double qvalue; 
  } qvalue_; 
  
  extern struct {
    int ipol; 
  } enregistre_; 
  
  extern struct {
    int ilin; 
    int mn; 
    int mnr; 
  } verif_; 
  
  extern struct {
    double dro; 
    double dro0; 
  } densidad_;


  extern struct{
    bool isinitialized; 
  } flags_; 

  extern struct{
    double  vPnu[4];
    double  vPmu[4];
    double  vq[4];
    double  vPn[4];
    double  vPp[4];
  } fourvectors_;
}

#define SelectLepton(a) selectlepton_(&a)

#endif 
