/* g_iondipole.c is meant to calculate the dipole moment of ionic systems
 * Written by Thomas Jansen 2018 */

/* Include various gromacs subroutines */

#include <math.h>
#include <ctype.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "princ.h"
#include "rmpbc.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "tpxio.h"
#include "physics.h"
#include "macros.h"
#include "gmx_fatal.h"
#include "index.h"
#include "gstat.h"
#include "string2.h"
#include "bondf.h"


#if OMP_PARALLEL
        #include "omp.h"
#endif

void calc_field(const char *fnTRX, atom_id **index, int gnx[], 
                t_topology *top, int nr_grps, int units,output_env_t oenv)
{
  matrix box;
  rvec d,lE;
  rvec *x;
  real *charges;
  t_pbc pbc;
  int       ePBC=-1;
  t_trxstatus *status;
  FILE *FH_field;
  FILE *FH_pos;
  FILE *FH_local;
  FILE *FH_charge;
  int molecules;
  int n,a,u,nr_frames;
  int natoms;
  real dist,idist,idist3,idist2;
  real t;
  real q;  
  real phi;
  rvec E;
  matrix gE;

  if (gnx[0] % units !=0){
    printf("The selected index group does not contain a whole number of molecules with %d atoms each.\n",units);
  }
  molecules=gnx[0]/units;

  FH_field=fopen("Field.bib","wb");
  FH_pos=fopen("Position.bib","wb");
  FH_local=fopen("Local.bib","wb"); // Contain information of box size etc.
  FH_charge=fopen("Charge.bib","wb");

  // Open Trajectory
  if ((natoms = read_first_x(oenv,&status,fnTRX,&t,&x,box)) == 0)
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");
    
  set_pbc(&pbc,ePBC,box); 
  // Save Charges
  charges=(float *)calloc(natoms,sizeof(float));
  for (n=0;n<natoms;n++){
    charges[n]=top->atoms.atom[n].q;
  }
  fwrite(&charges,sizeof(float),natoms,FH_charge);
  fclose(FH_charge);

  nr_frames=0;
  /*********** Start processing trajectory ***********/
  do 
  {
    // Write snapshot number to files
    fwrite(&nr_frames,sizeof(int),1,FH_field);
    fwrite(&nr_frames,sizeof(int),1,FH_pos);
    fwrite(&nr_frames,sizeof(int),1,FH_local);
    fwrite(&box,sizeof(float),9,FH_local);
    // Loop over molecules
    for (u=0;u<molecules;u++){
      // Initalize the field info
      phi=0;
      E[XX]=E[YY]=E[ZZ]=0;
      gE[XX][XX]=gE[YY][XX]=gE[XX][YY]=gE[ZZ][XX]=gE[XX][ZZ]=gE[YY][YY]=gE[YY][ZZ]=gE[ZZ][YY]=gE[ZZ][ZZ]=0;
      // Loop over atoms in molecule
      for (a=0;a<units;a++){
        // Write atom positions to file      
        fwrite(&x[index[0][u*units+a]],sizeof(float),3,FH_pos);
        // Loop over solvent atoms
        for (n=0;n<natoms;n++){
          if (n!=index[0][u*units+a]){
            q=top->atoms.atom[n].q;
            pbc_dx(&pbc,x[index[0][u*units+a]],x[n],d);
            dist=norm(d);
            idist=1.0/dist;
            idist2=idist*idist;
            idist3=idist2*idist;
            phi+=q*idist;
            svmul(q*idist3,d,lE);
            rvec_add(lE,E,E);
            // Add Electric field gradient here
            gE[XX][XX]+=(1.0-3.0*d[XX]*d[XX]*idist2)*q*idist3;
            gE[YY][YY]+=(1.0-3.0*d[YY]*d[YY]*idist2)*q*idist3;
            gE[ZZ][ZZ]+=(1.0-3.0*d[ZZ]*d[ZZ]*idist2)*q*idist3;
            gE[XX][YY]-=3.0*d[XX]*d[YY]*idist2*q*idist3;
            gE[XX][ZZ]-=3.0*d[XX]*d[ZZ]*idist2*q*idist3;
            gE[YY][ZZ]-=3.0*d[YY]*d[ZZ]*idist2*q*idist3;
          }
        }
        // Subtract contributions from within the same unit
        for (n=0;n<units;n++){
          if(n!=a){
            // Use charge with opposite sign
            q=-top->atoms.atom[index[0][u*units+n]].q;
            pbc_dx(&pbc,x[index[0][u*units+a]],x[index[0][u*units+n]],d);
            dist=norm(d);
            idist=1.0/dist;
            idist2=idist*idist;
            idist3=idist*idist2;
            phi+=q*idist;
            svmul(q*idist3,d,lE);
            rvec_add(lE,E,E);
            // Add Electric field gradient here
            gE[XX][XX]+=(1.0-3.0*d[XX]*d[XX]*idist2)*q*idist3;
            gE[YY][YY]+=(1.0-3.0*d[YY]*d[YY]*idist2)*q*idist3;
            gE[ZZ][ZZ]+=(1.0-3.0*d[ZZ]*d[ZZ]*idist2)*q*idist3;
            gE[XX][YY]-=3.0*d[XX]*d[YY]*idist2*q*idist3;
            gE[XX][ZZ]-=3.0*d[XX]*d[ZZ]*idist2*q*idist3;
            gE[YY][ZZ]-=3.0*d[YY]*d[ZZ]*idist2*q*idist3;
          }
        }
        gE[YY][XX]=gE[XX][YY];
        gE[ZZ][XX]=gE[XX][ZZ];
        gE[ZZ][YY]=gE[YY][ZZ];
        fwrite(&phi,sizeof(float),1,FH_field);
        fwrite(&E,sizeof(float),3,FH_field);
        fwrite(&gE,sizeof(float),6,FH_field);
      }
    }

    nr_frames++;
  } while (read_next_x(oenv,status,&t,natoms,x,box));
  fclose(FH_field);
  fclose(FH_pos);
  fclose(FH_local);
}
int main(int argc,char *argv[])
{

  const char *desc[] = {
    "Program to compute the electrical field on specified locations"
    "from the trajectory file."
  };

//  const char *bugs[] = {
//    "Discarding slices for integration should not be necessary."
//  };

  static int help=0;
  int cut=0;
      t_pargs pa [] = {
      { "-u",   FALSE, etINT, {&cut}, "Units"},
      { "-h",   FALSE, etINT, {&help}, 
        "Help." }};

  char      **grpname;                      /* groupnames                 */
  int       ngrps = 1,                      /* nr. of groups              */
            *ngx;                           /* sizes of groups            */
  int  ePBC;
  t_topology *top;                          /* topology                   */ 
  atom_id   **index;                        /* indices for all groups     */
  t_filenm  fnm[] = {                       /* files for g_order          */
    { efTRX, "-f", NULL,  ffREAD },         /* trajectory file            */
    { efNDX, NULL, NULL,  ffREAD },         /* index file                 */
    { efTPX, NULL, NULL,  ffREAD },         /* topology file              */
    { efXVG,"-of","field", ffWRITE },       /* xvgr output file           */
  };

  output_env_t oenv;

#define NFILE asize(fnm)
#define NPA asize(pa)
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL,&oenv);
  top = read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);

  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(ngx,ngrps);

  rd_index(ftp2fn(efNDX,NFILE,fnm),ngrps,ngx,index,grpname); 
  
  calc_field(ftp2fn(efTRX,NFILE,fnm),index,ngx,top,ngrps,cut,oenv);

  thanx(stderr);

  return 0;
}
