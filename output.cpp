#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "const.h"
#include "shared.hpp"
#include "ParticleSource.h"
#include "output.h"


void outputAbsorbtion(std::string filename, double circleCounter, int nParticles, int nBins, Bin absorbtionBins[])
{
  FILE    *df;
  double cumulated = 0.0;
  const char* formstr=     "%-14d%-17.10e%-17.10e%-17.10e%-17.10e%-17.10e%-17.10e\n";
  const char* titleformstr="# %-12s%-17s%-17s%-17s%-17s%-17s%-17s\n";

  df = fopen(filename.c_str(), "w"); if(df==NULL){ printf("Can't open output file. Maybe \"out\" folder does not exist?\n");exit(-1);}

  fprintf(df,"# depmsition distribution\n");		      
  fprintf(df,"# total number of landed particles = %f\n",circleCounter);		        
  fprintf(df, "# total fraction of landed particles = %f\n", circleCounter / double(nParticles * NUM_REDEP * SUBC));
  fprintf(df,titleformstr,"bin #","# part's","# Xe+","# Xe++","angle","surface","fraction");			    
  for (int i = 0; i < nBins; i++)
  {
    fprintf(df,formstr,
            i, absorbtionBins[i].nParticles, absorbtionBins[i].xe2p, absorbtionBins[i].xe2pp,
            i?.5*(absorbtionBins[i].start+absorbtionBins[i].end):absorbtionBins[i].start,
	    fabs(surf_sph_ring(absorbtionBins[i].start*M_PI/180, absorbtionBins[i].end*M_PI/180 )),
            absorbtionBins[i].nParticles / double(nParticles) / double(NUM_REDEP) / double(SUBC)
           );
    cumulated+=absorbtionBins[i].nParticles / double(nParticles)/double(NUM_REDEP)/double(SUBC);
  }

  fclose(df);
}

void outputEmission(std::string filename, int nParticles, int nBins, Bin emissionBins[])
{
  FILE    *df;

  const char* formstr=     "%-14d%-17.10e%-17.10e%-17.10e%-17.10e%-17.10e%-17.10e\n";
  const char* titleformstr="# %-12s%-17s%-17s%-17s%-17s%-17s%-17s\n";

  df = fopen(filename.c_str(), "w"); if(df==NULL){ printf("Can't open output file. Maybe \"out\" folder does not exist?\n");exit(-1);}

  fprintf(df,"# emissimn distribution\n");
  fprintf(df,titleformstr,"bin #","# part's","# Xe+","# Xe++","angle","surface","fraction");
  for (int i = 0; i < nBins; i++)
  {

    fprintf(df,formstr,
            i, emissionBins[i].nParticles, emissionBins[i].xe2p, emissionBins[i].xe2pp, 
	    .5*(emissionBins[i].start+emissionBins[i].end), 
	    fabs(surf_sph_ring( emissionBins[i].start*M_PI/180, emissionBins[i].end*M_PI/180)),
	    emissionBins[i].nParticles/nParticles 
	   );
  }
  fclose(df);
}

/**************************************** 2D 2D 2D **********************************************/
void outputAbsorbtion2DChan(std::string filename, double circleCounter, int nParticles, int nBinsphi, int nBinsz, Bin2D absorbtionBins[], double r)
{
  FILE    	*df;
  int 		k=0;
  double 	cumulated=0.0, dphi=0.0, dz=0.0;
  const char* 	formstr     ="%-14d%-17.10e%-17.10e%-17.10e%-17.10e%-17.10e%-17.10e%-17.10e\n";
  const char* 	titleformstr="# %-12s%-17s%-17s%-17s%-17s%-17s%-17s%-17s\n";

  df = fopen(filename.c_str(), "w"); if(df==NULL){ printf("Can't open output file. Maybe \"out\" folder does not exist?\n");exit(-1);}

  fprintf(df,"# deposition distribution inside thruster channel\n");		      
  fprintf(df,"# total number of landed particles = %f\n",circleCounter);		        
  fprintf(df, "# total fraction of landed particles = %f\n", circleCounter / double(nParticles * NUM_REDEP * SUBC));
  fprintf(df,titleformstr,"bin #","# part's","# Xe+","# Xe++","phi","z","surface","fraction");			    
  for (int i = 0; i < nBinsphi; i++)
  {
    for(int j=0; j<nBinsz; j++)
    {  
        k   = i*nBinsz + j;
	dphi= absorbtionBins[k].endI - absorbtionBins[k].startI;
	dz  = absorbtionBins[k].endJ - absorbtionBins[k].startJ;

        fprintf(df,formstr,
            k, absorbtionBins[k].nParticles, absorbtionBins[k].xe2p, absorbtionBins[k].xe2pp,
            absorbtionBins[k].startI+(dphi*0.5), absorbtionBins[k].startJ+(dz*0.5),
	    dz*2*M_PI*r*dphi/360,
                absorbtionBins[i].nParticles / double(nParticles) / double(NUM_REDEP) / double(SUBC)
           );
        cumulated+=absorbtionBins[i].nParticles / double(nParticles)/double(NUM_REDEP)/double(SUBC);
    }
    fprintf(df,"\n");
  }

  fclose(df);
}


void outputAbsorbtion2DBot(std::string filename, double circleCounter, int nParticles, int nBinsphi, int nBinsz, Bin2D absorbtionBins[])
{
  FILE    	*df;
  int 		k=0;
  double 	cumulated=0.0, dphi=0.0, dz=0.0;
  const char* 	formstr     ="%-14d%-17.10e%-17.10e%-17.10e%-17.10e%-17.10e%-17.10e%-17.10e\n";
  const char* 	titleformstr="# %-12s%-17s%-17s%-17s%-17s%-17s%-17s%-17s\n";

  df = fopen(filename.c_str(), "w"); if(df==NULL){ printf("Can't open output file. Maybe \"out\" folder does not exist?\n");exit(-1);}

  fprintf(df,"# deposition distribution inside thruster channel bottom\n");		      
  fprintf(df,"# total number of landed particles = %f\n",circleCounter);		        
  fprintf(df, "# total fraction of landed particles = %f\n", circleCounter / double(nParticles * NUM_REDEP * SUBC));
  fprintf(df,titleformstr,"bin #","# part's","# Xe+","# Xe++","phi","r","surface","fraction");			    
  for (int i = 0; i < nBinsphi; i++)
  {
    for(int j=0; j<nBinsz; j++)
    {  
        k   = i*nBinsz + j;
	dphi= absorbtionBins[k].endI - absorbtionBins[k].startI;
	dz  = absorbtionBins[k].endJ - absorbtionBins[k].startJ;
        fprintf(df,formstr,
            k, absorbtionBins[k].nParticles, absorbtionBins[k].xe2p, absorbtionBins[k].xe2pp,
            absorbtionBins[k].startI+(dphi*0.5), absorbtionBins[k].startJ+(dz*0.5),
	    ((absorbtionBins[k].endJ*absorbtionBins[k].endJ)-(absorbtionBins[k].startJ*absorbtionBins[k].startJ))*M_PI*dphi/360,
                absorbtionBins[i].nParticles / double(nParticles) / double(NUM_REDEP) / double(SUBC)
           );
        cumulated+=absorbtionBins[i].nParticles / double(nParticles)/double(NUM_REDEP)/double(SUBC);
    }
    fprintf(df,"\n");
  }

  fclose(df);
}