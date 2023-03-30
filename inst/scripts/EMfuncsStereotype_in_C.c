// Program to estimate the parameters of a Stereotype model using command "optim"
// Program created: 07-Jun-2012
// Last modification: 07-Jun-2012

#include <math.h>   // for math commands such as exp function
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <Rinternals.h>


void unpack_param_Stereo_nullmodel(
  double *parstart,
  int	 *q,
  double *mu,
  double *phi)
{  
  int k, index;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
      mu[k] = parstart[index];
      index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }  
}

void unpack_param_Stereo_roweffect(
  double *parstart,
  int	 *q,
  int	 *n,
  double *mu,
  double *phi,
  double *alpha)
{  
  int k, i, index;
  double sumalpha;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
      mu[k] = parstart[index];
      index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }  
  
  sumalpha=0;
  for (i=0; i<(*n-1); i++)
  {
    alpha[i]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*n-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (i=1; i<*n; n++)
//  {
//      alpha[i] = parstart[index];
//      index++;
//  } 

}

void unpack_param_Stereo_columneffect(
  double *parstart,
  int	 *q,
  int	 *m,
  double *mu,
  double *phi,
  double *beta)
{  
  int k, j, index;
  double sumbeta;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
      mu[k] = parstart[index];
      index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index]; 
    index++;
  }  
  
  sumbeta=0;
  for (j=0; j<(*m-1); j++)
  {
    beta[j]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*m-1)]=-sumbeta;

//  This is the corner-point parametrization constraint is the beta[0]=0
//  beta[0]=0;
//  for (j=1; j<*m; m++)
//  {
//      beta[j] = parstart[index];
//      index++;
//  } 

}


void unpack_param_Stereo_maineffects(
  double *parstart,
  int	 *q,
  int	 *n, 
  int	 *m,
  double *mu,
  double *phi,
  double *alpha,
  double *beta)
{  
  int k, i, j, index;
  double sumalpha, sumbeta;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
      mu[k] = parstart[index];
      index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }  
  
  
  sumalpha=0;
  for (i=0; i<(*n-1); i++)
  {
    alpha[i]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*n-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (i=1; i<*n; n++)
//  {
//      alpha[i] = parstart[index];
//      index++;
//  } 
  
  sumbeta=0;
  for (j=0; j<(*m-1); j++)
  {
    beta[j]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*m-1)]=-sumbeta;

//  This is the corner-point parametrization constraint is the beta[0]=0
//  beta[0]=0;
//  for (j=1; j<*m; m++)
//  {
//      beta[j] = parstart[index];
//      index++;
//  } 

}


void unpack_param_Stereo_RowCluster(
  double *parstart,
  int	 *m,
  int 	 *R,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *arraygamma)
{
  int j, r, k, index;
  double matrixgamma[*R][*m];
  double sumcol, sumrow, sumalpha, sumbeta; 
  int indexarraygamma;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
      mu[k] = parstart[index];
      index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (r=0; r<(*R-1); r++)
  {
    alpha[r]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*R-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (r=1; r<*R; r++)
//  {
//      alpha[r] = parstart[index];
//      index++;
//  } 

  sumbeta=0;
  for (j=0; j<(*m-1); j++)
  {
    beta[j]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*m-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (j=1; j<*m; j++)
//  {
//      beta[j] = parstart[index];
//      index++;
//  } 
  
  //Gamma is build up in several steps: sum-to-zero constraints
  
  //Step1: Initialize the matrix to 0's
  for (r=0; r<*R; r++)
  {
    for (j=0; j<*m; j++)
    {
      matrixgamma[r][j] = 0;
    } 
  }
  // Step2: fill the par.start values
  for (r=1; r<*R; r++)
  {
    for (j=0; j<(*m-1); j++)
    {
      matrixgamma[r][j] = parstart[index];
      index++;
    } 
  }
  // Step3: fill with the sum of columns and rows
  
  for (r=1; r<*R; r++)
  {
    sumrow=0.0;
    for (j=0; j<(*m-1); j++)
    {
      sumrow += matrixgamma[r][j];
    }
    matrixgamma[r][(*m-1)]=-sumrow;
  }
  
  for (j=0; j<*m; j++)
  {
    sumcol=0;
    for (r=1; r<*R; r++)
    {
      sumcol += matrixgamma[r][j];
    }
    matrixgamma[0][j]=-sumcol;
  }
  // step4: We create the array of gamma (by columns)
  indexarraygamma=0;
  for (j=0; j<*m; j++)
  {
    for (r=0; r<*R; r++)
    {
	arraygamma[indexarraygamma] = matrixgamma[r][j];
	indexarraygamma++;
    }
  }
}

void unpack_param_Stereo_RowCluster_withpiR(
  double *parstart,
  int	 *m,
  int 	 *R,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *arraygamma,
  double *piR)
{
  int j, r, k, index;
  double matrixgamma[*R][*m];
  double sumcol, sumrow, sumalpha, sumbeta, sumpiR; 
  int indexarraygamma;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
      mu[k] = parstart[index];
      index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (r=0; r<(*R-1); r++)
  {
    alpha[r]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*R-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (r=1; r<*R; r++)
//  {
//      alpha[r] = parstart[index];
//      index++;
//  } 

  sumbeta=0;
  for (j=0; j<(*m-1); j++)
  {
    beta[j]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*m-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (j=1; j<*m; j++)
//  {
//      beta[j] = parstart[index];
//      index++;
//  } 
  
  //Gamma is build up in several steps: sum-to-zero constraints
  
  //Step1: Initialize the matrix to 0's
  for (r=0; r<*R; r++)
  {
    for (j=0; j<*m; j++)
    {
      matrixgamma[r][j] = 0;
    } 
  }
  // Step2: fill the par.start values
  for (r=1; r<*R; r++)
  {
    for (j=0; j<(*m-1); j++)
    {
      matrixgamma[r][j] = parstart[index];
      index++;
    } 
  }
  // Step3: fill with the sum of columns and rows
  
  for (r=1; r<*R; r++)
  {
    sumrow=0.0;
    for (j=0; j<(*m-1); j++)
    {
      sumrow += matrixgamma[r][j];
    }
    matrixgamma[r][(*m-1)]=-sumrow;
  }
  
  for (j=0; j<*m; j++)
  {
    sumcol=0;
    for (r=1; r<*R; r++)
    {
      sumcol += matrixgamma[r][j];
    }
    matrixgamma[0][j]=-sumcol;
  }
  // step4: We create the array of gamma (by columns)
  indexarraygamma=0;
  for (j=0; j<*m; j++)
  {
    for (r=0; r<*R; r++)
    {
	arraygamma[indexarraygamma] = matrixgamma[r][j];
	indexarraygamma++;
    }
  }
  
  sumpiR=0;
  for (r=0; r<(*R-1); r++)
  {
    piR[r]=parstart[index];
    sumpiR += parstart[index];
    index++;
  }
  piR[(*R-1)]=1-sumpiR;
}

void unpack_param_Stereo_ColumnCluster(
  double *parstart,
  int	 *n,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *arraygamma)
{
  int i, c, k, index;
  double matrixgamma[*n][*C];
  double sumcol, sumrow, sumalpha, sumbeta; 
  int indexarraygamma;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (i=0; i<(*n-1); i++)
  {
    alpha[i]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*n-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (i=1; i<*n; i++)
//  {
//      alpha[i] = parstart[index];
//      index++;
//  } 

  sumbeta=0;
  for (c=0; c<(*C-1); c++)
  {
    beta[c]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*C-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (c=1; c<*C; c++)
//  {
//      beta[c] = parstart[index];
//      index++;
//  } 

  //Gamma is build up in several steps: sum-to-zero constraints
  
  //Step1: Initialize the matrix to 0's
  for (i=0; i<*n; i++)
  {
    for (c=0; c<*C; c++)
    {
      matrixgamma[i][c] = 0;
    } 
  }
  // Step2: fill the par.start values
  for (i=1; i<*n; i++)
  {
    for (c=0; c<(*C-1); c++)
    {
      matrixgamma[i][c] = parstart[index];
      index++;
    } 
  }
  // Step3: fill with the sum of columns and rows
  
  for (i=1; i<*n; i++)
  {
    sumrow=0.0;
    for (c=0; c<(*C-1); c++)
    {
      sumrow += matrixgamma[i][c];
    }
    matrixgamma[i][(*C-1)]=-sumrow;
  }
  
  for (c=0; c<*C; c++)
  {
    sumcol=0;
    for (i=1; i<*n; i++)
    {
      sumcol += matrixgamma[i][c];
    }
    matrixgamma[0][c]=-sumcol;
  }
  // step4: We create the array of gamma (by columns)
  indexarraygamma=0;
  for (c=0; c<*C; c++)
  {
    for (i=0; i<*n; i++)
    {
	arraygamma[indexarraygamma] = matrixgamma[i][c];
	indexarraygamma++;
    }
  } 
  
}

void unpack_param_Stereo_ColumnCluster_withkappaC(
  double *parstart,
  int	 *n,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *arraygamma,
  double *kappaC)
{
  int i, c, k, index;
  double matrixgamma[*n][*C];
  double sumcol, sumrow, sumalpha, sumbeta, sumkappaC; 
  int indexarraygamma;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (i=0; i<(*n-1); i++)
  {
    alpha[i]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*n-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (i=1; i<*n; i++)
//  {
//      alpha[i] = parstart[index];
//      index++;
//  } 

  sumbeta=0;
  for (c=0; c<(*C-1); c++)
  {
    beta[c]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*C-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (c=1; c<*C; c++)
//  {
//      beta[c] = parstart[index];
//      index++;
//  } 

  //Gamma is build up in several steps: sum-to-zero constraints
  
  //Step1: Initialize the matrix to 0's
  for (i=0; i<*n; i++)
  {
    for (c=0; c<*C; c++)
    {
      matrixgamma[i][c] = 0;
    } 
  }
  // Step2: fill the par.start values
  for (i=1; i<*n; i++)
  {
    for (c=0; c<(*C-1); c++)
    {
      matrixgamma[i][c] = parstart[index];
      index++;
    } 
  }
  // Step3: fill with the sum of columns and rows
  
  for (i=1; i<*n; i++)
  {
    sumrow=0.0;
    for (c=0; c<(*C-1); c++)
    {
      sumrow += matrixgamma[i][c];
    }
    matrixgamma[i][(*C-1)]=-sumrow;
  }
  
  for (c=0; c<*C; c++)
  {
    sumcol=0;
    for (i=1; i<*n; i++)
    {
      sumcol += matrixgamma[i][c];
    }
    matrixgamma[0][c]=-sumcol;
  }
  // step4: We create the array of gamma (by columns)
  indexarraygamma=0;
  for (c=0; c<*C; c++)
  {
    for (i=0; i<*n; i++)
    {
	arraygamma[indexarraygamma] = matrixgamma[i][c];
	indexarraygamma++;
    }
  } 
  sumkappaC=0;
  for (c=0; c<(*C-1); c++)
  {
    kappaC[c]=parstart[index];
    sumkappaC += parstart[index];
    index++;
  }
  kappaC[(*C-1)]=1-sumkappaC;
  
}

void unpack_param_Stereo_BiCluster(
  double *parstart,
  int	 *R,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *arraygamma)
{
  int r, c, k, index;
  double matrixgamma[*R][*C];
  double sumcol, sumrow, sumalpha, sumbeta; 
  int indexarraygamma;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (r=0; r<(*R-1); r++)
  {
    alpha[r]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*R-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (r=1; r<*R; r++)
//  {
//      alpha[r] = parstart[index];
//      index++;
//  } 
  sumbeta=0;
  for (c=0; c<(*C-1); c++)
  {
    beta[c]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*C-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (c=1; c<*C; c++)
//  {
//      beta[c] = parstart[index];
//      index++;
//  } 
    //Gamma is build up in several steps: sum-to-zero constraints
  
  //Step1: Initialize the matrix to 0's
  for (r=0; r<*R; r++)
  {
    for (c=0; c<*C; c++)
    {
      matrixgamma[r][c] = 0;
    } 
  }
  // Step2: fill the par.start values
  for (r=1; r<*R; r++)
  {
    for (c=0; c<(*C-1); c++)
    {
      matrixgamma[r][c] = parstart[index];
      index++;
    } 
  }
  // Step3: fill with the sum of columns and rows
  
  for (r=1; r<*R; r++)
  {
    sumrow=0.0;
    for (c=0; c<(*C-1); c++)
    {
      sumrow += matrixgamma[r][c];
    }
    matrixgamma[r][(*C-1)]=-sumrow;
  }
  
  for (c=0; c<*C; c++)
  {
    sumcol=0;
    for (r=1; r<*R; r++)
    {
      sumcol += matrixgamma[r][c];
    }
    matrixgamma[0][c]=-sumcol;
  }
  // step4: We create the array of gamma (by columns)
  indexarraygamma=0;
  for (c=0; c<*C; c++)
  {
    for (r=0; r<*R; r++)
    {
	arraygamma[indexarraygamma] = matrixgamma[r][c];
	indexarraygamma++;
    }
  }
  
}


void unpack_param_Stereo_RowCluster_without_iterations(
  double *parstart,
  int	 *m,
  int 	 *R,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta)
{
  int j, r, k, index;
  double sumalpha, sumbeta; 
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (r=0; r<(*R-1); r++)
  {
    alpha[r]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*R-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (r=1; r<*R; r++)
//  {
//      alpha[r] = parstart[index];
//      index++;
//  } 

  sumbeta=0;
  for (j=0; j<(*m-1); j++)
  {
    beta[j]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*m-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (j=1; j<*m; j++)
//  {
//      beta[j] = parstart[index];
//      index++;
//  } 
  
  
}

void unpack_param_Stereo_RowCluster_without_iterations_withpiR(
  double *parstart,
  int	 *m,
  int 	 *R,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *piR)
{
  int j, r, k, index;
  double sumalpha, sumbeta, sumpiR; 
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (r=0; r<(*R-1); r++)
  {
    alpha[r]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*R-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (r=1; r<*R; r++)
//  {
//      alpha[r] = parstart[index];
//      index++;
//  } 

  sumbeta=0;
  for (j=0; j<(*m-1); j++)
  {
    beta[j]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*m-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (j=1; j<*m; j++)
//  {
//      beta[j] = parstart[index];
//      index++;
//  } 
  sumpiR=0;
  for (r=0; r<(*R-1); r++)
  {
    piR[r]=parstart[index];
    sumpiR += parstart[index];
    index++;
  }
  piR[(*R-1)]=1-sumpiR;
  
}

void unpack_param_Stereo_ColumnCluster_without_iterations_withkappaC(
  double *parstart,
  int	 *n,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *kappaC)
{
  int i, c, k, index;
  double sumalpha, sumbeta, sumkappaC; 
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (i=0; i<(*n-1); i++)
  {
    alpha[i]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*n-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (i=1; i<*n; i++)
//  {
//      alpha[i] = parstart[index];
//      index++;
//  } 

  sumbeta=0;
  for (c=0; c<(*C-1); c++)
  {
    beta[c]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*C-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (c=1; c<*C; c++)
//  {
//      beta[c] = parstart[index];
//      index++;
//  } 
  sumkappaC=0;
  for (c=0; c<(*C-1); c++)
  {
    kappaC[c]=parstart[index];
    sumkappaC += parstart[index];
    index++;
  }
  kappaC[(*C-1)]=1-sumkappaC;
  
}


void unpack_param_Stereo_ColumnCluster_without_iterations(
  double *parstart,
  int	 *n,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta)
{
  int i, c, k, index;
  double sumalpha, sumbeta; 
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (i=0; i<(*n-1); i++)
  {
    alpha[i]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*n-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (i=1; i<*n; i++)
//  {
//      alpha[i] = parstart[index];
//      index++;
//  } 

  sumbeta=0;
  for (c=0; c<(*C-1); c++)
  {
    beta[c]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*C-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (c=1; c<*C; c++)
//  {
//      beta[c] = parstart[index];
//      index++;
//  } 
  
  
}

void unpack_param_Stereo_BiCluster_without_iterations(
  double *parstart,
  int	 *R,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta)
{
  int r, c, k, index;
  double sumalpha, sumbeta; 
  
  index=0;
 
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
 
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (r=0; r<(*R-1); r++)
  {
    alpha[r]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*R-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (r=1; r<*R; r++)
//  {
//      alpha[r] = parstart[index];
//      index++;
//  } 
  sumbeta=0;
  for (c=0; c<(*C-1); c++)
  {
    beta[c]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*C-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (c=1; c<*C; c++)
//  {
//      beta[c] = parstart[index];
//      index++;
//  } 

  
}


void unpack_param_Stereo_BiCluster_without_iterations_withpiR_and_kappaC(
  double *parstart,
  int	 *R,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *piR,
  double *kappaC)
{
  int r, c, k, index;
  double sumalpha, sumbeta, sumpiR, sumkappaC; 
  
  index=0;
 
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
 
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (r=0; r<(*R-1); r++)
  {
    alpha[r]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*R-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (r=1; r<*R; r++)
//  {
//      alpha[r] = parstart[index];
//      index++;
//  } 
  sumbeta=0;
  for (c=0; c<(*C-1); c++)
  {
    beta[c]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*C-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (c=1; c<*C; c++)
//  {
//      beta[c] = parstart[index];
//      index++;
//  } 
  sumpiR=0;
  for (r=0; r<(*R-1); r++)
  {
    piR[r]=parstart[index];
    sumpiR += parstart[index];
    index++;
  }
  piR[(*R-1)]=1-sumpiR;
  
  sumkappaC=0;
  for (c=0; c<(*C-1); c++)
  {
    kappaC[c]=parstart[index];
    sumkappaC += parstart[index];
    index++;
  }
  kappaC[(*C-1)]=1-sumkappaC;
  
}

void unpack_param_Stereo_BiCluster_withpiR_and_kappaC(
  double *parstart,
  int	 *R,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *arraygamma,
  double *piR,
  double *kappaC)
{
  int r, c, k, index;
  double matrixgamma[*R][*C];
  double sumcol, sumrow, sumalpha, sumbeta, sumpiR, sumkappaC; 
  int indexarraygamma;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (r=0; r<(*R-1); r++)
  {
    alpha[r]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*R-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (r=1; r<*R; r++)
//  {
//      alpha[r] = parstart[index];
//      index++;
//  } 
  sumbeta=0;
  for (c=0; c<(*C-1); c++)
  {
    beta[c]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*C-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (c=1; c<*C; c++)
//  {
//      beta[c] = parstart[index];
//      index++;
//  } 
    //Gamma is build up in several steps: sum-to-zero constraints
  
  //Step1: Initialize the matrix to 0's
  for (r=0; r<*R; r++)
  {
    for (c=0; c<*C; c++)
    {
      matrixgamma[r][c] = 0;
    } 
  }
  // Step2: fill the par.start values
  for (r=1; r<*R; r++)
  {
    for (c=0; c<(*C-1); c++)
    {
      matrixgamma[r][c] = parstart[index];
      index++;
    } 
  }
  // Step3: fill with the sum of columns and rows
  
  for (r=1; r<*R; r++)
  {
    sumrow=0.0;
    for (c=0; c<(*C-1); c++)
    {
      sumrow += matrixgamma[r][c];
    }
    matrixgamma[r][(*C-1)]=-sumrow;
  }
  
  for (c=0; c<*C; c++)
  {
    sumcol=0;
    for (r=1; r<*R; r++)
    {
      sumcol += matrixgamma[r][c];
    }
    matrixgamma[0][c]=-sumcol;
  }
  // step4: We create the array of gamma (by columns)
  indexarraygamma=0;
  for (c=0; c<*C; c++)
  {
    for (r=0; r<*R; r++)
    {
	arraygamma[indexarraygamma] = matrixgamma[r][c];
	indexarraygamma++;
    }
  }
  sumpiR=0;
  for (r=0; r<(*R-1); r++)
  {
    piR[r]=parstart[index];
    sumpiR += parstart[index];
    index++;
  }
  piR[(*R-1)]=1-sumpiR;
  
  sumkappaC=0;
  for (c=0; c<(*C-1); c++)
  {
    kappaC[c]=parstart[index];
    sumkappaC += parstart[index];
    index++;
  }
  kappaC[(*C-1)]=1-sumkappaC;
  
}


void unpack_param_Stereo_RowCluster_rRcC1(
  double *parstart,
  int 	 *R,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha)
{
  int r, k, index;
  double sumalpha; 
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (r=0; r<(*R-1); r++)
  {
    alpha[r]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*R-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (r=1; r<*R; r++)
//  {
//      alpha[r] = parstart[index];
//      index++;
//  } 

  
}

void unpack_param_Stereo_RowCluster_rRcC1_withpiR(
  double *parstart,
  int 	 *R,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *piR
)
{
  int r, k, index;
  double sumalpha, sumpiR; 
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumalpha=0;
  for (r=0; r<(*R-1); r++)
  {
    alpha[r]=parstart[index];
    sumalpha += parstart[index];
    index++;
  }
  alpha[(*R-1)]=-sumalpha;

//  This is the corner-point parametrization constraint is the alpha[0]=0
//  alpha[0]=0;
//  for (r=1; r<*R; r++)
//  {
//      alpha[r] = parstart[index];
//      index++;
//  } 
  sumpiR=0;
  for (r=0; r<(*R-1); r++)
  {
    piR[r]=parstart[index];
    sumpiR += parstart[index];
    index++;
  }
  piR[(*R-1)]=1-sumpiR;
  
}

void unpack_param_Stereo_ColumnCluster_rR1cC_withkappaC(
  double *parstart,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *beta,
  double *kappaC)
{
  int c, k, index;
  double sumbeta, sumkappaC;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumbeta=0;
  for (c=0; c<(*C-1); c++)
  {
    beta[c]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*C-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (c=1; c<*C; c++)
//  {
//      beta[c] = parstart[index];
//      index++;
//  } 
  sumkappaC=0;
  for (c=0; c<(*C-1); c++)
  {
    kappaC[c]=parstart[index];
    sumkappaC += parstart[index];
    index++;
  }
  kappaC[(*C-1)]=1-sumkappaC;
  
}



void unpack_param_Stereo_ColumnCluster_rR1cC(
  double *parstart,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *beta)
{
  int c, k, index;
  double sumbeta;
  
  index=0;
  
  mu[0]=0;
  for (k=1; k<*q; k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
  
  phi[0]=0;
  phi[(*q-1)]=1;
  for (k=1; k<(*q-1); k++)
  {
    phi[k] = parstart[index];
    index++;
  }
  
  sumbeta=0;
  for (c=0; c<(*C-1); c++)
  {
    beta[c]=parstart[index];
    sumbeta += parstart[index];
    index++;
  }
  beta[(*C-1)]=-sumbeta;
//  This is the corner-point parametrization constraint is the beta[0]=0  
//  beta[0]=0;
//  for (c=1; c<*C; c++)
//  {
//      beta[c] = parstart[index];
//      index++;
//  } 
  
  
}




void become_matrix(
  int *n,
  int *m,
  double *arraydata,
  double data[*n][*m])
{
  int i, k;
  
  for (i=0; i<*n; i++)
  {
    for (k=0; k<*m; k++)
    {
      data[i][k] = arraydata[(*n * k) + i];
    }
  }
}

// void become_matrix_3D(
//   int *n,
//   int *m,
//   int *g,
//   double *arraydata,
//   double data[*n][*m][*g])
// {
//   int i, k, l;
//   
//   for (i=0; i<*n; i++)
//   {
//     for (k=0; k<*m; k++)
//     {
//       for (l=0; l<*g; l++)
//       {
// 	data[i][k][l] = arraydata[(*n * *m * l)+(*n * k) + i];
// 
//       }
//     }
//   }
// }

void recover_mix_from_S(
    int *G,
    double *S) // This is the result parameter. 
{
    int g, l;
    double sumnewmix;
    double newmix[*G];
    
    
    newmix[0]=1/(1+exp(-S[0]));
   
    sumnewmix=newmix[0];
    for (g=1; g<(*G-1); g++)
    {
	 newmix[g]=1/(1+exp(-S[g]));
	 for (l=0; l<g; l++)
	 {
	   newmix[g] *= (1-(1/(1+exp(-S[l]))));
	 }
	 sumnewmix += newmix[g];
    }
    newmix[(*G-1)]=1-sumnewmix;
    
    for (g=0; g<*G; g++)
    {
      S[g]=newmix[g];     
    } 
}


void repar_problem( 
    int *q,
    double *phi) // This is the result parameter
{
    int i, k;
    double nu2, sum;
    double phiAux[*q];
    
    
    nu2=phi[1];
   
    phiAux[0]=phi[0];
    phiAux[1]=1/(1+exp(-nu2));
    phiAux[(*q-1)]=phi[(*q-1)];
    
    if (*q >= 4)
    {    
      for (k=2; k<(*q-1); k++)
      {
	 sum=nu2;
	 for (i=2; i<(k+1); i++)
	 {
	    sum = sum + exp(phi[i]);
	 }
	 phiAux[k]=1/(1+exp(-sum));
      }
    }
    
    for (i=0; i<*q; i++)
    {
      phi[i]=phiAux[i];     
    }  
}


//this is no a "external" function
double compute_prob_Stereo_BiCluster( 
  int r,
  int c,
  int k,
  int *n,
  int *q,
  int *R,
  int *C,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *arraygamma)
{
  int l;
  double sumatheta;
  double theta[*q];
  double gamma[*R][*C];
  
  // we create the matrix of our parameter gamma
  become_matrix(R,C,arraygamma,gamma);
  
  sumatheta=0;
  for (l=0; l<*q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*(alpha[r]+beta[c]+gamma[r][c])));
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}




//this is no a "external" function
double compute_prob_Stereo_BiCluster_without_iterations( 
  int r,
  int c,
  int k,
  int *n,
  int *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta)
{
  int l;
  double sumatheta;
  double theta[*q];
  

  sumatheta=0;
  for (l=0; l<*q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*(alpha[r]+beta[c])));
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}




//this is no a "external" function
double compute_prob_Stereo_RowCluster( 
  int r,
  int j,
  int k,
  int *R,
  int *m,
  int *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *arraygamma)
{
  int l;
  double sumatheta;
  double theta[*q];
  double gamma[*R][*m];
  
  // we create the matrix of our parameter gamma
  become_matrix(R,m,arraygamma,gamma);
  
  sumatheta=0;
  for (l=0; l<*q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*(alpha[r]+beta[j]+gamma[r][j])));
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}

//this is no a "external" function
double compute_prob_Stereo_ColumnCluster( 
  int c,
  int i,
  int k,
  int *C,
  int *n,
  int *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta,
  double *arraygamma)
{
  int l;
  double sumatheta;
  double theta[*q];
  double gamma[*n][*C];
  
  // we create the matrix of our parameter gamma
  become_matrix(n,C,arraygamma,gamma);

  sumatheta=0;
  for (l=0; l<*q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*(alpha[i]+beta[c]+gamma[i][c])));
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}



//this is no a "external" function
double compute_prob_Stereo_ColumnCluster_without_iterations( 
  int c,
  int i,
  int k,
  int *C,
  int *n,
  int *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta)
{
  int l;
  double sumatheta;
  double theta[*q];
  

  sumatheta=0;
  for (l=0; l<*q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*(alpha[i]+beta[c])));
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}

//this is no a "external" function
double compute_prob_Stereo_RowCluster_without_iterations( 
  int r,
  int j,
  int k,
  int *R,
  int *m,
  int *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta)
{
  int l;
  double sumatheta;
  double theta[*q];
  

  sumatheta=0;
  for (l=0; l<*q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*(alpha[r]+beta[j])));
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}


//this is no a "external" function
double compute_prob_Stereo_RowCluster_rRcC1( 
  int r,
  int k,
  int *R,
  int *q,
  double *mu,
  double *phi,
  double *alpha)
{
  int l;
  double sumatheta;
  double theta[*q];
  

  sumatheta=0;
  for (l=0; l<*q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*(alpha[r])));
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}

//this is no a "external" function
double compute_prob_Stereo_ColumnCluster_rR1cC( 
  int c,
  int k,
  int *C,
  int *q,
  double *mu,
  double *phi,
  double *beta)
{
  int l;
  double sumatheta;
  double theta[*q];
  

  sumatheta=0;
  for (l=0; l<*q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*beta[c]));
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}


//this is no a "external" function
double compute_prob_Stereo_nullmodel(
  int k,
  int *q,
  double *mu,
  double *phi)
{
  int l;
  double sumatheta;
  double theta[*q];
  
  sumatheta=0;
  for (l=0; l < *q ; l++)
  {
    theta[l]=exp(mu[l]+phi[l]);
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}

//this is no a "external" function
double compute_prob_Stereo_roweffectmodel(
  int i,
  int k,
  int *n,
  int *q,
  double *mu,
  double *phi,
  double *alpha)
{
  int l;
  double sumatheta;
  double theta[*q];
  
  sumatheta=0;
  for (l=0; l < *q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*alpha[i]));
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}

//this is no a "external" function
double compute_prob_Stereo_columneffectmodel(
  int j,
  int k,
  int *m,
  int *q,
  double *mu,
  double *phi,
  double *beta)
{
  int l;
  double sumatheta;
  double theta[*q];
  
  sumatheta=0;
  for (l=0; l < *q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*beta[j]));
    sumatheta += theta[l];
  }  
  return(theta[k]/sumatheta);
}

//this is no a "external" function
double compute_prob_Stereo_maineffectsmodel(
  int i,
  int j,
  int k,
  int *n,
  int *m,
  int *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta)
{
  int l;
  double sumatheta;
  double theta[*q];
  
  sumatheta=0;
  for (l=0; l < *q ; l++)
  {
    theta[l]=exp(mu[l]+(phi[l]*(alpha[i]+beta[j])));
    sumatheta += theta[l];
  }
  return(theta[k]/sumatheta);
}



void Q_Stereo_ColumnCluster_inC(
    double *parstart,
    double *arraydata,
    double *arrayx,
    int *n, 
    int *m, 
    int *q, 
    int *C,
    int *reparC,
    double *Qstereo)
{
  double data[*n][*m];
  double x[*m][*C];
  double mu[*q], phi[*q], alpha[*n], beta[*C], arraygamma[(*n * *C)];
  double theta; 
  int i, j, k, c;
  
  unpack_param_Stereo_ColumnCluster(parstart,n,C,q,mu,phi,alpha,beta,arraygamma);
 
  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  // we create the matrix of the expected value z
  become_matrix(m,C,arrayx,x);
 
  // Ask if we have to reparameterizer the phi's (ordered)
  if (*reparC == 1)
  {    
      repar_problem(q,phi);
  }
  //Compute the function to maximize
  *Qstereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	for (c=0; c < *C; c++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_ColumnCluster(c, i, k, C, n, q, mu, phi, alpha, beta, arraygamma);
	    *Qstereo += (log(theta)*x[j][c]);
	  } 
	} 
      }
    }   
  } 
  *Qstereo = (-1) * (*Qstereo); 
}





void Q_Stereo_RowCluster_inC(
    double *parstart,
    double *arraydata,
    double *arrayz,
    int *n, 
    int *m, 
    int *q, 
    int *R,
    int *repar,
    double *Qstereo)
{
  double data[*n][*m];
  double z[*n][*R];
  double mu[*q], phi[*q], alpha[*R], beta[*m], arraygamma[(*R * *m)];
  double theta;
  int i, j, k, r;
  
  unpack_param_Stereo_RowCluster(parstart,m,R,q,mu,phi,alpha,beta,arraygamma);
 
  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  // we create the matrix of the expected value z
  become_matrix(n,R,arrayz,z);
 
  // Ask if we have to reparameterizer the phi's (ordered)
  if (*repar == 1)
  {    
      repar_problem(q,phi);
  }

  //Compute the function to maximize
  *Qstereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
     for (k=0; k<*q; k++)
     {
	for (r=0; r<*R; r++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_RowCluster(r, j, k, R, m, q, mu, phi, alpha, beta, arraygamma);
	    *Qstereo += (log(theta)*z[i][r]);
	  } 
	}
       
      }
    }
    
  } 
  *Qstereo = (-1) * (*Qstereo); 
}


void Q_Stereo_ColumnCluster_without_iterations_inC(
    double *parstart,
    double *arraydata,
    double *arrayx,
    int *n, 
    int *m, 
    int *q, 
    int *C,
    int *reparC,
    double *Qstereo)
{
  double data[*n][*m];
  double x[*m][*C];
  double mu[*q], phi[*q], alpha[*n], beta[*C];
  double theta; 
  int i, j, k, c;
  
  unpack_param_Stereo_ColumnCluster_without_iterations(parstart,n,C,q,mu,phi,alpha,beta);
 
  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  // we create the matrix of the expected value z
  become_matrix(m,C,arrayx,x);
 
  // Ask if we have to reparameterizer the phi's (ordered)
  if (*reparC == 1)
  {    

      repar_problem(q,phi);
  }
  //Compute the function to maximize
  *Qstereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	for (c=0; c < *C; c++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_ColumnCluster_without_iterations(c, i, k, C, n, q, mu, phi, alpha, beta);
	    *Qstereo += (log(theta)*x[j][c]);
	  } 
	} 
      }
    }   
  } 
  *Qstereo = (-1) * (*Qstereo); 
}


void Q_Stereo_RowCluster_without_iterations_inC(
    double *parstart,
    double *arraydata,
    double *arrayz,
    int *n, 
    int *m, 
    int *q, 
    int *R,
    int *reparC,
    double *Qstereo)
{
  double data[*n][*m];
  double z[*n][*R];
  double mu[*q], phi[*q], alpha[*R], beta[*m];
  double theta; 
  int i, j, k, r;
  
  unpack_param_Stereo_RowCluster_without_iterations(parstart,m,R,q,mu,phi,alpha,beta);
 
  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  // we create the matrix of the expected value z
  become_matrix(n,R,arrayz,z);
 
  // Ask if we have to reparameterizer the phi's (ordered)
  if (*reparC == 1)
  {    
      repar_problem(q,phi);
  }
  //Compute the function to maximize
  *Qstereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	for (r=0; r < *R; r++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_RowCluster_without_iterations(r, j, k, R, m, q, mu, phi, alpha, beta);
	    *Qstereo += (log(theta)*z[i][r]);
	  } 
	} 
      }
    }   
  } 
  *Qstereo = (-1) * (*Qstereo); 
}

void Q_Stereo_RowCluster_rRcC1_inC(
    double *parstart,
    double *arraydata,
    double *arrayz,
    int *n, 
    int *m, 
    int *q, 
    int *R,
    int *reparC,
    double *Qstereo)
{
  double data[*n][*m];
  double z[*n][*R];
  double mu[*q], phi[*q], alpha[*R];
  double theta; 
  int i, j, k, r;
  
  unpack_param_Stereo_RowCluster_rRcC1(parstart,R,q,mu,phi,alpha);
 
  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  // we create the matrix of the expected value z
  become_matrix(n,R,arrayz,z);
 
  // Ask if we have to reparameterizer the phi's (ordered)
  if (*reparC == 1)
  {    
      repar_problem(q,phi);
  }
  //Compute the function to maximize
  *Qstereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	for (r=0; r < *R; r++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_RowCluster_rRcC1(r, k, R, q, mu, phi, alpha);
	    *Qstereo += (log(theta)*z[i][r]);
	  } 
	} 
      }
    }   
  } 
  *Qstereo = (-1) * (*Qstereo); 
}


void Q_Stereo_ColumnCluster_rR1cC_inC(
    double *parstart,
    double *arraydata,
    double *arrayx,
    int *n, 
    int *m, 
    int *q, 
    int *C,
    int *reparC,
    double *Qstereo)
{
  double data[*n][*m];
  double x[*m][*C];
  double mu[*q], phi[*q], beta[*C];
  double theta; 
  int i, j, k, c;
  
  unpack_param_Stereo_ColumnCluster_rR1cC(parstart,C,q,mu,phi,beta);
 
  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  // we create the matrix of the expected value x
  become_matrix(m,C,arrayx,x);
 
  // Ask if we have to reparameterizer the phi's (ordered)
  if (*reparC == 1)
  {    
      repar_problem(q,phi);
  }
  //Compute the function to maximize
  *Qstereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	for (c=0; c < *C; c++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_ColumnCluster_rR1cC(c, k, C, q, mu, phi, beta);
	    *Qstereo += (log(theta)*x[j][c]);
	  } 
	}  
      }
    }   
  } 
  *Qstereo = (-1) * (*Qstereo); 
}


void Q_Stereo_BiCluster_inC(
    double *parstart,
    double *arraydata,
    double *arrayz,
    double *arrayx,
    int *n, 
    int *m, 
    int *q, 
    int *R,
    int *C,
    int *reparC,
    double *Qstereo)
{
  double data[*n][*m];
  double z[*n][*R];
  double x[*m][*C];
  double mu[*q], phi[*q], alpha[*R], beta[*C], arraygamma[(*R * *C)];
  double theta; 
  int i, j, k, r, c;
  
  unpack_param_Stereo_BiCluster(parstart,R,C,q,mu,phi,alpha,beta,arraygamma);
 
  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  // we create the matrix of the expected value z
  become_matrix(n,R,arrayz,z);
  // we create the matrix of the expected value x
  become_matrix(m,C,arrayx,x);
 
  // Ask if we have to reparameterizer the phi's (ordered)
  if (*reparC == 1)
  {    
      repar_problem(q,phi);
  }
  //Compute the function to maximize
  *Qstereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	for (r=0; r < *R; r++)
	{  
	   for (c=0; c < *C; c++)
	   {
	     if (data[i][j] == (k+1))
	     {
	       theta=compute_prob_Stereo_BiCluster(r, c, k, n, q, R, C, mu, phi, alpha, beta, arraygamma);
	       *Qstereo += (log(theta)*z[i][r]*x[j][c]);
	     } 
	   } 
	}
      }
    }   
  } 
  *Qstereo = (-1) * (*Qstereo); 
}


void Q_Stereo_BiCluster_without_iterations_inC(
    double *parstart,
    double *arraydata,
    double *arrayz,
    double *arrayx,
    int *n, 
    int *m, 
    int *q, 
    int *R,
    int *C,
    int *reparC,
    double *Qstereo)
{
  double data[*n][*m];
  double z[*n][*R];
  double x[*m][*C];
  double mu[*q], phi[*q], alpha[*R], beta[*C];
  double theta; 
  int i, j, k, r, c;
  
  unpack_param_Stereo_BiCluster_without_iterations(parstart,R,C,q,mu,phi,alpha,beta);
 
  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  // we create the matrix of the expected value z
  become_matrix(n,R,arrayz,z);
  // we create the matrix of the expected value x
  become_matrix(m,C,arrayx,x);
 
  // Ask if we have to reparameterizer the phi's (ordered)
  if (*reparC == 1)
  {    
      repar_problem(q,phi);
  }
  //Compute the function to maximize
  *Qstereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	for (r=0; r < *R; r++)
	{  
	   for (c=0; c < *C; c++)
	   {
	     if (data[i][j] == (k+1))
	     {
	       theta=compute_prob_Stereo_BiCluster_without_iterations(r, c, k, n, q, mu, phi, alpha, beta);
	       *Qstereo += (log(theta)*z[i][r]*x[j][c]);
	     } 
	   } 
	}
      }
    }   
  } 
  *Qstereo = (-1) * (*Qstereo); 
}





void loglike_Stereo_nullmodel_inC(
    double *parstart,
    double *arraydata,
    int *n, 
    int *m, 
    int *q, 
    int *repar,
    double *loglikestereo)
{
  double data[*n][*m];
  double mu[*q], phi[*q];
  int i, j, k;
  
  unpack_param_Stereo_nullmodel(parstart,q,mu,phi);

  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);

  // Ask if we have to reparameterizer the phi's (ordered)
  if (*repar == 1)
  {    
      repar_problem(q,phi);
  }

  //Compute the function to maximize
  *loglikestereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	if (data[i][j] == (k+1))
	{
	    *loglikestereo += log(compute_prob_Stereo_nullmodel(k, q, mu, phi));
	}        
      }
    } 
  } 
  *loglikestereo = (-1) * (*loglikestereo); 
}

void loglike_Stereo_roweffect_inC(
    double *parstart,
    double *arraydata,
    int *n, 
    int *m, 
    int *q, 
    int *repar,
    double *loglikestereo)
{
  double data[*n][*m];
  double mu[*q], phi[*q], alpha[*n];
  int i, j, k;
   
  unpack_param_Stereo_roweffect(parstart,q,n,mu,phi,alpha);

  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);

  // Ask if we have to reparameterizer the phi's (ordered)
  if (*repar == 1)
  {    
      repar_problem(q,phi);
  }

  //Compute the function to maximize
  *loglikestereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	if (data[i][j] == (k+1))
	{
	    *loglikestereo += log(compute_prob_Stereo_roweffectmodel(i, k, n, q, mu, phi, alpha));
	}        
      }
    } 
  } 
  *loglikestereo = (-1) * (*loglikestereo); 
}


void loglike_Stereo_columneffect_inC(
    double *parstart,
    double *arraydata,
    int *n, 
    int *m, 
    int *q, 
    int *repar,
    double *loglikestereo)
{
  double data[*n][*m];
  double mu[*q], phi[*q], beta[*m];
  int i, j, k;
   
  unpack_param_Stereo_columneffect(parstart,q,m,mu,phi,beta);

  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);

  // Ask if we have to reparameterizer the phi's (ordered)
  if (*repar == 1)
  {    
      repar_problem(q,phi);
  }

  //Compute the function to maximize
  *loglikestereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	if (data[i][j] == (k+1))
	{
	    *loglikestereo += log(compute_prob_Stereo_columneffectmodel(j, k, m, q, mu, phi, beta));
	}        
      }
    } 
  } 
  *loglikestereo = (-1) * (*loglikestereo); 
}

void loglike_Stereo_maineffects_inC(
    double *parstart,
    double *arraydata,
    int *n, 
    int *m, 
    int *q, 
    int *repar,
    double *loglikestereo)
{
  double data[*n][*m];
  double mu[*q], phi[*q], alpha[*n], beta[*m];
  int i, j, k;
   
  unpack_param_Stereo_maineffects(parstart,q,n,m,mu,phi,alpha,beta);

  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);

  // Ask if we have to reparameterizer the phi's (ordered)
  if (*repar == 1)
  {    
      repar_problem(q,phi);
  }

  //Compute the function to maximize
  *loglikestereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	if (data[i][j] == (k+1))
	{
	  *loglikestereo += log(compute_prob_Stereo_maineffectsmodel(i, j, k, n, m, q, mu, phi, alpha, beta));  
	}        
      }
    } 
  } 
  
  *loglikestereo = (-1) * (*loglikestereo);   
}




//Build the logL incomplete by sum over all the the columns (kappa's)
//This function is used when we optimise the COMPLETE logL.
void logL_INCOMP_BiClusterNOiter_sumcols(
    double *parhatpack,
    double *arraydata,
    int *n, 
    int *m, 
    int *q,
    int *R,
    int *C,
    double *arraykappamat,
    double *arraybetamat,
    int *totalcombs,
    double *arraycombosmat,
    int *nrowscombosmat,
    double *mu,
    double *phi,
    double *alpha,
    double *beta,
    double *resultcombv,
    double *thetav,
    int *reparC,
    double *logLincompletehat)
{
  double data[*n][*m], betamat[*totalcombs][*m];
  double combosmat[*nrowscombosmat][*m];
  double kappaC[*C], piR[*R];
  double secondterm, theta,
         thetabypi, thetaacum, sumathetav, valuetheta, aux; 
  int i, r, j, k, l, numcomb, indexc;
   
  // First, we create the alpha and beta arrays
  unpack_param_Stereo_BiCluster_without_iterations_withpiR_and_kappaC(parhatpack,R,C,q,mu,phi,alpha,beta,piR,kappaC);
 // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  become_matrix(totalcombs,m,arraybetamat,betamat);
  become_matrix(totalcombs,m,arraycombosmat,combosmat); //nrowscombosmat

  if (*reparC == 1)
  {
    recover_mix_from_S(R,piR); //pi1=expit(s1), pi2=expit(s2)*(1-expit(s1)),...
    recover_mix_from_S(C,kappaC); 
  }  
  for (numcomb=0; numcomb<*totalcombs; numcomb++)
  {    
    secondterm = 1;
    for (i=0; i < *n; i++)
    {
      thetabypi = 0;
      for (r=0; r < *R; r++)
      {
	thetaacum = 1;
	for (j=0; j < *m; j++)
	{
	  for (k=0; k < *q; k++)
	  {
	    theta = 1;
	    if (data[i][j] == (k+1))
	    {
	      sumathetav = 0;
              for (l=0; l < *q; l++) 
              {
		indexc = (int)((combosmat[numcomb][j])-1);
		valuetheta = mu[l]+(phi[l]*(alpha[r]+beta[indexc]));		    
 //      valuetheta = mu[l]+(phi[l]*(alpha[r]+betamat[numcomb][j]));          
		thetav[l] = exp(valuetheta);
		sumathetav += thetav[l];
              }
              theta = thetav[k]/sumathetav; 
	      thetaacum *= theta;
	    }   
	  }
	}
	thetabypi = thetabypi + (thetaacum * piR[r]);	
      }  
      secondterm *= thetabypi;
    }
    
    for (i=0; i<*totalcombs; i++)
    { 
	resultcombv[i] = secondterm;
       	for (j=0; j < *m; j++)
	{
	   indexc = (int)((combosmat[i][j])-1);
	   for (l=0; l < *C; l++)
	   {
	     if (indexc == l) resultcombv[i] *= kappaC[l];
	   }  
	}   
    }	
  }

  aux = 0;
  for (numcomb = 0; numcomb < *totalcombs; numcomb++)
  {
    aux = (aux + resultcombv[numcomb]);
  } 
  *logLincompletehat = log(aux); 

}  


//Build the logL incomplete by sum over all the the rows (pi's)
//This function is used when we optimise the COMPLETE logL.
void logL_INCOMP_BiClusterNOiter_sumrows(
    double *parhatpack,
    double *arraydata,
    int *n, 
    int *m, 
    int *q,
    int *R,
    int *C,
    double *arraypimat,
    double *arrayalphamat,
    int *totalcombs,
    double *arraycombosmat,
    int *nrowscombosmat,
    double *mu,
    double *phi,
    double *alpha,
    double *beta,
    double *resultcombv,
    double *thetav,
    int *reparC,
    double *logLincompletehat)
{
  double data[*n][*m], alphamat[*totalcombs][*n]; 
  double combosmat[*nrowscombosmat][*n];
  double kappaC[*C], piR[*R];
  double secondterm, theta,
         thetabykappa, thetaacum, sumathetav, valuetheta, aux; 
  int i, c, j, k, l, numcomb,  indexr;
 
  // First, w  POR AQUIe create the alpha and beta arrays
  unpack_param_Stereo_BiCluster_without_iterations_withpiR_and_kappaC(parhatpack,R,C,q,mu,phi,alpha,beta,piR,kappaC);
 
 // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  become_matrix(totalcombs,n,arrayalphamat,alphamat);
  
  become_matrix(totalcombs,n,arraycombosmat,combosmat); //nrowscombosmat

  if (*reparC == 1)
  {
    recover_mix_from_S(R,piR);
    recover_mix_from_S(C,kappaC); //kappa1=expit(s1), kappa2=expit(s2)*(1-expit(s1)),...
  }  
  
  for (numcomb=0; numcomb<*totalcombs; numcomb++)
  {  
    secondterm = 1;
    for (j=0; j < *m; j++)
    {
      thetabykappa = 0;
      for (c=0; c < *C; c++)
      {
	thetaacum = 1;
	for (i=0; i < *n; i++)
	{
	  for (k=0; k < *q; k++)
	  {
	    theta = 1;
	    if (data[i][j] == (k+1))
	    { 
	      sumathetav=0;
              for (l=0; l < *q; l++) 
              {
		indexr = (int)((combosmat[numcomb][i])-1);
		valuetheta = mu[l]+(phi[l]*(alpha[indexr]+beta[c]));		    
 //		valuetheta=mu[l]+(phi[l]*(alphamat[numcomb][i]+beta[c]));
		thetav[l] = exp(valuetheta);
		sumathetav += thetav[l];
              }
              theta = thetav[k]/sumathetav; 
              thetaacum *= theta;
	    }   
	  }
	}
	thetabykappa += (thetaacum * kappaC[c]);	
      }  
      secondterm *= thetabykappa;
    }

    for (i=0; i<*totalcombs; i++)
    { 
	resultcombv[i] = secondterm;
       	for (j=0; j < *n; j++)
	{
	   indexr = (int)((combosmat[i][j])-1);
	   for (l=0; l < *R; l++)
	   {
	     if (indexr == l) resultcombv[i] *= piR[l];
	   }  
	}   
    }	
  }  
  aux = 0;
  for (numcomb = 0; numcomb < *totalcombs; numcomb++)
  {
    aux += resultcombv[numcomb];
  } 
  *logLincompletehat = log(aux); 
}  

// void logLikelihood_Incomplete_Stereo_BiCluster_without_iterations_sumcols(
//     double *parhatpack,
//     double *piR,
//     double *arraydata,
//     int *n, 
//     int *m, 
//     int *q,
//     int *R,
//     int *C,
//     double *arraykappamat,
//     double *arraybetamat,
//     int *totalcombs,
//     double *kappaprodv,
//     double *mu,
//     double *phi,
//     double *alpha,
//     double *beta,
//     double *resultcombv,
//     double *thetav,
//     double *logLincompletehat)
// {
//   double data[*n][*m], betamat[*totalcombs][*m]; 
//   double secondterm, theta, logLikelihoodincompletestereo, thetabypi, thetaacum, sumathetav, valuetheta, aux; 
//   int i, r, c, j, k, l, numcomb, index;
//    
//   // First, we create the alpha and beta arrays
//   unpack_param_Stereo_BiCluster_without_iterations(parhatpack,R,C,q,mu,phi,alpha,beta);
//  
//  // we create the matrix of our data
//   become_matrix(n,m,arraydata,data);
//   become_matrix(totalcombs,m,arraybetamat,betamat);
// 
//   for (numcomb=0; numcomb<*totalcombs; numcomb++)
//   {  
//     secondterm = 1;
//     for (i=0; i < *n; i++)
//     {
//       thetabypi = 0;
//       for (r=0; r < *R; r++)
//       {
// 	thetaacum = 1;
// 	for (j=0; j < *m; j++)
// 	{
// 	  for (k=0; k < *q; k++)
// 	  {
// 	    theta = 1;
// 	    if (data[i][j] == (k+1))
// 	    {
// 	      sumathetav = 0;
//               for (l=0; l < *q; l++) 
//               {
//                 valuetheta = mu[l]+(phi[l]*(alpha[r]+betamat[numcomb][j])); 
// 		thetav[l] = exp(valuetheta);
// 		sumathetav += thetav[l];
//               }
//               theta = thetav[k]/sumathetav; 
// 	      thetaacum *= theta;
// 	    }   
// 	  }
// 	}
// 	thetabypi = thetabypi + (thetaacum * piR[r]);	
//       }  
//       secondterm *= thetabypi;
//     }
//     resultcombv[numcomb] = kappaprodv[numcomb]*secondterm;
//   }
// 
//   aux = 0;
//   for (numcomb = 0; numcomb < *totalcombs; numcomb++)
//   {
//     aux = (aux + resultcombv[numcomb]);
//   } 
//   *logLincompletehat = log(aux);
// }  


// void logLikelihood_Incomplete_Stereo_BiCluster_without_iterations_sumrows(
//     double *parhatpack,
//     double *kappaC,
//     double *arraydata,
//     int *n, 
//     int *m, 
//     int *q,
//     int *R,
//     int *C,
//     double *arraypimat,
//     double *arrayalphamat,
//     int *totalcombs,
//     double *piprodv,
//     double *mu,
//     double *phi,
//     double *alpha,
//     double *beta,
//     double *resultcombv,
//     double *thetav,
//     double *logLincompletehat)
// {
//   double data[*n][*m], alphamat[*totalcombs][*n]; 
//   double secondterm, theta, logLikelihoodincompletestereo, thetabykappa, thetaacum, sumathetav, valuetheta, aux; 
//   int i, r, c, j, k, l, numcomb, index;
//    
// 
//   unpack_param_Stereo_BiCluster_without_iterations(parhatpack,R,C,q,mu,phi,alpha,beta);
//  
//  // we create the matrix of our data
//   become_matrix(n,m,arraydata,data);
//   become_matrix(totalcombs,n,arrayalphamat,alphamat);
// 
// 
//   for (numcomb=0; numcomb<*totalcombs; numcomb++)
//   {  
//     secondterm = 1;
//     for (j=0; j < *m; j++)
//     {
//       thetabykappa = 0;
//       for (c=0; c < *C; c++)
//       {
// 	thetaacum = 1;
// 	for (i=0; i < *n; i++)
// 	{
// 	  for (k=0; k < *q; k++)
// 	  {
// 	    theta = 1;
// 	    if (data[i][j] == (k+1))
// 	    { 
// 	      sumathetav=0;
//               for (l=0; l < *q; l++) 
//               {
// 		valuetheta=mu[l]+(phi[l]*(alphamat[numcomb][i]+beta[c]));
// 		thetav[l] = exp(valuetheta);
// 		sumathetav += thetav[l];
//               }
//               theta = thetav[k]/sumathetav; 
//               thetaacum *= theta;
// 	    }   
// 	  }
// 	}
// 	thetabykappa += (thetaacum * kappaC[c]);	
//       }  
//       secondterm *= thetabykappa;
//     }
//     resultcombv[numcomb] = piprodv[numcomb]*secondterm;
//   }
// 
//   aux = 0;
//   for (numcomb = 0; numcomb < *totalcombs; numcomb++)
//   {
//     aux += resultcombv[numcomb];
//   } 
//   *logLincompletehat = log(aux); 
// }  

//Build the logL incomplete by sum over all the the columns (kappa's)
//This function is used when we optimise the COMPLETE logL.
void logL_INCOMP_BiCluster_sumcols(
    double *parhatpack,
    double *arraydata,
    int *n, 
    int *m, 
    int *q,
    int *R,
    int *C,
    double *arraykappamat,
    double *arraybetamat,
    double *arraygammamat,
    int *totalcombs,
    double *arraycombosmat,
    int *nrowscombosmat,
    double *mu,
    double *phi,
    double *alpha,
    double *beta,
    double *resultcombv,
    double *thetav,
    int *reparC,
    double *logLincompletehat)
{ 
  double data[*n][*m], betamat[*totalcombs][*m], gamma[*R][*C];
  double combosmat[*nrowscombosmat][*m];
  double arraygamma[(*R * *C)];
  double kappaC[*C], piR[*R];
  double secondterm, theta, thetabypi,
         thetaacum, sumathetav, valuetheta, aux; 
  int i, r, j, k, l, numcomb, indexc;

  // First, we create the alpha and beta arrays
  unpack_param_Stereo_BiCluster_withpiR_and_kappaC(parhatpack,R,C,q,mu,phi,alpha,beta,arraygamma,piR,kappaC);

 // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  become_matrix(totalcombs,m,arraybetamat,betamat);
  
  become_matrix(R,C,arraygamma,gamma); 
  become_matrix(totalcombs,m,arraycombosmat,combosmat); //nrowscombosmat
  
  if (*reparC == 1)
  {
    recover_mix_from_S(R,piR); //pi1=expit(s1), pi2=expit(s2)*(1-expit(s1)),...
    recover_mix_from_S(C,kappaC);  
  }  
  
  for (numcomb=0; numcomb<*totalcombs; numcomb++)
  {    
    secondterm = 1;
    for (i=0; i < *n; i++)
    {
      thetabypi = 0;
      for (r=0; r < *R; r++)
      {
	thetaacum = 1;
	for (j=0; j < *m; j++)
	{
	  for (k=0; k < *q; k++)
	  {
	    theta = 1;
	    if (data[i][j] == (k+1))
	    {
	      sumathetav = 0;
              for (l=0; l < *q; l++) 
              {
		indexc = (int)((combosmat[numcomb][j])-1);
		valuetheta = mu[l]+(phi[l]*(alpha[r]+beta[indexc]+gamma[r][indexc]));		    
//                valuetheta=mu[l]+(phi[l]*(alpha[r]+betamat[numcomb][j]+arraygammamat[(*totalcombs * *R * j)+(*totalcombs * r) + numcomb]));
		thetav[l] = exp(valuetheta);
		sumathetav += thetav[l];
              }
              theta = thetav[k]/sumathetav; 
	      thetaacum *= theta;
	    }   
	  }
	}
	thetabypi = thetabypi + (thetaacum * piR[r]);	
      }  
      secondterm *= thetabypi;
    }
    for (i=0; i<*totalcombs; i++)
    { 
	resultcombv[i] = secondterm;
       	for (j=0; j < *m; j++)
	{
	   indexc = (int)((combosmat[i][j])-1);
	   for (l=0; l < *C; l++)
	   {
	     if (indexc == l) resultcombv[i] *= kappaC[l];
	   }  
	}   
    }	
  }
  	

  aux = 0;
  for (numcomb = 0; numcomb < *totalcombs; numcomb++)
  {
    aux = (aux + resultcombv[numcomb]);
  } 
  *logLincompletehat = log(aux);
 
}  


// void logLikelihood_Incomplete_Stereo_BiCluster_sumcols(
//     double *parstart,
//     double *piR,
//     double *arraydata,
//     int *n, 
//     int *m, 
//     int *q,
//     int *R,
//     int *C,
//     double *arraykappamat,
//     double *arraybetamat,
//     double *arraygammamat,
//     int *totalcombs,
//     double *kappaprodv,
//     double *mu,
//     double *phi,
//     double *alpha,
//     double *beta,
//     double *resultcombv,
//     double *thetav,
//     double *logLincompletehat)
// {
//   double data[*n][*m], betamat[*totalcombs][*m]; 
//   double arraygamma[(*R * *C)];
//   double secondterm, theta, logLikelihoodincompletestereo, thetabypi, thetaacum, sumathetav, valuetheta, aux; 
//   int i, j, k, r, c, l, numcomb, index;
//    
//   // First, we create the alpha and beta arrays
//   unpack_param_Stereo_BiCluster(parstart,R,C,q,mu,phi,alpha,beta,arraygamma);
//   // We create the matrix of our data
//   become_matrix(n,m,arraydata,data);
//   become_matrix(totalcombs,m,arraybetamat,betamat);
// 
//   for (numcomb=0; numcomb<*totalcombs; numcomb++)
//   {  
//     secondterm = 1;
//     for (i=0; i < *n; i++)
//     {
//       thetabypi = 0;
//       for (r=0; r < *R; r++)
//       {
// 	thetaacum = 1;
// 	for (j=0; j < *m; j++)
// 	{
// 	  for (k=0; k < *q; k++)
// 	  {
// 	    theta = 1;
// 	    if (data[i][j] == (k+1))
// 	    {
// 	     sumathetav = 0;
//              for (l=0; l < *q; l++) 
//              {
// 		valuetheta=mu[l]+(phi[l]*(alpha[r]+betamat[numcomb][j]+arraygammamat[(*totalcombs * *R * j)+(*totalcombs * r) + numcomb]));
// 		thetav[l] = exp(valuetheta);
// 		sumathetav += thetav[l];
//               }
//               theta = thetav[k]/sumathetav; 
//               thetaacum *= theta;
// 	    }   
// 	  }
// 	}
// 	thetabypi = thetabypi + (thetaacum * piR[r]);
//       }  
//       secondterm *= thetabypi;
//     }
//     resultcombv[numcomb] = kappaprodv[numcomb]*secondterm;
//   }
// 
//   aux = 0;
//   for (numcomb = 0; numcomb < *totalcombs; numcomb++)
//   {
//     aux = (aux + resultcombv[numcomb]);
//   } 
//   *logLincompletehat = log(aux);
// }   

//Build the logL incomplete by sum over all the the rows (pi's)
//This function is used when we optimise the COMPLETE logL.
void logL_INCOMP_BiCluster_sumrows(
    double *parhatpack,
    double *arraydata,
    int *n, 
    int *m, 
    int *q,
    int *R,
    int *C,
    double *arraypimat,
    double *arrayalphamat,
    double *arraygammamat,
    int *totalcombs,
    double *arraycombosmat,
    int *nrowscombosmat,
    double *mu,
    double *phi,
    double *alpha,
    double *beta,
    double *resultcombv,
    double *thetav,
    int *reparC,
    double *logLincompletehat)
{
  double data[*n][*m], alphamat[*totalcombs][*n], gamma[*R][*C];
  double combosmat[*nrowscombosmat][*n];
  double arraygamma[(*R * *C)];
  double kappaC[*C], piR[*R];
  double secondterm, theta, thetabykappa,
  thetaacum, sumathetav, valuetheta, aux; 
  int i, c, j, k, l, numcomb,  indexr;
   
  // First, we create the alpha and beta arrays
  unpack_param_Stereo_BiCluster_withpiR_and_kappaC(parhatpack,R,C,q,mu,phi,alpha,beta,arraygamma,piR,kappaC);
 // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  become_matrix(totalcombs,n,arrayalphamat,alphamat);
  
  become_matrix(R,C,arraygamma,gamma); 
  become_matrix(totalcombs,n,arraycombosmat,combosmat); //nrowscombosmat

  if (*reparC == 1)
  {
    recover_mix_from_S(R,piR); 
    recover_mix_from_S(C,kappaC); //kappa1=expit(s1), kappa2=expit(s2)*(1-expit(s1)),...
  }  
  
  for (numcomb=0; numcomb<*totalcombs; numcomb++)
  {  
    secondterm = 1;
    for (j=0; j < *m; j++)
    {
      thetabykappa = 0;
      for (c=0; c < *C; c++)
      {
	thetaacum = 1;
	for (i=0; i < *n; i++)
	{
	  for (k=0; k < *q; k++)
	  {
	    theta = 1;
	    if (data[i][j] == (k+1))
	    { 
	      sumathetav=0;
              for (l=0; l < *q; l++) 
              {
		indexr = (int)((combosmat[numcomb][i])-1);
		valuetheta = mu[l]+(phi[l]*(alpha[indexr]+beta[c]+gamma[indexr][c]));	
//		valuetheta=mu[l]+(phi[l]*(alphamat[numcomb][i]+beta[c]+arraygammamat[(*totalcombs * *C * i)+(*totalcombs * c) + numcomb]));
		thetav[l] = exp(valuetheta);
		sumathetav += thetav[l];
              }
              theta = thetav[k]/sumathetav; 
              thetaacum *= theta;
	    }   
	  }
	}
	thetabykappa += (thetaacum * kappaC[c]);	
      }  
      secondterm *= thetabykappa;
    }
    for (i=0; i<*totalcombs; i++)
    { 
	resultcombv[i] = secondterm;
       	for (j=0; j < *n; j++)
	{
	   indexr = (int)((combosmat[i][j])-1);
	   for (l=0; l < *R; l++)
	   {
	     if (indexr == l) resultcombv[i] *= kappaC[l];
	   }  
	}   
    }	
  }

  aux = 0;
  for (numcomb = 0; numcomb < *totalcombs; numcomb++)
  {
    aux += resultcombv[numcomb];
  } 
  *logLincompletehat = log(aux); 
}  


// void logLikelihood_Incomplete_Stereo_BiCluster_sumrows(
//     double *parstart,
//     double *kappaC,
//     double *arraydata,
//     int *n, 
//     int *m, 
//     int *q,
//     int *R,
//     int *C,
//     double *arraypimat,
//     double *arrayalphamat,
//     double *arraygammamat,
//     int *totalcombs,
//     double *piprodv,
//     double *mu,
//     double *phi,
//     double *alpha,
//     double *beta,
//     double *resultcombv,
//     double *thetav,
//     double *logLincompletehat)
// {
//   double data[*n][*m], alphamat[*totalcombs][*n]; 
//   double arraygamma[(*R * *C)];
//   double secondterm, theta, logLikelihoodincompletestereo, thetabykappa, thetaacum, sumathetav, valuetheta, aux; 
//   int i, r, c,  j, k, l, index, numcomb;
//    
//   // First, we create the alpha and beta arrays
//   unpack_param_Stereo_BiCluster(parstart,R,C,q,mu,phi,alpha,beta,arraygamma);
// 
//   //We create the matrix of our data
//   become_matrix(n,m,arraydata,data);
//   become_matrix(totalcombs,n,arrayalphamat,alphamat);
// 
//   for (numcomb=0; numcomb<*totalcombs; numcomb++)
//   {  
//     secondterm = 1;
//     for (j=0; j < *m; j++)
//     {
//       thetabykappa = 0;
//       for (c=0; c < *C; c++)
//       {
// 	thetaacum = 1;
// 	for (i=0; i < *n; i++)
// 	{
// 	  for (k=0; k < *q; k++)
// 	  {
// 	    theta = 1;
// 	    if (data[i][j] == (k+1))
// 	    {
// 	     sumathetav=0;
//              for (l=0; l < *q; l++) 
//              {
//                valuetheta=mu[l]+(phi[l]*(alphamat[numcomb][i]+beta[c]+arraygammamat[(*totalcombs * *C * i)+(*totalcombs * c) + numcomb]));
// 	       thetav[l] = exp(valuetheta);
// 	       sumathetav += thetav[l];
//               }
//               theta = thetav[k]/sumathetav; 
//               thetaacum *= theta;
// 	    }   
// 	  }
// 	}
// 	thetabykappa += (thetaacum * kappaC[c]);
//       }  
//       secondterm *= thetabykappa;
//     }
//     resultcombv[numcomb] = piprodv[numcomb]*secondterm;
//   }
// 
//   aux = 0;
//   for (numcomb = 0; numcomb < *totalcombs; numcomb++)
//   {
//     aux = (aux + resultcombv[numcomb]);
//   } 
//   *logLincompletehat = log(aux); 
// }   


void logLikelihood_Incomplete_Stereo_RowClustering_rRcC1_inC(
    double *parstart,
    double *arraydata,
    int *n, 
    int *m, 
    int *q, 
    int *R,
    int *reparC,
    double *logL)
{
  double data[*n][*m];
  double mu[*q], phi[*q], alpha[*R], piR[*R];
  double theta, thetabypi, thetaacum; 
  int i, j, k, r;
  
  unpack_param_Stereo_RowCluster_rRcC1_withpiR(parstart,R,q,mu,phi,alpha,piR);

  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
 
  // Ask if we have to reparameterizer
  // input: S  output: piR
  if (*reparC == 1)
  {    
      recover_mix_from_S(R,piR);
  }

  //Compute the function to maximize
  *logL = 0.0;
  for (i=0; i < *n; i++)
  {
    thetabypi = 0;
    for (r=0; r < *R; r++)
    {
      thetaacum = 1;
      for (j=0; j < *m; j++)
      {
	for (k=0; k < *q; k++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_RowCluster_rRcC1(r, k, R, q, mu, phi, alpha);
	    thetaacum *= theta;
	  } 
	} 
      }	
      thetabypi += (thetaacum * piR[r]); 
    }   
    *logL += log(thetabypi);  
  } 
   *logL = (-1) * (*logL); 
  
}



void logLikelihood_Incomplete_Stereo_RowClustering_rRcm_without_iterations_inC(
    double *parstart,
    double *arraydata,
    int *n, 
    int *m, 
    int *q, 
    int *R,
    int *reparC,
    double *logL)
{
  double data[*n][*m];
  double mu[*q], phi[*q], alpha[*R], beta[*m], piR[*R];
  double theta, thetabypi, thetaacum; 
  int i, j, k, r;
  
  unpack_param_Stereo_RowCluster_without_iterations_withpiR(parstart,m,R,q,mu,phi,alpha,beta,piR);

  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
 
  // Ask if we have to reparameterizer
  // input: S  output: piR
  if (*reparC == 1)
  {    
      recover_mix_from_S(R,piR);
  }

  //Compute the function to maximize
  *logL = 0.0;
  for (i=0; i < *n; i++)
  {
    thetabypi = 0;
    for (r=0; r < *R; r++)
    {
      thetaacum = 1;
      for (j=0; j < *m; j++)
      {
	for (k=0; k < *q; k++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_RowCluster_without_iterations(r, j, k, R, m, q, mu, phi, alpha, beta);
	    thetaacum *= theta;
	  } 
	} 
      }
      thetabypi += (thetaacum * piR[r]); 
    }   
    *logL += log(thetabypi);  
  } 
   *logL = (-1) * (*logL); 
  
}


void logLikelihood_Incomplete_Stereo_RowClustering_rRcm_inC(
    double *parstart,
    double *arraydata,
    int *n, 
    int *m, 
    int *q, 
    int *R,
    int *reparC,
    double *logL)
{
  double data[*n][*m];
  double mu[*q], phi[*q], alpha[*R], beta[*m], piR[*R];
  double arraygamma[*R * *m];
  double theta, thetabypi, thetaacum; 
  int i, j, k, r;
  
  
  unpack_param_Stereo_RowCluster_withpiR(parstart,m,R,q,mu,phi,alpha,beta,arraygamma,piR);

  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
 
  // Ask if we have to reparameterizer
  // input: S  output: piR
  if (*reparC == 1)
  {    
      recover_mix_from_S(R,piR);
  }

  //Compute the function to maximize
  *logL = 0.0;
  for (i=0; i < *n; i++)
  {
    thetabypi = 0;
    for (r=0; r < *R; r++)
    {
      thetaacum = 1;
      for (j=0; j < *m; j++)
      {
	for (k=0; k < *q; k++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_RowCluster(r, j, k, R, m, q, mu, phi, alpha, beta,arraygamma);
	    thetaacum *= theta;
	  } 
	} 
      }
      thetabypi += (thetaacum * piR[r]); 
    }   
    *logL += log(thetabypi);  
  } 
   *logL = (-1) * (*logL); 
  
}

void logLikelihood_Incomplete_Stereo_ColClustering_rR1cC_inC(
    double *parstart,
    double *arraydata,
    int *n, 
    int *m, 
    int *q, 
    int *C,
    int *reparC,
    double *logL)
{
  double data[*n][*m];
  double mu[*q], phi[*q], beta[*C], kappaC[*C];
  double theta, thetabykappa, thetaacum; 
  int i, j, k, c;
  
  unpack_param_Stereo_ColumnCluster_rR1cC_withkappaC(parstart,C,q,mu,phi,beta,kappaC);

  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
 
  // Ask if we have to reparameterizer
  // input: S  output: piR
  if (*reparC == 1)
  {    
      recover_mix_from_S(C,kappaC);
  }

  //Compute the function to maximize
  *logL = 0.0;
  for (i=0; i < *n; i++)
  {
    thetabykappa = 0;
    for (c=0; c < *C; c++)
    {
      thetaacum = 1;
      for (j=0; j < *m; j++)
      {
	for (k=0; k < *q; k++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_ColumnCluster_rR1cC(c, k, C, q, mu, phi, beta);
	    thetaacum *= theta;
	  } 
	} 
      }
      thetabykappa += (thetaacum * kappaC[c]); 
    }   
    *logL += log(thetabykappa);  
  } 
   *logL = (-1) * (*logL); 
  
}



void logLikelihood_Incomplete_Stereo_ColClustering_rncC_without_iterations_inC(
    double *parstart,
    double *arraydata,
    int *n, 
    int *m, 
    int *q, 
    int *C,
    int *reparC,
    double *logL)
{
  double data[*n][*m];
  double mu[*q], phi[*q], alpha[*n], beta[*C], kappaC[*C];
  double theta, thetabykappa, thetaacum; 
  int i, j, k, c;
  
  unpack_param_Stereo_ColumnCluster_without_iterations_withkappaC(parstart,n,C,q,mu,phi,alpha,beta,kappaC);

  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
 
  // Ask if we have to reparameterizer
  // input: S  output: piR
  if (*reparC == 1)
  {    
      recover_mix_from_S(C,kappaC);
  }

  //Compute the function to maximize
  *logL = 0.0;
  for (j=0; j < *m; j++)
  {
    thetabykappa = 0;
    for (c=0; c < *C; c++)
    {
      thetaacum = 1;
      for (i=0; i < *n; i++)
      {
	for (k=0; k < *q; k++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_ColumnCluster_without_iterations(c, i, k, C, n, q, mu, phi, alpha, beta);
	    thetaacum *= theta;
	  } 
	} 
      }
      thetabykappa += (thetaacum * kappaC[c]); 
    }   
    *logL += log(thetabykappa);  
  } 
   *logL = (-1) * (*logL); 
  
}

void logLikelihood_Incomplete_Stereo_ColClustering_rncC_inC(
    double *parstart,
    double *arraydata,
    int *n, 
    int *m, 
    int *q, 
    int *C,
    int *reparC,
    double *logL)
{
  double data[*n][*m];
  double mu[*q], phi[*q], alpha[*n], beta[*C], kappaC[*C];
  double arraygamma[*n * *C];
  double theta, thetabykappa, thetaacum; 
  int i, j, k, c;
    
  unpack_param_Stereo_ColumnCluster_withkappaC(parstart,n,C,q,mu,phi,alpha,beta,arraygamma,kappaC);

  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
 
  // Ask if we have to reparameterizer
  // input: S  output: piR
  if (*reparC == 1)
  {    
      recover_mix_from_S(C,kappaC);
  }

  //Compute the function to maximize
  *logL = 0.0;
  for (j=0; j < *m; j++)
  {
    thetabykappa = 0;
    for (c=0; c < *C; c++)
    {
      thetaacum = 1;
      for (i=0; i < *n; i++)
      {
	for (k=0; k < *q; k++)
	{
	  if (data[i][j] == (k+1))
	  {
	    theta=compute_prob_Stereo_ColumnCluster(c, i, k, C, n, q, mu, phi, alpha, beta, arraygamma);
	    thetaacum *= theta;
	  } 
	} 
      }
      thetabykappa += (thetaacum * kappaC[c]); 
    }   
    *logL += log(thetabykappa);  
  } 
  *logL = (-1) * (*logL); 
}


//**********************************************************
//**** FUNCTIONS FOR POM ***********************************
//**********************************************************


void unpack_param_BiCluster_without_iterations_POM(
  double *parstart,
  int	 *R,
  int 	 *C,
  int	 *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta)
{
  int r, c, k, index;
  
  index=0;
 
  for (k=0; k<(*q-1); k++)
  {
    mu[k] = parstart[index];
    index++;
  } 
 
  alpha[0] = 0;
  for (r=1; r<*R; r++)
  {
    alpha[r] = parstart[index];
    index++;
  }


  beta[0]=0;
  for (c=1; c<*C; c++)
  {
     beta[c] = parstart[index];
     index++;
  } 

  
}

//this is no a "external" function

	   
double compute_prob_BiCluster_without_iterations_POM( 
  int r,
  int c,
  int k,
  int *n,
  int *q,
  double *mu,
  double *phi,
  double *alpha,
  double *beta)
{
  double theta;
  double thetanum, thetaden, thetanum2, thetaden2;
  

  if (k == 0)
  {    
    thetanum=exp(mu[k]-alpha[r]-beta[c]);
    thetaden=1+thetanum;
    
    theta=thetanum/thetaden;
  }
     
  if ( (k > 0) && (k < (*q-1)) )
  {    
    thetanum=exp(mu[k]-alpha[r]-beta[c]);
    thetaden=1+thetanum;
    
    thetanum2=exp(mu[k-1]-alpha[r]-beta[c]);
    thetaden2=1+thetanum2;
    
    
    theta=(thetanum/thetaden)-(thetanum2/thetaden2);
  }   
     
  if (k == (*q-1))
  {    
      thetanum=exp(mu[k-1]-alpha[r]-beta[c]);
      thetaden=1+thetanum;
      
      theta=1-(thetanum/thetaden);
  }
 
 
  return(theta);
}


void Q_Stereo_BiCluster_without_iterations_inC_POM(
    double *parstart,
    double *arraydata,
    double *arrayz,
    double *arrayx,
    int *n, 
    int *m, 
    int *q, 
    int *R,
    int *C,
    int *reparC,
    double *Qstereo)
{
  double data[*n][*m];
  double z[*n][*R];
  double x[*m][*C];
  double mu[*q], phi[*q], alpha[*R], beta[*C];
  double theta; 
  int i, j, k, r, c;
  
  unpack_param_BiCluster_without_iterations_POM(parstart,R,C,q,mu,phi,alpha,beta);
 
  // we create the matrix of our data
  become_matrix(n,m,arraydata,data);
  // we create the matrix of the expected value z
  become_matrix(n,R,arrayz,z);
  // we create the matrix of the expected value x
  become_matrix(m,C,arrayx,x);
 
  // Ask if we have to reparameterizer the phi's (ordered)
  if (*reparC == 1)
  {    
      repar_problem(q,phi);
  }
  //Compute the function to maximize
  *Qstereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	for (r=0; r < *R; r++)
	{  
	   for (c=0; c < *C; c++)
	   {
	     if (data[i][j] == (k+1))
	     {
	       theta=compute_prob_BiCluster_without_iterations_POM(r, c, k, n, q, mu, phi, alpha, beta);
	       *Qstereo += (log(theta)*z[i][r]*x[j][c]);
	     } 
	   } 
	}
      }
    }   
  } 
  *Qstereo = (-1) * (*Qstereo); 
}

///////////////////////
// FAR functions //////
///////////////////////



void become_matrix_3D(
  int *n,
  int *m,
  int *S,
  double *arraydata,
  double data[*n][*m][*S])
{
  int i, j, k;
  
  for (i=0; i<*n; i++)
  {
    for (j=0; j<*m; j++)
    {
      for (k=0; k<*S; k++)
      {	
        data[i][j][k] = arraydata[(*n * *m * k) + (*n * j) + i];
      }	
    }
  }
}


void Q_Stereo_RowCluster_without_iterations_inC_FAR(
    double *parstart,
    double *arraydata,
    double *arrayz,
    int *n, 
    int *m, 
    int *S,
    int *q, 
    int *R,
    int *reparC,
    double *Qstereo)
{
  double data[*n][*m][*S];
  double z[*n][*R];
  double mu[*q], phi[*q], alpha[*R], beta[*m];
  double theta; 
  int i, j, k, t, r;
  
  unpack_param_Stereo_RowCluster_without_iterations(parstart,m,R,q,mu,phi,alpha,beta);
 
  // we create the matrix of our data
  become_matrix_3D(n,m,S,arraydata,data);
  // we create the matrix of the expected value z
  become_matrix(n,R,arrayz,z);
 
  // Ask if we have to reparameterizer the phi's (ordered)
  if (*reparC == 1)
  {    
      repar_problem(q,phi);
  }
  //Compute the function to maximize
  *Qstereo = 0.0;
  for (i=0; i < *n; i++)
  {
    for (j=0; j < *m; j++)
    {
      for (k=0; k < *q; k++)
      {
	for (t=0; t < *S; t++)
	{  
	  for (r=0; r < *R; r++)
	  {
	    if (data[i][j][t] == (k+1))
	    {
	      theta=compute_prob_Stereo_RowCluster_without_iterations(r, j, k, R, m, q, mu, phi, alpha, beta);
	      *Qstereo += (log(theta)*z[i][r]);
	    } 
	  }
	}    
      }
    }   
  } 
  *Qstereo = (-1) * (*Qstereo); 

}

