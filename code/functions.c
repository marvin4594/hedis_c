#include <math.h>



double sign(double a)
{
  return (a>0)?1:-1; 
}





double min(double a, double b)
{
  return (a<b)?a:b; 
}





double max(double a, double b)
{
  return (a>b)?a:b; 
}





double minmod(double a, double b)
{

  double value;
  
  if (fabs(a) < fabs(b))
  {
    value = a;
  }
  else
  {
    value = b;
  }

  if (a*b < 0.0) value = 0.0;

  return value;
  
} 





double H2F1(double z)
{
  //computes the hypergeometrical function 2F1

  int n;
  double a, b, c, R, result;


  n = 0;
  a = 0.5;
  b = -1.0;
  c = 0.69314718-0.25*log(1.0-z);
  result = 0.0;


  if (z > 0.5)
  {

    do
    {

      result += pow(a,2)*(b+c)*pow((1.0-z),n);
      n++;
      a *= (1.0 + 1.0/(2.0*n));
      b += 1.0/(2.0*n) - 1.0/(2.0*n+1.0);
      R = fabs((pow(a,2)*b*pow((1.0-z),n))/(1.0-pow((1.0+1.0/(2.0*(n+1.0))),2)*(1.0-z)))*(fabs(b)+fabs(c));

    } while (R != 0.0);

    result *= (128.0/3.141592654);

  }
  else {
  
    do {

      result += pow(a,2)*((n+1.0)/(n+2.0))*pow(z,n);
      R = a*(8.0*(n+2.0)/(n+3.0))*(pow(z,(n+1.0))/(1.0-z));
      n++;
      a *= ((n+0.5)/(n+1.0));

    } while (R != 0.0);

    result *= 8.0;

  }

  return result;

}
