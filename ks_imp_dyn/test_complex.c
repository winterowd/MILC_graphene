// testing the complex functions used in d_congrad_fn_u1.c and other
// files that have been modified

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"

//test function
void test(complex *vec);

int main(int argc, char **argv) {

  int i;
  double mag_sq;
  fcomplex *a;
  fcomplex b[2];

  printf("You have entered %d arguments:\n", argc);
  for(i=0; i<argc; i++, *argv++) {
    printf("%s\n", *argv);
  }
  
  for(i=0; i<2; i++) {
    b[i].real = (float) 1;
    b[i].imag = (float) 1;
  }

  test( b );

  //test modification done in function test
  printf("In main function the values in b: \n");
  for(i=0; i<2; i++) {
    printf("b[%d].real = %f, b[%d].imag = %f\n", i, b[i].real, i, b[i].imag);
  }
  //printf("\n Real part of a: %f; Imaginary part of a: %f \n", a.real, a.imag);
  //printf("\n Mag Sq of b[1]: %f \n", mag_sq);

  //CMULREAL(*b, 3., *b);
  printf("Value of b[0]: %f, %f\n", b[0].real, b[0].imag);
    

  return 1;

}//main

void test(complex *vec) {
  int i;
  double norm_sq = 0.0, norm_sq2;
  printf("In function test the values in b are: \n");
  for(i=0; i<2; i++) {
    CMULREAL( vec[i], 2.0, vec[i]); //multiply each element by 2
    CONJG( vec[i], vec[i]); //take complex conjugate
    norm_sq2 += vec[i].real*vec[i].real + vec[i].imag*vec[i].imag;
    norm_sq += (double) cabs_sq( (vec+i) );
    printf("b[%d].real = %f, b[%d].imag = %f\n", i, vec[i].real, i, vec[i].imag);
  }
  printf("vector's norm squared: %f or %f\n", norm_sq, norm_sq2);
  CMULREAL(*vec, 3., *vec);

}//test
