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
    b[i].real = (float) i;
    b[i].imag = (float) i;
  }

  test( b );

  //test modification done in function test
  printf("In main function the values in b: \n");
  for(i=0; i<2; i++) {
    printf("b[%d].real = %f, b[%d].imag = %f\n", i, (b+i)->real, i, (b+i)->imag);
  }
  //printf("\n Real part of a: %f; Imaginary part of a: %f \n", a.real, a.imag);
  //printf("\n Mag Sq of b[1]: %f \n", mag_sq);

  return 1;

}//main

void test(complex *vec) {
  int i;
  complex *temp = vec;
  printf("In function test the values in b are: \n");
  for(i=0; i<2; i++) {
    CMULREAL( *(vec+i), 2.0, *(vec+i)); //multiply each element by 2
    CONJG( vec[i], vec[i]); //take complex conjugate
    //(vec+i)->real = 0.0; (vec+i)->imag = 0.0; //set elements to zero
    temp = vec + i;
    printf("b[%d].real = %f, b[%d].imag = %f\n Mag = %f\n", 
	   i, (vec+i)->real, 
	   i, (vec+i)->imag,
	   (float)cabs_sq(temp) );
  }

}//test
