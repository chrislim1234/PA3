#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include "bmplib.h"

using namespace std;

//============================Add function prototypes here======================

void convolve(unsigned char[][SIZE][3], unsigned char[][SIZE][3] , int, double[][11]);
void sobel(unsigned char[][SIZE][3], unsigned char[][SIZE][3]);
void gaussian(int[][11], int, double);
void gaussian_filter(unsigned char[][SIZE][3], unsigned char[][SIZE][3], int, double);
void dummy(unsigned char[][SIZE][RGB], unsigned char[][SIZE][RGB]);
void unsharp(unsigned char[][SIZE][RGB], unsigned char[][SIZE][RGB], int, double, double);

//============================Do not change code in main()======================

#ifndef AUTOTEST

int main(int argc, char* argv[])
{
    //First check argc
    if(argc < 3)
    {
        //we need at least ./filter <input file> <filter name> to continue
        cout << "usage: ./filter <input file> <filter name> <filter parameters>";
        cout << " <output file name>" << endl;
        return -1;
    }
    //then check to see if we can open the input file
    unsigned char input[SIZE][SIZE][RGB];
    unsigned char output[SIZE][SIZE][RGB];
    char* outfile;
    int N;
    double sigma, alpha;
    
    // read file contents into input array
    int status = readRGBBMP(argv[1], input);
    if(status != 0)
    {
        cout << "unable to open " << argv[1] << " for input." << endl;
        return -1;
    }
    //Input file is good, now look at next argument
    if( strcmp("sobel", argv[2]) == 0)
    {
        sobel(output, input);
        outfile = argv[3];
    }
    else if( strcmp("blur", argv[2]) == 0)
    {
        if(argc < 6)
        {
            cout << "not enough arguments for blur." << endl;
            return -1;
        }
        N = atoi(argv[3]);
        sigma = atof(argv[4]);
        outfile = argv[5];
        gaussian_filter(output, input, N, sigma);
    }
    else if( strcmp("unsharp", argv[2]) == 0)
    {
        if(argc < 7)
        {
            cout << "not enough arguments for unsharp." << endl;
            return -1;
        }
        N = atoi(argv[3]);
        sigma = atof(argv[4]);
        alpha = atof(argv[5]);
        outfile = argv[6];
        unsharp(output, input, N, sigma, alpha);
        
    }
    else if( strcmp("dummy", argv[2]) == 0)
    {
        //do dummy
        dummy(output, input);
        outfile = argv[3];
    }
    else
    {
        cout << "unknown filter type." << endl;
        return -1;
    }
    
    if(writeRGBBMP(outfile, output) != 0)
    {
        cout << "error writing file " << argv[3] << endl;
    }
    
}

#endif

//=========================End Do not change code in main()=====================


// Creates an identity kernel (dummy kernel) that will simply
// copy input to output image and applies it via convolve()
//
// ** This function is complete and need not be changed.
// Use this as an example of how to create a kernel array, fill it in
// appropriately and then use it in a call to convolve.
void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
    double k[11][11];
    for (int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            k[i][j] = 0;
        }
    }
    k[1][1] = 1;
    convolve(out, in, 3, k);
}


// Convolves an input image with an NxN kernel to produce the output kernel
// You will need to complete this function by following the
//  instructions in the comments
void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB],
              int N, double kernel[][11])
{
    
    int padded[SIZE+10][SIZE+10][RGB];  // Use for input image with appropriate
    // padding
    int temp[SIZE][SIZE][RGB];          // Use for the unclamped output pixel
    // values then copy from temp to out,
    // applying clamping
    
    //first set all of padded to 0 (black)
    for (int i=0; i<SIZE+10; i++) {
        for (int j=0; j<SIZE+10; j++) {
            for(int c=0; c<RGB; c++) {
                padded[i][j][c] = 0;
            }
        }
    }
    
    //now copy input into padding to appropriate locations
    for (int i=0; i<SIZE; i++) {
        for (int j=0; j<SIZE; j++) {
            for(int c=0; c<RGB; c++) {
                padded[i][j][c] =in[i][j][c];
            }
        }
    }
    
    //initialize temp pixels to 0 (black)
    for (int i=0; i<SIZE; i++) {
        for (int j=0; j<SIZE; j++) {
            for(int c=0; c<RGB; c++) {
                temp[i][j][c] = 0;
            }
        }
    }
    
    
    //now perform convolve (using convolution equation on each pixel of the
    // actual image) placing the results in temp (i.e. unclamped result)
    //Here we give you the structure of the convolve for-loops, you need
    //to figure out the loop limits
    
    for(int k=0;k<RGB;k++) {
      for(int y=0;y<SIZE ;y++) {
        for(int x=0 ;x<SIZE ;x++) {
              for(int i=-N/2 ; i<=N/2 ; i++) {
                for(int j=-N/2 ; j<=N/2 ; j++) {
                  if ((y+i)>=0 && (x+j)>=0 && (y+i)<=SIZE && (x+j)<=255 ) {
                temp[y][x][k] += padded[y+i][x+j][k]*kernel[N/2+i][N/2+j];
              }
       }
      }
      }
      }
      }
                
                
                //now clamp and copy to output
                // You may need to cast to avoid warnings from the compiler:
                // (i.e. out[i][j][k] = (unsigned char) temp[i][j][k];)
                
                for (int i=0; i<SIZE; i++) {
                    for (int j=0; j<SIZE; j++) {
                        for(int c=0; c<RGB; c++) {
                            int pix = temp[i][j][c];
                            
                            if (pix<0) {
                                pix=0;
                            }
                            else if (pix>255) {
                                pix=255;
                            }
                            out[i][j][c] = (unsigned char) pix;
                        }
                    }
                }
            }
                


void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
    double k[11][11];
    double s_h1[3][3] = { {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1} };
    double s_h2[3][3] = { {1, 0, -1},
        {2, 0, -2},
        {1, 0, -1} };
    
    unsigned char h1_soble[SIZE][SIZE][RGB]; //hold intemediate images
    unsigned char h2_soble[SIZE][SIZE][RGB];
    
    for (int i = 0; i < 11; i++)
    {
        for(int j=0; j < 11; j++)
        {
            k[i][j] = 0;
        }
    }
    
    
    // Copy in 1st 3x3 horizontal sobel kernel (i.e. copy s_h1 into k)
    
    for (int i = 0; i < 3; i++)
    {
        for(int j=0; j < 3; j++)
        {
            
            k[i][j] = s_h1[i][j];
            
        }
    }
    
    // Call convolve to apply horizontal sobel kernel with result in h1_soble
    //unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB],
    //int N, double kernel[][11]
    
    convolve(h1_soble,in,3,k);
    
    // Copy in 2nd 3x3 horizontal sobel kernel (i.e. copy s_h2 into k)
    for (int i = 0; i < 3; i++)
    {
        for(int j=0; j < 3; j++)
        {
            
            k[i][j] = s_h2[i][j];
            
        }
    }
    
    
    // Call convolve to apply horizontal sobel kernel with result in h2_soble
    
    convolve(h2_soble,in,3,k);
    
    // Add the two results (applying clamping) to produce the final output
    for (int i=0; i<SIZE; i++) {
        for (int j=0; j<SIZE; j++) {
            for(int c=0; c<RGB; c++) {
                double sum = h1_soble[i][j][c] + h2_soble[i][j][c];
                    out[i][j][c] = sum;
                if (sum<0) {
                    out[i][j][c]=0;
                }
                else if (sum>255) {
                    out[i][j][c]=255;
                }
                out[i][j][c]= (unsigned char) out[i][j][c];
            }
        }
    }
}

// Add the rest of your functions here (i.e. gaussian, gaussian_filter, unsharp)
void gaussian(double k[][11],int N, double sigma)
{
    double total = 0;
    //create k, and make it all black
    for (int i = 0; i < 11; i++)
    {
        for(int j=0; j < 11; j++)
        {
            k[i][j] = 0;
        }
    }
    //calculate guassin formula and find total
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            k[i][j] = 1.0 * exp(-(((pow((i-N/2),2)/(2*pow(sigma,2)))+(pow((j-N/2),2)/(2*pow(sigma,2))))));
            total=total+k[i][j];
        }
    }
    //divide by total so that they add up to 1 at the end
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            k[i][j]=k[i][j]/total;
        }
    }
    
}

void gaussian_filter(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma)
{
    //create k
    double k[11][11];
    for (int i=0; i<11; i++) {
        for (int j=0; j<11; j++) {
            k[i][j]=0;
        }
    }
    //put into gaussian
    gaussian(k,N,sigma);
    
    //create new image
    convolve(out,in,N,k);
    
}

void unsharp(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int
             N, double sigma, double alpha)
{
    //Create a blur array
    unsigned char blur[SIZE][SIZE][RGB];
    
    //create k
    double k[11][11];
    for (int i=0; i<11; i++) {
        for (int j=0; j<11; j++) {
            k[i][j]=0;
        }
    }
    //create blur() function
    //𝐵 = 𝑏𝑙𝑢𝑟(𝐼𝑀)
    gaussian_filter(blur,in,N,sigma);
    //create new image
    
    //create D and S arrays
    double d[SIZE][SIZE][RGB];
    double s[SIZE][SIZE][RGB];
    
    //𝐷 = 𝐼𝑀 − 𝐵
    for (int i=0; i<SIZE; i++) {
        for (int j=0; j<SIZE; j++) {
            for (int c=0; c<RGB; c++) {
                d[i][j][c]= in[i][j][c]-blur[i][j][c];
            }
        }
    }
    
    //𝑆 = 𝐼𝑀 + 𝛼𝐷
    for (int ii=0; ii<SIZE; ii++) {
        for (int jj=0; jj<SIZE; jj++) {
            for (int cc=0; cc<RGB; cc++) {
                s[ii][jj][cc]= in[ii][jj][cc]+ alpha * d[ii][jj][cc];
            }
        }
    }
        
        for (int a=0; a<SIZE; a++) {
            for (int b=0; b<SIZE; b++) {
                for (int c=0; c<RGB; c++) {
                    double pix =  s[a][b][c];
                    
                    if (pix<0) {
                        pix=0;
                    }
                    else if (pix>255) {
                        pix=255;
                    }
                    out[a][b][c] = (unsigned char) pix;
                }
            }
        }
        
    
}
    
    


