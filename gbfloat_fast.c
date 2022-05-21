#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <algorithm>
#include <sys/time.h>
#include <time.h>
#include <immintrin.h>
#include <omp.h>
// #include <emmintrin.h>
#include <immintrin.h>
//inplement dymanic

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define PI 3.14159

#define BLOCK_SIZE 16

typedef struct FVec
{
    unsigned int length;
    unsigned int min_length;
    unsigned int min_deta;
    float* data;
    float* sum;
} FVec;

typedef struct Image
{
    unsigned int dimX, dimY, numChannels;
    float* data;
} Image;

void normalize_FVec(FVec v)
{
    // float sum = 0.0;
    unsigned int i,j;
    int ext = v.length / 2;
    v.sum[0] = v.data[ext];
    for (i = ext+1,j=1; i < v.length; i++,j++)
    {
        v.sum[j] = v.sum[j-1] + v.data[i]*2;
    }
    // for (i = 0; i <= ext; i++)
    // {
    //      v.data[i] /= v.sum[v.length - ext - 1 ] ;
    //      printf("%lf ",v.sum[i]);
    // }
}

float gd(float a, float b, float x,double temp)
{
    float c = (x-b) / a;
    return exp((-.5) * c * c) / (a * temp);
}

float o_gd(float a, float x){
    float c = x / a;
    return exp((-.5) * c * c) / a;
}

/*can be optimised*/
FVec make_gv(float a, float x0, float x1, unsigned int length, unsigned int min_length)
{
    FVec v;
    v.length = length;
    v.min_length = min_length;
    if(v.min_length > v.length){
        v.min_deta = 0;
    }else{
        v.min_deta = ((v.length - v.min_length) / 2);
    }
    v.data = malloc(length * sizeof(float));
    v.sum = malloc((length / 2 + 1)* sizeof(float));
    float step = (x1 - x0) / ((float)length);
    int offset = length/2;

    // multi threads 

    float temp =(float) sqrt(2 * PI);
    __m256 c_temp = _mm256_set1_ps(temp);
    // int max_threads = omp_get_max_threads();
    // omp_set_num_threads(max_threads);
    // #pragma omp parallel for
    for (int i = 0; i < length/8*8; i+=8){
        __m256 opt1 = _mm256_set_ps(o_gd(a,(i+7-offset)*step),o_gd(a,(i+6-offset)*step),\
        o_gd(a,(i+5-offset)*step),o_gd(a,(i+4-offset)*step),\
        o_gd(a,(i+3-offset)*step),o_gd(a,(i+2-offset)*step),\
        o_gd(a,(i+1-offset)*step),o_gd(a,(i-offset)*step));
        __m256 data1  = _mm256_div_ps(opt1,c_temp);


        _mm256_storeu_ps(v.data+i,data1);

        // v.data[i] = gd(a, 0.0f, (i-offset)*step,temp);
    }
    
    for(int i=length/8*8;i<length;++i){
        v.data[i] = gd(a, 0.0f, (i-offset)*step,temp);
    }
    normalize_FVec(v);
    return v;
}

void print_fvec(FVec v)
{
    unsigned int i;
    printf("\n");
    for (i = 0; i < v.length; i++)
    {
        printf("%f ", v.data[i]);
    }
    printf("\n");
}

Image img_sc(Image a)
{
    Image b = a;
    b.data = malloc(b.dimX * b.dimY * b.numChannels * sizeof(float));
    return b;
}

float* get_pixel(Image img, int x, int y)
{
    if (x < 0)
    {
        x = 0;
    }
    if (x >= img.dimX)
    {
        x = img.dimX - 1;
    }
    if (y < 0)
    {
        y = 0;
    }
    if (y >= img.dimY)
    {
        y = img.dimY - 1;
    }
    return img.data + img.numChannels * (y * img.dimX + x);
}

Image gb_h(Image a, FVec gv)
{
    Image b = img_sc(a);

    int ext = gv.length / 2;

    // int max_threads = omp_get_max_threads();
    // omp_set_num_threads(max_threads);
    // # pragma omp parallel for
    // for (y = 0; y < a.dimY; y+=BLOCK_SIZE)
    // {
    //     for (x = 0; x < a.dimX; x+=BLOCK_SIZE)
    //     {
    //         for(int Y=y;Y<y+BLOCK_SIZE&&Y<a.dimY;++Y)
    //         {
    //             for(int X=x;X<x+BLOCK_SIZE&&X<a.dimX;++X)
    //             {
    //                 pc = get_pixel(b, X, Y);
    //                 unsigned int deta = fmin(fmin(a.dimY-Y-1, Y),fmin(a.dimX-X-1, X));
    //                 deta = fmin(deta, gv.min_deta);
    //                 float Sum[3] = {0,0,0};
    //                 for (i = deta; i < gv.length-deta; ++i)
    //                 {
    //                     offset = i - ext;
    //                     Sum[0] += gv.data[i]/gv.sum[ext - deta] * (float)get_pixel(a, X + offset, Y)[0];
    //                     Sum[1] += gv.data[i]/gv.sum[ext - deta] * (float)get_pixel(a, X + offset, Y)[1];
    //                     Sum[2] += gv.data[i]/gv.sum[ext - deta] * (float)get_pixel(a, X + offset, Y)[2];
    //                 }
    //                 pc[0] = Sum[0];
    //                 pc[1] = Sum[1];
    //                 pc[2] = Sum[2];
    //             }
    //         }
    //     }
    // }

// parallel
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    //omp_set_num_threads(2);
    # pragma omp parallel for
    for (int y = 0; y < a.dimY; ++y)
    {
        for (int x = 0; x < a.dimX; ++x)
        {
            float * pc = get_pixel(b, x, y);
            unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
            deta = fmin(deta, gv.min_deta);
            float Sum[3] = {0,0,0};
            int offset;
            float s_data = gv.sum[ext - deta];

            float Sum0[3][8] = {0.0f}; 
            
            __m256 Sum0256 = _mm256_setzero_ps();
            __m256 Sum1256 = _mm256_setzero_ps();
            __m256 Sum2256 = _mm256_setzero_ps();

            for (int i = deta; i < (gv.length-2*deta)/8*8+deta; i+=8)
            {
                offset = i - ext;

                float *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;
                /*
                float *add8;
                float *add9;
                float *add10;
                float *add11;
                float *add12;
                float *add13;
                float *add14;
                float *add15;
                
                float *add16;
                float *add17;
                float *add18;
                float *add19;
                float *add20;
                float *add21;
                float *add22;
                float *add23;

                float *add24;
                float *add25;
                float *add26;
                float *add27;
                float *add28;
                float *add29;
                float *add30;
                float *add31;
                
                */
                if(x+offset+7<=0)
                    add0=add1=add2=add3=add4=add5=add6=add7=a.data+3*(y*a.dimX);
                else if(x+offset>=(int)a.dimX-1)
                    add0=add1=add2=add3=add4=add5=add6=add7=a.data+3*(y*a.dimX+(int)a.dimX-1);
                if (x + offset >= 0 && x + offset + 7 <= a.dimX - 1)
                {
                    add0 =  a.data+3*(y*a.dimX+x+offset);
                    add1 = add0 + 3;
                    add2 = add0 + 6;
                    add3 = add0 + 9;
                    add4 = add0 + 12;
                    add5 = add0 + 15;
                    add6 = add0 + 18;
                    add7 = add0 + 21;
                }
                else
                {
                    add0 = get_pixel(a, x + offset, y);
                    add1 = get_pixel(a, x + offset + 1, y);
                    add2 = get_pixel(a, x + offset + 2, y);
                    add3 = get_pixel(a, x + offset + 3, y);
                    add4 = get_pixel(a, x + offset + 4, y);
                    add5 = get_pixel(a, x + offset + 5, y);
                    add6 = get_pixel(a, x + offset + 6, y);
                    add7 = get_pixel(a, x + offset + 7, y);
                }
                /*
                if(x+offset+15<=0)
                    add8=add9=add10=add11=add12=add13=add14=add15=a.data+3*(y*a.dimX+8);
                else if(x+offset+8>=(int)a.dimX-1)
                    add8=add9=add10=add11=add12=add13=add14=add15=a.data+3*(y*a.dimX+(int)a.dimX-1);
                if (x + offset+8 >= 0 && x + offset + 15 <= a.dimX - 1)
                {
                    add8 = a.data+3*(y*a.dimX+x+offset+8);
                    add9 = add8 + 3;
                    add10 = add8  + 6;
                    add11 = add8 + 9;
                    add12 = add8  + 12;
                    add13 = add8  + 15;
                    add14 = add8  + 18;
                    add15 = add8  + 21;
                }
                else
                {
                    add8 = get_pixel(a, x + offset+8, y);
                    add9 = get_pixel(a, x + offset + 9, y);
                    add10 = get_pixel(a, x + offset + 10, y);
                    add11 = get_pixel(a, x + offset + 11, y);
                    add12 = get_pixel(a, x + offset + 12, y);
                    add13 = get_pixel(a, x + offset + 13, y);
                    add14 = get_pixel(a, x + offset + 14, y);
                    add15 = get_pixel(a, x + offset + 15, y);
                }
                
                if(x+offset+23<=0)
                    add16=add17=add18=add19=add20=add21=add22=add23=a.data+3*(y*a.dimX+16);
                else if(x+offset+16>=(int)a.dimX-1)
                    add16=add17=add18=add19=add20=add21=add22=add23=a.data+3*(y*a.dimX+(int)a.dimX-1);
                if (x + offset+16 >= 0 && x + offset + 23 <= a.dimX - 1)
                {
                    add16 = a.data+3*(y*a.dimX+x+offset+16);
                    add17 = add16+ 3;
                    add18 = add16  + 6;
                    add19 = add16+ 9;
                    add20 = add16  + 12;
                    add21 = add16  + 15;
                    add22 = add16  + 18;
                    add23 = add16  + 21;
                }
                else
                {
                    add16 = get_pixel(a, x + offset+16, y);
                    add17 = get_pixel(a, x + offset + 17, y);
                    add18 = get_pixel(a, x + offset + 18, y);
                    add19 = get_pixel(a, x + offset + 19, y);
                    add20 = get_pixel(a, x + offset + 20, y);
                    add21 = get_pixel(a, x + offset + 21, y);
                    add22 = get_pixel(a, x + offset + 22, y);
                    add23 = get_pixel(a, x + offset + 23, y);
                }
               
                if(x+offset+31<=0)
                    add24=add25=add26=add27=add28=add29=add30=add31=a.data+3*(y*a.dimX+24);
                else if(x+offset+24>=(int)a.dimX-1)
                    add24=add25=add26=add27=add28=add29=add30=add31=a.data+3*(y*a.dimX+(int)a.dimX-1);
                if (x + offset+24 >= 0 && x + offset + 31 <= a.dimX - 1)
                {
                    add24 = a.data+3*(y*a.dimX+x+offset+24);
                    add25 = add24+ 3;
                    add26 = add24 + 6;
                    add27 = add24 + 9;
                    add28 = add24  + 12;
                    add29 = add24  + 15;
                    add30 = add24  + 18;
                    add31 = add24  + 21;
                }
                else
                {
                    add24 = get_pixel(a, x + offset+24, y);
                    add25 = get_pixel(a, x + offset + 25, y);
                    add26 = get_pixel(a, x + offset + 26, y);
                    add27 = get_pixel(a, x + offset + 27, y);
                    add28 = get_pixel(a, x + offset + 28, y);
                    add29 = get_pixel(a, x + offset + 29, y);
                    add30 = get_pixel(a, x + offset + 30, y);
                    add31 = get_pixel(a, x + offset + 31, y);
                }
                */
                __m256 Data = _mm256_loadu_ps(gv.data + i);
                // __m256 Data2 = _mm256_loadu_ps(gv.data+i+8);
                // __m256 Data3 = _mm256_loadu_ps(gv.data+i+16);
                // __m256 Data4 = _mm256_loadu_ps(gv.data+i+24);
                

                __m256 Chan0256 = _mm256_setr_ps(add0[0], add1[0], add2[0], add3[0], add4[0],
                                                add5[0], add6[0], add7[0]);
                __m256 Chan1256 = _mm256_setr_ps(add0[1], add1[1], add2[1], add3[1], add4[1],
                                                add5[1], add6[1], add7[1]);
                __m256 Chan2256 = _mm256_setr_ps(add0[2], add1[2], add2[2], add3[2], add4[2],
                                                add5[2], add6[2], add7[2]);

                // __m256 R2Chan0256 = _mm256_setr_ps(add8[0], add9[0], add10[0], add11[0], add12[0],
                //                                 add13[0], add14[0], add15[0]);
                // __m256 R2Chan1256 = _mm256_setr_ps(add8[1], add9[1], add10[1], add11[1], add12[1],
                //                                 add13[1], add14[1], add15[1]);
                // __m256 R2Chan2256 = _mm256_setr_ps(add8[2], add9[2], add10[2], add11[2], add12[2],
                //                                 add13[2], add14[2], add15[2]);


                // __m256 R3Chan0256 = _mm256_setr_ps(add16[0], add17[0], add18[0], add19[0], add20[0],
                //                                 add21[0], add22[0], add23[0]);
                // __m256 R3Chan1256 = _mm256_setr_ps(add16[1], add17[1], add18[1], add19[1], add20[1],
                //                                 add21[1], add22[1], add23[1]);
                // __m256 R3Chan2256 = _mm256_setr_ps(add16[2], add17[2], add18[2], addps(_mm256_mul_ps(R2Chan0256, Data2), Sum0256);

                // Sum1256 = _mm256_add_ps(_mm256_mul_ps(R2Chan1256, Data2), Sum1256);
                // Sum2256 = _mm256_add_ps(_mm256_mul_ps(R2Chan2256, Data2), Sum2256);

                // Sum0256 = _mm256_add_ps(_mm256_mul_ps(R3Chan0256, Data3), Sum0256);
                // Sum1256 = _mm256_add_ps(_mm256_mul_ps(R3Chan1256, Data3), Sum1256);
                // Sum2256 = _mm256_add_ps(_mm256_mul_ps(R3Chan2256, Data3), Sum2256);

                // Sum0256 = _mm256_add_ps(_mm256_mul_ps(R4Chan0256, Data4), Sum0256);
                // Sum1256 = _mm256_add_ps(_mm256_mul_ps(R4Chan1256, Data4), Sum1256);
                // Sum2256 = _mm256_add_ps(

                // Sum0256 = _mm256_add_ps(_mm256_mul_ps(R3Chan0256, Data3), Sum0256);
                // Sum1256 = _mm256_add_ps(_mm256_mul_ps(R3Chan1256, Data3), Sum1256);
                // Sum2256 = _mm256_add_ps(_mm256_mul_ps(R3Chan2256, Data3), Sum2256);

                // Sum0256 = _mm256_add_ps(_mm256_mul_ps(R4Chan0256, Data4), Sum0256);
                // Sum1256 = _mm256_add_ps(_mm256_mul_ps(R4Chan1256, Data4), Sum1256);
                // Sum2256 = _mm256_add_ps(

                // Sum0256 = _mm256_add_ps(_mm256_mul_ps(R3Chan0256, Data3), Sum0256);
                // Sum1256 = _mm256_add_ps(_mm256_mul_ps(R3Chan1256, Data3), Sum1256);
                // Sum2256 = _mm256_add_ps(_mm256_mul_ps(R3Chan2256, Data3), Sum2256);

                Sum0256 = _mm256_add_ps(_mm256_mul_ps(Chan0256, Data), Sum0256);
                Sum1256 = _mm256_add_ps(_mm256_mul_ps(Chan1256, Data), Sum1256);
                Sum2256 = _mm256_add_ps(_mm256_mul_ps(Chan2256, Data), Sum2256);
            }
                _mm256_storeu_ps(Sum0[0],Sum0256);
                _mm256_storeu_ps(Sum0[1],Sum1256);
                _mm256_storeu_ps(Sum0[2],Sum2256);

                Sum[0] += Sum0[0][0]+Sum0[0][1]+Sum0[0][2]+Sum0[0][3]+Sum0[0][4]+Sum0[0][5]+Sum0[0][6]+Sum0[0][7];
                Sum[1] += Sum0[1][0]+Sum0[1][1]+Sum0[1][2]+Sum0[1][3]+Sum0[1][4]+Sum0[1][5]+Sum0[1][6]+Sum0[1][7];
                Sum[2] += Sum0[2][0]+Sum0[2][1]+Sum0[2][2]+Sum0[2][3]+Sum0[2][4]+Sum0[2][5]+Sum0[2][6]+Sum0[2][7];
                
            for (int i = (gv.length-2*deta)/8*8+deta;i<gv.length-deta; ++i){
                offset = i - ext;
                float data = gv.data[i];
                float * add = get_pixel(a, x + offset, y);
                Sum[0] += data * add[0];
                Sum[1] += data * add[1];
                Sum[2] += data * add[2];
            }

            pc[0] = Sum[0]/s_data;
            pc[1] = Sum[1]/s_data;
            pc[2] = Sum[2]/s_data;
        }
    }

    return b;
}

Image gb_v(Image a, FVec gv)
{
    Image b = img_sc(a);

    int ext = gv.length / 2;

// cache blocking 
//     int max_threads = omp_get_max_threads();
//     omp_set_num_threads(max_threads);
//     # pragma omp parallel for
//     for (unsigned int x = 0; x < a.dimX; x+=BLOCK_SIZE)
//     {
//         for (unsigned int y = 0; y < a.dimY; y+=BLOCK_SIZE)
//         {
//             for(int Y=y;Y<y+BLOCK_SIZE&&Y<a.dimY;++Y)
//             {
//                 for(int X=x;X<x+BLOCK_SIZE&&X<a.dimX;++X)
//                 {
//                     float* pc = get_pixel(b, X, Y);
//                     unsigned int deta = fmin(fmin(a.dimY-Y-1, Y),fmin(a.dimX-X-1, X));
//                     deta = fmin(deta, gv.min_deta);
//                     float Sum[3] = {0,0,0};
//                     int offset;
//                     for (int i = deta; i < gv.length-deta; ++i)
//                     {
//                         offset = i - ext;
//                         Sum[0] += gv.data[i] /gv.sum[ext - deta] * (float)get_pixel(a, X, Y + offset)[0];
//                         Sum[1] += gv.data[i] /gv.sum[ext - deta] * (float)get_pixel(a, X, Y + offset)[1];
//                         Sum[2] += gv.data[i] /gv.sum[ext - deta] * (float)get_pixel(a, X, Y + offset)[2];
//                     }
//                     pc[0] = Sum[0];
//                     pc[1] = Sum[1];
//                     pc[2] = Sum[2];
//                 }
//             }
//         }
//     }

// another cache blocking (just for testing )
    // int max_threads = omp_get_max_threads();
    // omp_set_num_threads(max_threads);
    // # pragma omp parallel for
    // for (x = 0; x < a.dimX; x+=BLOCK_SIZE)
    // {
    //     for (y = 0; y < a.dimY; y+=BLOCK_SIZE)
    //     {
    //         for(int Y=y;Y<y+BLOCK_SIZE&&Y<a.dimY;++Y)
    //         {
    //             for(int X=x;X<x+BLOCK_SIZE&&X<a.dimX;++X)
    //             {
    //                 pc = get_pixel(b, X, Y);
    //                 unsigned int deta = fmin(fmin(a.dimY-Y-1, Y),fmin(a.dimX-X-1, X));
    //                 deta = fmin(deta, gv.min_deta);
    //                 float Sum[3] = {0,0,0};
    //                 for (i = deta; i < gv.length-deta; i+=BLOCK_SIZE)
    //                 {
    //                     for(int I = i;I<i+BLOCK_SIZE && I<gv.length-deta;++I){
    //                         offset = I - ext;
    //                         Sum[0] += gv.data[I] /gv.sum[ext - deta] * (float)get_pixel(a, X, Y + offset)[0];
    //                         Sum[1] += gv.data[I] /gv.sum[ext - deta] * (float)get_pixel(a, X, Y + offset)[1];
    //                         Sum[2] += gv.data[I] /gv.sum[ext - deta] * (float)get_pixel(a, X, Y + offset)[2];
    //                     }
    //                 }
    //                 pc[0] = Sum[0];
    //                 pc[1] = Sum[1];
    //                 pc[2] = Sum[2];
    //             }
    //         }
    //     }
    // }

// parallel and SIMD
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
#pragma omp parallel for
    for (int x = 0; x < a.dimX; ++x)
    {
        for (int y = 0; y < a.dimY; ++y)
        {
            float *pc = get_pixel(b, x, y);
            unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
            deta = fmin(deta, gv.min_deta);
            float Sum[3] = {0,0,0};
            int offset;
            float s_data = gv.sum[ext - deta];

            float Sum0[3][8] = {0.0f}; 
            
            __m256 Sum0256 = _mm256_setzero_ps();
            __m256 Sum1256 = _mm256_setzero_ps();
            __m256 Sum2256 = _mm256_setzero_ps();    

            for (int i = deta; i < (gv.length-2*deta)/8*8+deta; i+=8)
            {
                offset = i - ext;

                float *add0;
                float *add1;
                float *add2;
                float *add3;
                float *add4;
                float *add5;
                float *add6;
                float *add7;
                if(y+offset+7<=0)
                    add0=add1=add2=add3=add4=add5=add6=add7=a.data+ 3*x;
                else if(y+offset>=(int)a.dimY-1)
                    add0=add1=add2=add3=add4=add5=add6=add7=a.data+ 3*(((int)a.dimY-1)*a.dimX+x);
                else if (y + offset >= 0 && y + offset + 7 <= (int)a.dimY - 1)
                {
                    add0 = a.data+ 3*((y+offset)*a.dimX+x);
                    add1 = add0 + 3* a.dimX;
                    add2 = add0 + 6 * a.dimX;
                    add3 = add0 + 9 * a.dimX;
                    add4 = add0 + 12 * a.dimX;
                    add5 = add0 + 15* a.dimX;
                    add6 = add0 + 18* a.dimX;
                    add7 = add0 + 21* a.dimX;
                }
                else
                {
                    add0 = get_pixel(a, x, y + offset);
                    add1 = get_pixel(a, x, y + offset + 1);
                    add2 = get_pixel(a, x, y + offset + 2);
                    add3 = get_pixel(a, x, y + offset + 3);
                    add4 = get_pixel(a, x, y + offset + 4);
                    add5 = get_pixel(a, x, y + offset + 5);
                    add6 = get_pixel(a, x, y + offset + 6);
                    add7 = get_pixel(a, x, y + offset + 7);
                }

                __m256 Data = _mm256_loadu_ps(gv.data + i);

                __m256 Chan0256 = _mm256_setr_ps(add0[0], add1[0], add2[0], add3[0], add4[0], add5[0], add6[0], add7[0]);
                __m256 Chan1256 = _mm256_setr_ps(add0[1], add1[1], add2[1], add3[1], add4[1], add5[1], add6[1], add7[1]);
                __m256 Chan2256 = _mm256_setr_ps(add0[2], add1[2], add2[2], add3[2], add4[2], add5[2], add6[2], add7[2]);

                Sum0256 = _mm256_add_ps(_mm256_mul_ps(Chan0256, Data), Sum0256);
                Sum1256 = _mm256_add_ps(_mm256_mul_ps(Chan1256, Data), Sum1256);
                Sum2256 = _mm256_add_ps(_mm256_mul_ps(Chan2256, Data), Sum2256);
            }

            _mm256_storeu_ps(Sum0[0],Sum0256);
            _mm256_storeu_ps(Sum0[1],Sum1256);
            _mm256_storeu_ps(Sum0[2],Sum2256);
            Sum[0] += Sum0[0][0]+Sum0[0][1]+Sum0[0][2]+Sum0[0][3]+Sum0[0][4]+Sum0[0][5]+Sum0[0][6]+Sum0[0][7];
            Sum[1] += Sum0[1][0]+Sum0[1][1]+Sum0[1][2]+Sum0[1][3]+Sum0[1][4]+Sum0[1][5]+Sum0[1][6]+Sum0[1][7];
            Sum[2] += Sum0[2][0]+Sum0[2][1]+Sum0[2][2]+Sum0[2][3]+Sum0[2][4]+Sum0[2][5]+Sum0[2][6]+Sum0[2][7];

            for (int i = (gv.length-2*deta)/8*8+deta;i<gv.length-deta; ++i)
            {
                offset = i - ext;
                float data = gv.data[i];
                float * add = get_pixel(a, x , y+ offset);
                Sum[0] += data * add[0];
                Sum[1] += data * add[1];
                Sum[2] += data * add[2];
            }

            pc[0] = Sum[0]/s_data;
            pc[1] = Sum[1]/s_data;
            pc[2] = Sum[2]/s_data;
        }
    }

    return b;
}

Image apply_gb(Image a, FVec gv)
{
    struct timeval start_time, stop_time, elapsed_time; 
    gettimeofday(&start_time,NULL);
    Image b = gb_h(a, gv);
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); 
    printf("gb_h time: %f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    gettimeofday(&start_time,NULL);
    Image c = gb_v(b, gv);
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); 
    printf("gb_v time: %f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    free(b.data);
    return c;
}

int main(int argc, char** argv)
{
    struct timeval start_time, stop_time, elapsed_time; 
    gettimeofday(&start_time,NULL);
    if (argc < 6)
    {
        printf("Usage: ./gb.exe <inputjpg> <outputname> <float: a> <float: x0> <float: x1> <unsigned int: dim>\n");
        exit(0);
    }

    float a, x0, x1;
    unsigned int dim, min_dim;

    sscanf(argv[3], "%f", &a);
    sscanf(argv[4], "%f", &x0);
    sscanf(argv[5], "%f", &x1);
    sscanf(argv[6], "%u", &dim);
    sscanf(argv[7], "%u", &min_dim);

    // struct timeval start_time1, stop_time1, elapsed_time1; 
    // gettimeofday(&start_time1,NULL);
    FVec v = make_gv(a, x0, x1, dim, min_dim);
    // gettimeofday(&stop_time1,NULL);
    // timersub(&stop_time1, &start_time1, &elapsed_time1); 
    // printf("make_gv time :%f \n", elapsed_time1.tv_sec+elapsed_time1.tv_usec/1000000.0);

    // print_fvec(v);
    Image img;
    img.data = stbi_loadf(argv[1], &(img.dimX), &(img.dimY), &(img.numChannels), 0);

    Image imgOut = apply_gb(img, v);
    stbi_write_jpg(argv[2], imgOut.dimX, imgOut.dimY, imgOut.numChannels, imgOut.data, 90);
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); 
    printf("%f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
    free(imgOut.data);
    free(v.data);
    free(v.sum);
    return 0;
}