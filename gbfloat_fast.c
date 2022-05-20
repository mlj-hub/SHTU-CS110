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
    # pragma omp parallel for
    for (unsigned int y = 0; y < a.dimY; ++y)
    {
        for (unsigned int x = 0; x < a.dimX; ++x)
        {
            float * pc = get_pixel(b, x, y);
            unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
            deta = fmin(deta, gv.min_deta);
            float Sum[3] = {0,0,0};
            int offset;
            float s_data = gv.sum[ext - deta];

            for (int i = deta; i < (gv.length-2*deta)/4*4+deta; i+=4)
            {
                offset = i - ext;

                float * add1 = get_pixel(a, x + offset, y);
                float * add2 = get_pixel(a, x + offset+1, y);
                float * add3 = get_pixel(a, x + offset+2, y);
                float * add4 = get_pixel(a, x + offset+3, y);

                float opt1 = gv.data[i];
                float opt2 = gv.data[i+1];
                float opt3 = gv.data[i+2];
                float opt4 = gv.data[i+3];

                Sum[0] += opt1 * add1[0];
                Sum[1] += opt1 * add1[1];
                Sum[2] += opt1 * add1[2];

                Sum[0] += opt2 * add2[0];
                Sum[1] += opt2 * add2[1];
                Sum[2] += opt2 * add2[2];
                
                Sum[0] += opt3 * add3[0];
                Sum[1] += opt3 * add3[1];
                Sum[2] += opt3 * add3[2];
                
                Sum[0] += opt4 * add4[0];
                Sum[1] += opt4 * add4[1];
                Sum[2] += opt4 * add4[2];
            }

            for (int i = (gv.length-2*deta)/4*4+deta;i<gv.length-deta; ++i){
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

// SIMD and parallel
    // int max_threads = omp_get_max_threads();
    // omp_set_num_threads(max_threads);
    // # pragma omp parallel for
    // for (y = 0; y < a.dimY; ++y)
    // {
    //     for (x = 0; x < a.dimX; ++x)
    //     {
    //         pc = get_pixel(b, x, y);
    //         unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
    //         deta = fmin(deta, gv.min_deta);
    //         float Sum[3] = {0,0,0};
    //         float test_sum = 0;
    //         __m256 c_temp = _mm256_set1_ps(gv.sum[ext - deta]);
    //         float temp_result[8];
    //         for (i = deta; i < (gv.length-deta)/8*8; i+=8)
    //         {
    //             offset = i - ext;
    //             __m256 gvd = _mm256_loadu_ps(gv.data+i);

    //             __m256 result = _mm256_div_ps(gvd,c_temp);
    //             __m256 opt0 = _mm256_setr_ps((float)get_pixel(a, x + offset, y)[0],(float)get_pixel(a, x + offset+1, y)[0],\
    //             (float)get_pixel(a, x + offset+2, y)[0],(float)get_pixel(a, x + offset+3, y)[0],\
    //             (float)get_pixel(a, x + offset+4, y)[0],(float)get_pixel(a, x + offset+5, y)[0],\
    //             (float)get_pixel(a, x + offset+6, y)[0],(float)get_pixel(a, x + offset+7, y)[0]);
    //             result = _mm256_mul_ps(result,opt0);
    //             _mm256_storeu_ps(temp_result,result);
    //             for(int j = 0;j<8;++j){
    //                 // Sum[1]+=temp_result[j];
    //                 Sum[0]+=gv.data[i+j]/gv.sum[ext - deta] * (float)get_pixel(a, x + offset+j, y)[0];
    //             }

    //             result = _mm256_div_ps(gvd,c_temp);
    //             __m256 opt1 = _mm256_setr_ps((float)get_pixel(a, x + offset, y)[1],(float)get_pixel(a, x + offset+1, y)[1],\
    //             (float)get_pixel(a, x + offset+2, y)[1],(float)get_pixel(a, x + offset+3, y)[1],\
    //             (float)get_pixel(a, x + offset+4, y)[1],(float)get_pixel(a, x + offset+5, y)[1],\
    //             (float)get_pixel(a, x + offset+6, y)[1],(float)get_pixel(a, x + offset+7, y)[1]);
    //             result = _mm256_mul_ps(result,opt1);
    //             _mm256_storeu_ps(temp_result,result); 
    //             for(int j = 0;j<8;++j){
    //                 // Sum[1]+=temp_result[j];
    //                 Sum[1]+=gv.data[i+j]/gv.sum[ext - deta] * (float)get_pixel(a, x + offset+j, y)[1];
    //             }

    //             // static int cnt = 0;
    //             // for(int j = 0;j<8&&cnt<20;++j){
    //             //     printf("%f ",temp_result[j]);
    //             // }
    //             // if(cnt<20)
    //             //     printf("---- SIMD sum: %f \n",Sum[1]);

    //             // for(int j = 0;j<8&&cnt<20;++j){
    //             //     test_sum+=gv.data[i+j]/gv.sum[ext - deta] * (float)get_pixel(a, x + offset+j, y)[1];
    //             //     printf("%f ",gv.data[i+j]/gv.sum[ext - deta] * (float)get_pixel(a, x + offset+j, y)[1]);
    //             // }
    //             // if(cnt<20)
    //             //     printf("     Norm sum: %f \n",test_sum);

    //             // cnt++;

    //             result = _mm256_div_ps(gvd,c_temp);
    //             __m256 opt2 = _mm256_setr_ps((float)get_pixel(a, x + offset, y)[2],(float)get_pixel(a, x + offset+1, y)[2],\
    //             (float)get_pixel(a, x + offset+2, y)[2],(float)get_pixel(a, x + offset+3, y)[2],\
    //             (float)get_pixel(a, x + offset+4, y)[2],(float)get_pixel(a, x + offset+5, y)[2],\
    //             (float)get_pixel(a, x + offset+6, y)[2],(float)get_pixel(a, x + offset+7, y)[2]);
    //             result = _mm256_mul_ps(result,opt2);
    //             _mm256_storeu_ps(temp_result,result);
    //             for(int j = 0;j<8;++j){
    //                 // Sum[1]+=temp_result[j];
    //                 Sum[2]+=gv.data[i+j]/gv.sum[ext - deta] * (float)get_pixel(a, x + offset+j, y)[2];
    //             }
    //         }


    //         for (i = (gv.length-deta)/8*8; i < gv.length-deta; ++i)
    //         {
    //             offset = i - ext;
    //             Sum[0] += gv.data[i]/gv.sum[ext - deta] * (float)get_pixel(a, x + offset, y)[0];
    //             Sum[1] += gv.data[i]/gv.sum[ext - deta] * (float)get_pixel(a, x + offset, y)[1];
    //             Sum[2] += gv.data[i]/gv.sum[ext - deta] * (float)get_pixel(a, x + offset, y)[2];
    //         }

    //         pc[0] = Sum[0];
    //         pc[1] = Sum[1];
    //         pc[2] = Sum[2];
    //     }
    // }

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

// parallel
    // int max_threads = omp_get_max_threads();
    // omp_set_num_threads(max_threads);
    // # pragma omp parallel for
    // for (x = 0; x < a.dimX; ++x)
    // {
    //     for (y = 0; y < a.dimY; ++y)
    //     {
    //         pc = get_pixel(b, x, y);
    //         unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
    //         deta = fmin(deta, gv.min_deta);
    //         float Sum[3] = {0,0,0};
    //         for (i = deta; i < gv.length-deta; ++i)
    //         {
    //             offset = i - ext;
    //             Sum[0] += gv.data[i] /gv.sum[ext - deta] * (float)get_pixel(a, x, y + offset)[0];
    //             Sum[1] += gv.data[i] /gv.sum[ext - deta] * (float)get_pixel(a, x, y + offset)[1];
    //             Sum[2] += gv.data[i] /gv.sum[ext - deta] * (float)get_pixel(a, x, y + offset)[2];
    //         }
    //         pc[0] = Sum[0];
    //         pc[1] = Sum[1];
    //         pc[2] = Sum[2];
    //     }
    // }

// parallel
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    # pragma omp parallel for
    for (unsigned int x = 0; x < a.dimX; ++x)
    {
        for (unsigned int y = 0; y < a.dimY; ++y)
        {
            float *pc = get_pixel(b, x, y);
            unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
            deta = fmin(deta, gv.min_deta);
            float Sum[3] = {0,0,0};
            int offset;
            float s_data = gv.sum[ext - deta];

            for (int i = deta; i < (gv.length-2*deta)/4*4+deta; i+=4)
            {
                offset = i - ext;
            
                float * add1 = get_pixel(a, x , y + offset);
                float * add2 = get_pixel(a, x , y + offset+1);
                float * add3 = get_pixel(a, x , y + offset+2);
                float * add4 = get_pixel(a, x , y + offset+3);

                float opt1 = gv.data[i];
                float opt2 = gv.data[i+1];
                float opt3 = gv.data[i+2];
                float opt4 = gv.data[i+3];

                Sum[0] += opt1 * add1[0];
                Sum[1] += opt1 * add1[1];
                Sum[2] += opt1 * add1[2];
                
                Sum[0] += opt2 * add2[0];
                Sum[1] += opt2 * add2[1];
                Sum[2] += opt2 * add2[2];
                
                Sum[0] += opt3 * add3[0];
                Sum[1] += opt3 * add3[1];
                Sum[2] += opt3 * add3[2];
                
                Sum[0] += opt4 * add4[0];
                Sum[1] += opt4 * add4[1];
                Sum[2] += opt4 * add4[2];
            }

            for (int i = (gv.length-2*deta)/4*4+deta;i<gv.length-deta; ++i)
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
    // struct timeval start_time, stop_time, elapsed_time; 
    // gettimeofday(&start_time,NULL);
    Image b = gb_h(a, gv);
    // gettimeofday(&stop_time,NULL);
    // timersub(&stop_time, &start_time, &elapsed_time); 
    // printf("gb_h time: %f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    // gettimeofday(&start_time,NULL);
    Image c = gb_v(b, gv);
    // gettimeofday(&stop_time,NULL);
    // timersub(&stop_time, &start_time, &elapsed_time); 
    // printf("gb_v time: %f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

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