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
#define BIT_SIZE 256
#define UNROLLING 8

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


Image img_copy(Image a)
{
    Image b = a;
    b.dimX = a.dimY;
    b.dimY = a.dimX;
    b.data = malloc(b.dimX * b.dimY * b.numChannels * sizeof(float));
    return b;
}

Image transpose(Image a){
    Image b = img_copy(a);
    int max_threads = omp_get_max_threads();

    omp_set_num_threads(max_threads);
    // #pragma omp parallel for
    #pragma omp parallel for schedule(dynamic)
    for(int x=0;x<(int)a.dimX;x+=BLOCK_SIZE){
        for(int y=0;y<(int)a.dimY;y+=BLOCK_SIZE){
            for(int X = x;X<x+BLOCK_SIZE&&X<(int)a.dimX;++X){
                for(int Y = y;Y<y+BLOCK_SIZE&&Y<(int)a.dimY;++Y){
                    int data1 = 3*(X*b.dimX+Y); 
                    int data2 = 3*(Y*a.dimX+X);
                    (b.data+data1)[0]=(a.data+data2)[0];
                    (b.data+data1)[1]=(a.data+data2)[1];
                    (b.data+data1)[2]=(a.data+data2)[2];
                }
            }
        }
    }
    return b;
}

void get_RGB(Image a,Image * r,Image * g,Image*b){
    *r = *b = *g =a;
    r->numChannels=g->numChannels=b->numChannels=1;
    r->data = malloc(a.dimX*a.dimY*sizeof(float));
    g->data = malloc(a.dimX*a.dimY*sizeof(float));
    b->data = malloc(a.dimX*a.dimY*sizeof(float));
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    #pragma omp parallel for
    // #pragma omp parallel for schedule(dynamic)
    for(int y=0;y<(int)a.dimY;y+=BLOCK_SIZE)
    {
        for(int x=0;x<(int)a.dimX;x+=BLOCK_SIZE)
        {
            for(int Y = y;Y<y+BLOCK_SIZE&&Y<(int)a.dimY;++Y)
            {
                for(int X = x;X<x+BLOCK_SIZE&&X<(int)a.dimX;++X)
                {
                    int temp = Y*a.dimX +X;
                    int mul = 3*temp;
                    (r->data+temp)[0]=(a.data+mul)[0];
                    (g->data+temp)[0]=(a.data+mul)[1];
                    (b->data+temp)[0]=(a.data+mul)[2];
                }
            }
        }
    }
}

void get_RGB_n(Image a,Image * r,Image * g,Image*b){
    *r = *b = *g =a;
    int size = (a.dimX+16)*a.dimY*sizeof(float);
    r->numChannels=g->numChannels=b->numChannels=1;
    r->dimX=g->dimX=b->dimX=a.dimX+16;
    r->data = malloc(size);
    g->data = malloc(size);
    b->data = malloc(size);
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    #pragma omp parallel for
    // #pragma omp parallel for schedule(dynamic)
    for(int y=0;y<(int)a.dimY;y+=BLOCK_SIZE)
    {
        for(int x=0;x<(int)a.dimX;x+=BLOCK_SIZE)
        {
            for(int Y = y;Y<y+BLOCK_SIZE&&Y<(int)a.dimY;++Y)
            {
                for(int X = x;X<x+BLOCK_SIZE&&X<(int)a.dimX;++X)
                {
                    int temp = Y*a.dimX +X;
                    int temp2 = Y*r->dimX+X+8;
                    int mul = 3*temp;
                    (r->data+temp2)[0]=(a.data+mul)[0];
                    (g->data+temp2)[0]=(a.data+mul)[1];
                    (b->data+temp2)[0]=(a.data+mul)[2];
                }
            }
        }
    }
    int offset1 = 8+a.dimX;
    int offset2 = 7+a.dimX;
    max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    #pragma omp parallel for
    // #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<r->dimY*(r->dimX);i+=r->dimX){
        memset(r->data+i,*(r->data+i+8),8*sizeof(float));
        memset(g->data+i,*(g->data+i+8),8*sizeof(float));
        memset(b->data+i,*(b->data+i+8),8*sizeof(float));
        memset(r->data+i+offset1,*(r->data+i+offset2),8*sizeof(float));
        memset(g->data+i+offset1,*(g->data+i+offset2),8*sizeof(float));
        memset(b->data+i+offset1,*(b->data+i+offset2),8*sizeof(float));
    }
}

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
    for (int i = 0; i < (int)length/8*8; i+=8){
        __m256 opt1 = _mm256_set_ps(o_gd(a,(i+7-offset)*step),o_gd(a,(i+6-offset)*step),\
        o_gd(a,(i+5-offset)*step),o_gd(a,(i+4-offset)*step),\
        o_gd(a,(i+3-offset)*step),o_gd(a,(i+2-offset)*step),\
        o_gd(a,(i+1-offset)*step),o_gd(a,(i-offset)*step));
        __m256 data1  = _mm256_div_ps(opt1,c_temp);
        _mm256_storeu_ps(v.data+i,data1);
    }
    
    for(int i=length/8*8;i<(int)length;++i){
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
    if (x >=(int) img.dimX)
    {
        x = img.dimX - 1;
    }
    if (y < 0)
    {
        y = 0;
    }
    if (y >= (int)img.dimY)
    {
        y = img.dimY - 1;
    }
    return img.data + img.numChannels * (y * img.dimX + x);
}

int get_pixel_h(int dimX,int x,int y){
    if (x < 0)
        x = 0;
    if (x >= dimX)
        x = dimX - 1;
    return y*dimX+x;
}

Image gb_h(Image a, Image r,Image g,Image b,FVec gv)
{
    Image o = img_sc(a);

    int ext = gv.length / 2;

// parallel
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    // # pragma omp parallel for schedule(dynamic)
    # pragma omp parallel for
    for (int y = 0; y < (int)a.dimY; ++y)
    {
        for (int x = 0; x < (int)a.dimX; ++x)
        {
            float * pc = get_pixel(o, x, y);
            unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
            deta = fmin(deta, gv.min_deta);
            float Sum[3] = {0,0,0};
            int offset;
            float s_data = gv.sum[ext - deta];

            float Sum0[3][8] = {{0.0f,.0f,.0f},{0.0f,.0f,.0f},{0.0f,.0f,.0f}}; 
            
            __m256 Sum0256 = _mm256_setzero_ps();
            __m256 Sum1256 = _mm256_setzero_ps();
            __m256 Sum2256 = _mm256_setzero_ps();
            //#pragma UNROLL(8)
            for (unsigned int i = (int)deta; i < (gv.length-2*deta)/UNROLLING*UNROLLING+deta; i+=UNROLLING)
            {
                offset = i - ext;
                int pixel0,pixel1,pixel2,pixel3,pixel4,pixel5,pixel6,pixel7;
                __m256 Data = _mm256_loadu_ps(gv.data + i);

                __m256 DataR;
                __m256 DataG;
                __m256 DataB;
                /*
                __m256 DataR2;
                __m256 DataG2;
                __m256 DataB2;
                */        
                if(x+offset+7<=0){
                    int temp = y*(int)r.dimX;
                    DataR=_mm256_set1_ps((r.data+temp)[0]);
                    DataG=_mm256_set1_ps((g.data+temp)[0]);
                    DataB=_mm256_set1_ps((b.data+temp)[0]);
                }
                else if (x+offset>=(int)a.dimX-1){
                    int temp = (y+1)*(int)r.dimX-1;
                    DataR=_mm256_set1_ps((r.data+temp)[0]);
                    DataG=_mm256_set1_ps((g.data+temp)[0]);
                    DataB=_mm256_set1_ps((b.data+temp)[0]);
                }
                else
                {
                    int temp = y*(int)r.dimX+x+8+offset;
                    DataR = _mm256_loadu_ps(r.data+temp);
                    DataG = _mm256_loadu_ps(g.data+temp);
                    DataB = _mm256_loadu_ps(b.data+temp);
                }
                // else
                // {
                //     pixel0=get_pixel_h(a.dimX,x+offset,y);
                //     pixel1=get_pixel_h(a.dimX,x+offset+1,y);
                //     pixel2=get_pixel_h(a.dimX,x+offset+2,y);
                //     pixel3=get_pixel_h(a.dimX,x+offset+3,y);
                //     pixel4=get_pixel_h(a.dimX,x+offset+4,y);
                //     pixel5=get_pixel_h(a.dimX,x+offset+5,y);
                //     pixel6=get_pixel_h(a.dimX,x+offset+6,y);
                //     pixel7=get_pixel_h(a.dimX,x+offset+7,y);
                //     DataR = _mm256_setr_ps((r.data+pixel0)[0],(r.data+pixel1)[0],(r.data+pixel2)[0],\
                //     (r.data+pixel3)[0],(r.data+pixel4)[0],(r.data+pixel5)[0],(r.data+pixel6)[0],\
                //     (r.data+pixel7)[0]);
                //     DataG = _mm256_setr_ps((g.data+pixel0)[0],(g.data+pixel1)[0],(g.data+pixel2)[0],\
                //     (g.data+pixel3)[0],(g.data+pixel4)[0],(g.data+pixel5)[0],(g.data+pixel6)[0],\
                //     (g.data+pixel7)[0]);
                //     DataB = _mm256_setr_ps((b.data+pixel0)[0],(b.data+pixel1)[0],(b.data+pixel2)[0],\
                //     (b.data+pixel3)[0],(b.data+pixel4)[0],(b.data+pixel5)[0],(b.data+pixel6)[0],\
                //     (b.data+pixel7)[0]);
                // }
                DataR = _mm256_mul_ps(DataR, Data);
                DataG = _mm256_mul_ps(DataG, Data);
                DataB = _mm256_mul_ps(DataB, Data);
                
                Sum0256 = _mm256_add_ps(DataR, Sum0256);
                Sum1256 = _mm256_add_ps(DataG, Sum1256);
                Sum2256 = _mm256_add_ps(DataB, Sum2256);
            }
                _mm256_storeu_ps(Sum0[0],Sum0256);
                _mm256_storeu_ps(Sum0[1],Sum1256);
                _mm256_storeu_ps(Sum0[2],Sum2256);

                Sum[0] += Sum0[0][0]+Sum0[0][1]+Sum0[0][2]+Sum0[0][3]+Sum0[0][4]+Sum0[0][5]+Sum0[0][6]+Sum0[0][7];
                Sum[1] += Sum0[1][0]+Sum0[1][1]+Sum0[1][2]+Sum0[1][3]+Sum0[1][4]+Sum0[1][5]+Sum0[1][6]+Sum0[1][7];
                Sum[2] += Sum0[2][0]+Sum0[2][1]+Sum0[2][2]+Sum0[2][3]+Sum0[2][4]+Sum0[2][5]+Sum0[2][6]+Sum0[2][7];
                
            for (int i = (gv.length-2*deta)/UNROLLING*UNROLLING+deta;i<(int)gv.length-(int)deta; ++i){
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

    return o;
}

Image apply_gb(Image a, FVec gv)
{
    // struct timeval start_time, stop_time, elapsed_time; 
    // gettimeofday(&start_time,NULL);

    Image r,g,b;
    get_RGB(a,&r,&g,&b);
    
    // gettimeofday(&stop_time,NULL);
    // timersub(&stop_time, &start_time, &elapsed_time); 
    // printf("RGB change time1: %f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    // gettimeofday(&start_time,NULL);

    Image o = gb_h(a,r,g,b,gv);

    // gettimeofday(&stop_time,NULL);
    // timersub(&stop_time, &start_time, &elapsed_time); 
    // printf("gb_h time: %f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    // gettimeofday(&start_time,NULL);

    Image test = transpose(o);
    // gettimeofday(&stop_time,NULL);
    // timersub(&stop_time, &start_time, &elapsed_time); 
    // printf("transpose time: %f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);


    // gettimeofday(&start_time,NULL);

    Image rt,gt,bt;
    get_RGB(test,&rt,&gt,&bt);

    // gettimeofday(&stop_time,NULL);
    // timersub(&stop_time, &start_time, &elapsed_time); 
    // printf("RGB change time2: %f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    // gettimeofday(&start_time,NULL);

    test = gb_h(test,rt,gt,bt,gv);

    // gettimeofday(&stop_time,NULL);
    // timersub(&stop_time, &start_time, &elapsed_time); 
    // printf("gb_v: %f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    // gettimeofday(&start_time,NULL);

    Image c = transpose(test);

    // gettimeofday(&stop_time,NULL);
    // timersub(&stop_time, &start_time, &elapsed_time); 
    // printf("transpose time: %f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    free(b.data);
    free(test.data);
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

    // for (int i=0;i<19;i++){
    //     FVec v = make_gv(a, x0, x1, dim, min_dim);
    //     Image img;
    //     img.data = stbi_loadf(argv[1], &(img.dimX), &(img.dimY), &(img.numChannels), 0);
    //     Image imgOut = apply_gb(img, v);
    //     stbi_write_jpg(argv[2], imgOut.dimX, imgOut.dimY, imgOut.numChannels, imgOut.data, 90);
    //     printf("-------------------\n");
    // }

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
    printf("%f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0f);
    free(imgOut.data);
    free(v.data);
    free(v.sum);
    return 0;
}