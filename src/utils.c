#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

uint32_t cmd_num = 0; /*indicate number of the input cmd*/

/*handle the input file*/
uint32_t * input_handler(FILE * input_file){
    /*uint8_t flag = 0;*/
    uint32_t * cmd_data = (uint32_t *)calloc(sizeof(uint32_t),MAX_CMD_NUMBER);  /*array which store the cmd in uint32_t form*/
    /*used to store the character read from the input file*/
    char data=0;  
    /*i for bit shift operation and index indicate the writing position of cmd_data*/
    uint32_t i,index; 
    /*used to store the cmd*/
    uint32_t cmd=0;
    i=31;index=0;
    /*loop to read in all cmd */
    while(data!=EOF){
        data = fgetc(input_file);
        /*handle with \n, space, \t or etc */
        if(data!='0' && data!='1'){
            if(cmd)
            /*store cmd*/
                cmd_data[index++]=cmd;
            /*reset i*/
            i=31;
            /*rest cmd*/
            cmd = 0;
        }
        /*read in bits, shift them then adding to cmd*/
        else{
            cmd|=((data-'0')<<i);
            i--;
        }
        /*read character from file*/
    }
    /*update the number of cmd*/
    cmd_num = index;
    return cmd_data;
}

void output_handler(cmd_info_t  * cmd_info,FILE * output_file){
    /*initialize parameters*/
    uint32_t i = 0;
    uint32_t cmd;
    /*transerve the cmd*/
    for(i=0;i<cmd_num;i++){
        /*initialize the parameters*/
        uint32_t j=0;
        uint8_t bit;
        /*get cmd*/
        cmd = cmd_info[i].cmd;
        /*if imcompressible, write 32bits*/
        if(cmd_info[i].state !=COMPRESSIBLE){
            for(j=0;j<32;j++){
                /*left shift j bits*/
                bit = ((cmd<<j)&HIGHEST_BIT)>>31;
                fputc(bit+'0',output_file);
            }
            /*a new line*/
            fputc('\n',output_file);
        }
        /*if compressible. write 16 bits*/
        else{
            /*cmd shift left 16*/
            cmd<<=16;
            for(j=0;j<16;j++){
                /*shift left j*/
                bit = ((cmd<<j)&HIGHEST_BIT)>>31;
                fputc(bit+'0',output_file);
            }
            /*a new line*/
            fputc('\n',output_file);
        }
    }
    /*free cmd_info*/
    free(cmd_info);
}