#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

uint32_t cmd_num = 0; /*indicate number of the input cmd*/

/*handle the input file*/
uint32_t * input_handler(FILE * input_file){
    uint32_t * cmd_data = (uint32_t *)calloc(sizeof(uint32_t),MAX_CMD_NUMBER);  /*array which store the cmd in uint32_t form*/
    char data;  /*used to store the character read from the input file*/
    uint32_t i,index; /*i for bit shift operation and index indicate the writing position of cmd_data*/
    uint32_t cmd=0;/*used to store the cmd*/
    i=31;index=0;
    data = fgetc(input_file);
    /*loop to read in all cmd */
    while(data!=EOF){
        /*handle with \n, space, \t or etc */
        if(data!='0' && data!='1'){
            if(cmd)
                cmd_data[index++]=cmd;
            i=31;
            cmd = 0;
        }
        /*read in bits, shift them then adding to cmd*/
        else{
            cmd|=((data-'0')<<i);
            i--;
        }
        /*read character from file*/
        data = fgetc(input_file);
    }
    /*update the number of cmd*/
    cmd_num = index;
    return cmd_data;
}

void output_handler(cmd_info_t  * cmd_info){
    (void ) cmd_info;
}