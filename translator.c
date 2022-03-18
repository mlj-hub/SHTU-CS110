/*  Project 1.1: RISC-V instructions to RISC-V compressed instructions in C89.
    The following is the starter code provided for you. To finish the task, you 
    should define and implement your own functions in translator.c, compression.c, 
    utils.c and their header files.
    Please read the problem description before you start.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "src/compression.h"
#include "src/utils.h"

#include "translator.h"

uint32_t cmd_num = 0; /*indicate number of the input cmd*/
uint32_t cmd_data[MAX_CMD_NUMBER]; /*array which store the cmd in uint32_t form*/


/*handle the input file*/
void input_handler(FILE * input_file){
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
}

/*check if file can be correctly opened */
static int open_files(FILE** input, FILE** output, const char* input_name, const char* output_name){ 
    *input = fopen(input_name, "r");
    if (!*input){ /* open input file failed */
        printf("Error: unable to open input file: %s\n", input_name);
        return -1;
    }

    *output = fopen(output_name, "w");
    if (!*output){ /* open output file failed */
        printf("Error: unable to open output file: %s\n", output_name);
        fclose(*input);
        return -1;
    }
    return 0; /* no problem opening files */
}

static int close_files(FILE** input, FILE** output){
    fclose(*input);
    fclose(*output); /* close the files at the end */
    return 0;
}

static void print_usage_and_exit() {
    printf("Usage:\n");
    printf("Run program with translator <input file> <output file>\n"); /* print the correct usage of the program */
    exit(0);
}


/*Run the translator 
*/
int translate(const char*in, const char*out){
    FILE *input, *output;
    int err = 0;
    if (in){    /* correct input file name */
        if(open_files(&input, &output, in, out) != 0)
            exit(1);
        /* handle input file */
        input_handler(input);
        /* write correct result to the output file */
        /* ... */
        close_files(&input, &output);
    }
    return err;
}

/* main func */
int main(int argc, char **argv){
    char* input_fname, *output_fname;
    int err;

    if (argc != 3) /* need correct arguments */
        print_usage_and_exit();

    input_fname = argv[1];
    output_fname = argv[2];

    err = translate(input_fname, output_fname); /* main translation process */
    if (err)
        printf("One or more errors encountered during translation operation.\n"); /* something wrong */
    else
        printf("Translation process completed successfully.\n"); /* correctly output */

    return 0;
}