#ifndef UTILS_H
#define UTILS_H
#include <stdint.h>
#include "compression.h"
#define MAX_CMD_NUMBER  50  /*max number of lines each input file*/

/*function declare*/
uint32_t * input_handler(FILE * input_file);
void output_handler(cmd_info_t *);

extern uint32_t cmd_num;
#endif