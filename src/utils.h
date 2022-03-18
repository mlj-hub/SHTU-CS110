#ifndef UTILS_H
#define UTILS_H
#include <stdint.h>

#define MAX_CMD_NUMBER  50  /*max number of lines each input file*/

/*function declare*/
void input_handler(FILE * input_file);

extern uint32_t cmd_data[];
extern uint32_t cmd_num;
#endif