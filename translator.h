#ifndef TRANSLATOR_H
#define TRANSLATOR_H

#define MAX_CMD_NUMBER  50  /*max number of lines each input file*/

int translate(const char*in, const char*out);

extern uint32_t cmd_data[];
extern uint32_t cmd_num;

#endif