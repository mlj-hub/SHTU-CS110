#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "compression.h"
#include "utils.h"

void R_compress(cmd_info_t * cmd_info){
    (void) cmd_info;
}

cmd_state_t R_check(uint32_t cmd){
    (void) cmd;
    return 0;
}

void handle_compressible(cmd_info_t * cmd_info){
    (void) cmd_info;
}

void handle_unsure(cmd_info_t * cmd_info){
    (void) cmd_info;
}

void format_compressible_check(cmd_info_t * cmd_info,uint32_t * cmd_data){
    uint32_t i =0;
    for(i=0;i<cmd_num;i++){
        cmd_info[i].cmd = cmd_data[i];
        switch(cmd_data[i] & OPCODE){
            case R:
                cmd_info[i].format = R;
                cmd_info[i].state = R_check(cmd_info[i].cmd);
                break;
            case I:
                cmd_info[i].format = I;
                break;
            case LI:
                cmd_info[i].format = LI;
                break;
            case U:
                cmd_info[i].format = U;
                break;
            case S:
                cmd_info[i].format = S;
                break;
            case SB:
                cmd_info[i].format = SB;
                break;
            case UJ:
                cmd_info[i].format = UJ;
                break;
            case JI:
                cmd_info[i].format = JI;
                break;
            default:
                cmd_info[i].format = UNVAILD;
        }
    }
    return ;
}

cmd_info_t * compress(uint32_t * cmd_data){
    cmd_info_t * cmd_info;
    cmd_info = (cmd_info_t *)calloc(sizeof(cmd_info_t),cmd_num);

    format_compressible_check(cmd_info,cmd_data);
    free(cmd_data);
    return NULL;
}

/* Your code here... */

