#ifndef _GAME_CONFIG_H_
#define _GAME_CONFIG_H_
#include "lcd/lcd.h"
#include <stdint.h>

#define TREX_START_Y  40
#define TREX_END_Y  59
#define GROUND_START_Y 60
#define OBSTACLE_START_Y 40
#define OBSTACLE_END_Y  59
#define SCORE_START_X  150
#define SCORE_START_Y 10

typedef enum{
    MOVE0,
    MOVE1,
    SQUAT0,
    SQUAT1,
    JUMP
}trex_state_e;

typedef struct{
    trex_state_e state;
    u16 start_x;
    u16 start_y;
    u16 end_x;
    u16 end_y;
}trex_t;

typedef enum{
    INIT,
    PLAYING,
    OVER,
}game_state_e;

typedef struct{
    int32_t score;
    game_state_e state;
}game_t;

#endif