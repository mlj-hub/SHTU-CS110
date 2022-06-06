#include "lcd/lcd.h"
#include <string.h>
#include "utils.h"
#include "img.h"
#include "game_config.h"

game_t game;
trex_t trex;

void Inp_init(void)
{
    gpio_init(GPIOA, GPIO_MODE_IN_FLOATING, GPIO_OSPEED_50MHZ, GPIO_PIN_8);
}

void Adc_init(void) 
{
    rcu_periph_clock_enable(RCU_GPIOA);
    gpio_init(GPIOA, GPIO_MODE_AIN, GPIO_OSPEED_50MHZ, GPIO_PIN_0|GPIO_PIN_1);
    RCU_CFG0|=(0b10<<14)|(1<<28);
    rcu_periph_clock_enable(RCU_ADC0);
    ADC_CTL1(ADC0)|=ADC_CTL1_ADCON;
}

void IO_init(void)
{
    Inp_init(); // inport init
    Adc_init(); // A/D init
    Lcd_Init(); // LCD init
}

void update_ground(){

}

void update_trex(){

}

void update_obstacle(){

}

void update_score(){
    int cnt=1;
    int32_t temp=game.score;
    while(temp%10!=temp){
        cnt++;
        temp/=10;
    }
    LCD_ShowNum(SCORE_START_X-cnt*8,SCORE_START_Y,game.score,cnt,WHITE);
}

void update(){
    
}

void get_state(){

}



int main(void)
{
    IO_init();         // init OLED
    // YOUR CODE HERE
    LCD_Clear(BLACK);

    int g1_start_x=0,g2_start_x=160;
    int show_trex = 0;
    int32_t score=0;
    delay_1ms(100);
    while(1){

        if(show_trex==8){
            LCD_ShowPic(10,40,29,59,trex1);
        }
        if(show_trex==16){
            LCD_ShowPic(10,40,29,59,trex2);
        }
        if(g1_start_x==0){
            LCD_ShowPicture_Part(0,GROUND_START_Y,g1,160-g2_start_x,0,160,10,160);
            LCD_ShowPicture_Part(g2_start_x,GROUND_START_Y,g2,0,0,160-g2_start_x,10,160);
            g2_start_x--;
            if(g2_start_x==0)
                g1_start_x = 160;
        }
        if(g2_start_x==0){
            LCD_ShowPicture_Part(0,GROUND_START_Y,g2,160-g1_start_x,0,160,10,160);
            LCD_ShowPicture_Part(g1_start_x,GROUND_START_Y,g1,0,0,160-g1_start_x,10,160);
            g1_start_x--;
            if(g1_start_x==0)
                g2_start_x = 160;
        }
        if(show_trex<16)
            show_trex++;
        else
            show_trex=0;
        delay_1ms(50);
        game.score++;
        update_score();
    }
}
