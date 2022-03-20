#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "compression.h"
#include "utils.h"
/*right*/
void CR_compress(cmd_info_t * cmd_info){
    /*initiaize of the parameters in CR */
    uint32_t CR_funct4 = 0;
    uint32_t CR_rdORrs1 = 0;
    uint32_t CR_rs2 = 0;
    uint32_t CR_op= 0;
    /*if the RISCV32 command is R,then in CR_compress the command is [add] */
    if (cmd_info->format == R)
    {
        /*get the parameters of the old value according to the format*/
        uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
        uint32_t Old_rs1 = (cmd_info->cmd>>15)&REGISTER;
        uint32_t Old_rs2 = (cmd_info->cmd>>20)&REGISTER;
        /*c.add*/
        if (Old_rd == Old_rs1)
        {
            /*the cmd is [c.add]*/
            CR_funct4 = 0x9;
            CR_rdORrs1 = Old_rd;
            /*set rs2*/
            CR_rs2 = Old_rs2;
            CR_op = 0x2;
        }
        /*c.mv*/
        else
        {
            /*get parameters*/
            CR_funct4 = 0x8;
            CR_rdORrs1 = Old_rd;
            /*set rs2*/
            CR_rs2 = Old_rs2;
            CR_op = 0x2;
            /* the cmd is [c.mv]*/
        }
    }
    /* if the RISCV32 command is I,opcode is JI then in CR_compress the command is [jalr]*/
    else if(cmd_info->format == JI)
    {
        /*initialize the parameters*/
        uint32_t Old_rs1 = (cmd_info->cmd>>15)&REGISTER;
        uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
        if (Old_rd == 0x0)
        {
            /* the cmd is [c.jr] */
            CR_funct4 = 0x8;
            CR_rdORrs1 = Old_rs1;
            /*set rs2*/
            CR_rs2 = 0;
            CR_op = 0x2;
        }
        else
        {
            /*get parameters*/
            CR_funct4 = 0x9;
            CR_rdORrs1 = Old_rs1;
            /*set rs2*/
            CR_rs2 = 0;
            CR_op = 0x2;
            /* the cmd is [c.jalr] */
        }
    }
    /*change the 32-bit code in the cmd_info*/
    cmd_info->cmd = CR_op + (CR_rs2<<2) + (CR_rdORrs1<<7)+(CR_funct4<<12);
    return;

}
/*right*/
void CS_compress(cmd_info_t * cmd_info)
{
    /*if the RISCV32 command is S,then in CS_compress the command is [sw] */
    if (cmd_info->format == S)
    {
        /*parameters of the RISCV-16  Command CS_T1*/
        uint32_t CS_T1_funct3 = 0;
        uint32_t CS_T1_IMM3 = 0;
        uint32_t CS_T1_RS1 = 0;
        /*get imm2*/
        uint32_t CS_T1_IMM2 = 0;
        uint32_t CS_T1_RS2 = 0;
        uint32_t CS_T1_OP = 0;
        /*parameters of the RISCV-32 Command CS_T2*/
        uint32_t Old_imm12 = 0;
        uint32_t Old_imm5 = 0;
        uint32_t Old_imm7 = 0;
        /*get the parameters of old cmd */
        uint32_t Old_rs1 = (cmd_info->cmd>>15)&REGISTER;
        uint32_t Old_rs2 = (cmd_info->cmd>>20)&REGISTER;
        /*get imm from two part of the cmd*/
        Old_imm5 = (cmd_info->cmd>>7)&IMM5;
        Old_imm7= ( ((int32_t)cmd_info->cmd) >>25 )&SIGN_ALL;
        Old_imm12 = (Old_imm7<<5)|Old_imm5;
        /*get imm from two part of the cmd*/
        /* the  cmd is [c.sub]*/
        CS_T1_funct3 = 6;
        CS_T1_RS1 = Old_rs1-8;
        CS_T1_RS2 = Old_rs2-8;
        /*get the offset[5:3] of given imm*/
        CS_T1_IMM3 = (Old_imm12 >> 3) & 0x7;
        /*get the offset[2|6] of given imm*/
        CS_T1_IMM2 += (((Old_imm12 >> 2) & 0x1) <<1) ;
        CS_T1_IMM2 += ((Old_imm12 >> 6) & 0x1);
        /* change the cmd*/
        cmd_info->cmd = (CS_T1_OP) + (CS_T1_RS2 << 2) + (CS_T1_IMM2 << 5) + \
                        (CS_T1_RS1<<7) + (CS_T1_IMM3 << 10) + (CS_T1_funct3<< 13);
        return;

    }
    /*if the RISCV32 command is R,then in CS_compress the command into CS Type2 */
    else if (cmd_info->format == R)
    {
        /*parameters of the RISCV-16 CS Command [c.sw]*/
        uint32_t Old_funct3 = (cmd_info->cmd>>12)&FUNCT3;
        uint32_t Old_rs1 = (cmd_info->cmd>>15)&REGISTER;
        uint32_t Old_rs2 = (cmd_info->cmd>>20)&REGISTER;
        /*get funct7*/
        uint32_t Old_funct7 = (cmd_info->cmd>>25)&FUNCT7;
        /*get parameters from the green card*/
        uint32_t CS_T2_funct6 = 0x23;
        uint32_t CS_T2_RDorRS1 = Old_rs1 - 8;
        uint32_t CS_T2_funct2 = 0;
        /*set registers to the required ones*/
        uint32_t CS_T2_RS2 = Old_rs2 - 8;
        uint32_t CS_T2_OP = 1;
        /*c.add or c.or or c.xor*/
        if (Old_funct7 == 0)
        {
            /* the RISCV16 cmmmand is [c.and] [c.or] [c.xor]*/
            switch(Old_funct3)
            {
                case 7:
                /* and */
                CS_T2_funct2 = 3;
                break;
                case 6:
                CS_T2_funct2 = 2;
                /* or */
                break;
                case 4:
                CS_T2_funct2 = 1;
                /* xor*/
                break;
            }
        }
        else
        {
            /* the RISCV16 command is [c.sub] */
            CS_T2_funct2 = 0;
        }
        /*change the value of cmd*/
        cmd_info->cmd = CS_T2_OP + (CS_T2_RS2<<2) + \
            (CS_T2_funct2<<5) + (CS_T2_RDorRS1<<7) + \
            (CS_T2_funct6<<10);
    }



}
/*right*/
void CL_compress(cmd_info_t *cmd_info)
{
    /* the parameter of I*/
    uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
    uint32_t Old_rs1 = (cmd_info->cmd>>15)&REGISTER;
    int32_t  Old_imm12 = (cmd_info->cmd>>20)&IMM12;
    /* the parameter of CL*/
    uint32_t CL_funct3 = 2;
    uint32_t CL_IMM3 = (Old_imm12>>3) & 0x7;
    /*set registers to the requided ones*/
    uint32_t CL_RS1 = Old_rs1 - 8;
    uint32_t CL_IMM2 = 0;
    /*set registers to the requided ones*/
    uint32_t CL_RD= Old_rd - 8;
    uint32_t CL_OP = 0;
    /*get imm according to the green card*/
    CL_IMM2 += (((Old_imm12 >> 2) & 0x1) <<1) ;
    CL_IMM2 += ((Old_imm12 >> 6) & 0x1);
    cmd_info->cmd = CL_OP + (CL_RD << 2) + (CL_IMM2<<5) +(CL_RS1 << 7) + (CL_IMM3 << 10) + (CL_funct3 << 13);
}
/*check add,and,or,xor,sub*/  /*right*/
void CI_compress(cmd_info_t *cmd_info)
{
    /*CI parameters*/
    uint32_t CI_funct3 = 0;
    uint32_t CI_IMM1 = 0;
    uint32_t CI_RS1orRD = 0 ;
    /*set imm5*/
    uint32_t CI_IMM5 = 0;
    uint32_t CI_OP = 0;
    /*the input command is lui,the old type is U*/
    if(cmd_info->format == U)
    {
        /*change the value of CI ,the commmand is [c.LUI]*/
        uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
        int32_t Old_imm20 = (cmd_info->cmd>>12)&IMM20;
        /*get parameters according to the green card*/
        CI_funct3 = 3;
        CI_IMM1 = (Old_imm20 >> 5) & 0x1;
        /*set parameters according to the format*/
        CI_RS1orRD = Old_rd;
        CI_IMM5 = Old_imm20& IMM5;
        CI_OP = 1;
    }
    /* the old type is I,the input is addi or slli*/
    else if (cmd_info->format == I)
    {
        /*get parameters according to the format*/
        uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
        uint32_t Old_rs1 = (cmd_info->cmd>>15)&REGISTER;
        /*set funct3*/
        uint32_t Old_funct3 = (cmd_info->cmd>>12)&FUNCT3;
        int32_t  Old_imm12 = (cmd_info->cmd>>20)&IMM12;
        /*addi*/
        if (Old_funct3 == 0)
        {
            /* the old cmd is addi */
            CI_IMM1 = (Old_imm12>>5) & 0x1;
            CI_RS1orRD = Old_rd;
            /*get imm5*/
            CI_IMM5 =  Old_imm12 & IMM5;
            CI_OP = 1;
            /*c.li*/
            if (Old_rs1 == 0)
                /* the commmand is [c.li]*/
                CI_funct3 = 2;
            else
                /* the command is [c.slli]*/
                CI_funct3 = 0;
        }
        else
        {
            /* the old cmd is slli*/
            CI_funct3 = 0;
            /*shamt[5] should be zero, we directly assigned 0 to imm1*/
            CI_IMM1 = 0;
            CI_RS1orRD = Old_rd;
            /*get imm5 from the old parameters*/
            CI_IMM5 =  Old_imm12 & IMM5;
            CI_OP = 2;

        }
    }
    /*changt the value of the command*/
    cmd_info->cmd = CI_OP + (CI_IMM5 <<2) + (CI_RS1orRD << 7) +(CI_IMM1<<12) + (CI_funct3 << 13);
    return;

}
/*right*/
void CB_T2_Compress(cmd_info_t *cmd_info)
{
        /*CI parameters*/
    uint32_t CB_T2_funct3 = 0x4;
    uint32_t CB_T2_IMM1 = 0;
    uint32_t CB_T2_funct2 = 0;
    /*initilaize the parameters*/
    uint32_t CB_T2_RS1orRD = 0 ;
    uint32_t CB_T2_IMM5 = 0;
    uint32_t CB_T2_OP = 1;
    /* given parameters,type I*/
    uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
    uint32_t Old_funct3 = (cmd_info->cmd>>12)&FUNCT3;
    /*get imm*/
    uint32_t Old_imm12 = (cmd_info->cmd>>20)&IMM12;
    /*get shamt*/
    uint32_t Old_shamt = (cmd_info->cmd>>20)&REGISTER;
    uint32_t Old_funct7 = (cmd_info->cmd>>25)&FUNCT7;
    if (Old_funct3 == 0x5)
    {
        /*c.srli*/
        if (Old_funct7 == 0)
        {
            /* shamt 5 should be 0*/
            CB_T2_IMM1 = 0;
            CB_T2_funct2 = 0;
            /*change the value of Rd*/
            CB_T2_RS1orRD = Old_rd-8 ;
            CB_T2_IMM5 = Old_shamt;
            /* the RISCV 16 cmd is [c.srli]*/
        }
        else
        {
            /* shamt 5 should be 0*/
            CB_T2_IMM1 = 0;
            CB_T2_funct2 = 1;
            /*change the value of Rd*/
            CB_T2_RS1orRD = Old_rd-8 ;
            CB_T2_IMM5 = Old_shamt;
            /* the RISCV 16 cmd is [c.srai]*/
        }
    }
    else
    {   
            /* shamt 5 should be 0*/
            CB_T2_IMM1 = (Old_imm12>>5) & 0x1;
            CB_T2_funct2 = 2;
            /*change the value of Rd*/
            CB_T2_RS1orRD = Old_rd-8 ;
            CB_T2_IMM5 = Old_imm12 & IMM5;
        /*the RISCV16 cmd is [c.andi]*/
    }
    cmd_info->cmd = CB_T2_OP + (CB_T2_IMM5<<2) + (CB_T2_RS1orRD << 7) +(CB_T2_funct2<< 10) + (CB_T2_IMM1<<12)+(CB_T2_funct3 << 13);
}
/*R_Check*/ /*right*/
void R_check(cmd_info_t * cmd_info){
    uint32_t cmd = cmd_info->cmd;
    /*get parameters according to the format*/
    uint32_t rd = (cmd>>7)&REGISTER;
    uint32_t funct3 = (cmd>>12)&FUNCT3;
    /*get register*/
    uint32_t rs1 = (cmd>>15)&REGISTER;
    uint32_t rs2 = (cmd>>20)&REGISTER;
    uint32_t funct7 = (cmd>>25)&FUNCT7;
    /*according to funct3 to check what cmd it actually is*/
    switch(funct3){
        /*add or sub*/
        case 0x0:
            /*add*/
            if(funct7 == 0x0){
                /*c.add and c.mv*/
                if((rd == rs1 || !rs1) && rd && rs2){
                    /*conditions if compressible*/
                    cmd_info->state = COMPRESSIBLE;
                    cmd_info->c_format = CR;
                    return ;
                }
            }
            /*c.sub*/
            else{
                if(rd==rs1 && rd>=8 && rd<=15){
                    /*conditions if compressible*/
                    cmd_info->state = COMPRESSIBLE;
                    cmd_info->c_format = CS_T2;
                    /*return*/
                    return;
                }
            }
            /*cannot be compressed*/
            cmd_info-> state = INCOMPRESSIBLE;
            break;
        /*xor,and,or share the same condition when compressible*/
        case 0x4:
        case 0x6:
        case 0x7:
            /*compressible only when rd is same with rs1*/
            if(rd==rs1 && rd<=15 && rd>=8){
                cmd_info -> state = COMPRESSIBLE;
                cmd_info -> c_format = CS_T2;
                /*return */
                return ;
            }
            /*otherwise, incompressible*/
            cmd_info->state = INCOMPRESSIBLE;
            break;
        default:
            /*otherwise, incompressible*/
            cmd_info->state = INCOMPRESSIBLE;
    }
}
/*right*/
void I_check(cmd_info_t * cmd_info){
    /*initialize the parameters*/
    uint32_t cmd = cmd_info -> cmd;
    uint32_t rd = (cmd>>7)&REGISTER;
    uint32_t rs1 = (cmd>>15)&REGISTER;
    /*get funct3*/
    uint32_t funct3 = (cmd>>12)&FUNCT3;
    int32_t  imm12 = ((int32_t)cmd>>20)&SIGN_ALL;
    switch(funct3){
        /*addi*/
        case 0x0:
            /*c.li rd!=0,rs1=0,-32<=imm12<=31*/
            if(rd && !rs1 && imm12<=31 && imm12>=-32){
                /*conditions when compressible*/
                cmd_info->state = COMPRESSIBLE;
                cmd_info->c_format = CI;
                return;
            }
            /*c.addi rd=rs1!=0, imm!=0*/
            else if(rd==rs1 && rd && imm12 && imm12<=31 && imm12>=-32){
                /*conditions when compressible*/
                cmd_info->state = COMPRESSIBLE;
                cmd_info->c_format = CI;
                return;
            }
            /*otherwise it is incompressible*/
            cmd_info ->state = INCOMPRESSIBLE;
            break;
        /*slli*/
        case 0x1:
            if(((cmd>>20)&0x20) ==0 && rd==rs1 && rd!=0){
                /*conditions when compressible*/
                cmd_info -> c_format = CI;
                cmd_info -> state = COMPRESSIBLE;
                return;
            }
            /*otherwise, incompressible*/
            cmd_info->state = INCOMPRESSIBLE;
            break;
        /*srli or srai*/
        case 0x5:
            if(rd == rs1 && ((cmd>>20)&0x20) ==0 && rd>=8 && rd<=15){
                /*conditions when compressible*/
                cmd_info -> c_format = CB_T2;
                cmd_info -> state = COMPRESSIBLE;
                return;
            }
            /*otherwise, incompressible*/
            cmd_info->state = INCOMPRESSIBLE;
            break;
        /*andi*/
        case 0x7:
            /*c.andi rd=rs1, -32<=imm12<=31*/
            if(rd==rs1 && -32<=imm12 && imm12<=31 && rd>=8 && rd<=15){
                /*conditions when compressible*/
                cmd_info->state = COMPRESSIBLE;
                cmd_info->c_format = CB_T2;
                return;
            }
            /*otherwsie, incompressible*/
            cmd_info -> state = INCOMPRESSIBLE;
            break;
        default:
            /*otherwise, incompressible*/
            cmd_info -> state = INCOMPRESSIBLE;
    }
}

/*check lw*/ /*right*/
void LI_check(cmd_info_t * cmd_info){
    uint32_t cmd = cmd_info->cmd;
    /*get parameters according to format*/
    uint32_t rd = (cmd>>7)&REGISTER;
    uint32_t rs1 = (cmd>>15)&REGISTER;
    int32_t  imm12 = ((int32_t)cmd>>20)&SIGN_ALL;
    uint32_t funct3 = (cmd>>12)&FUNCT3;
    /*check whether it is lw*/
    if(funct3 != 0x2)
        cmd_info->state = INCOMPRESSIBLE;
    /*check whether rs1 is between x8-x15*/
    else if(rs1>=16 || rs1<=7)
        cmd_info->state = INCOMPRESSIBLE;
    /*check whether rd is between x8-x15*/
    else if(rd>=16 || rd<=7)
        cmd_info->state = INCOMPRESSIBLE;
    /*incompressible if offset is negative  */
    else if(imm12 < 0)
        cmd_info->state = INCOMPRESSIBLE;
    /*incompressible if the least two significant bits is non-zero*/
    else if(imm12 % 4!=0)
        cmd_info->state = INCOMPRESSIBLE;
    /*incompressible if imm is larger than 64*/
    else if ((imm12>>2)>=32)
        cmd_info->state = INCOMPRESSIBLE;
    /*otherwise compressible*/
    else{
        cmd_info->state = COMPRESSIBLE;
        cmd_info->c_format = CL;
    }
}

/*check sw*/
void S_check(cmd_info_t * cmd_info){
    uint32_t cmd = cmd_info ->cmd;
    /*get parameters according to format*/
    int32_t imm5=0;
    int32_t imm12 = 0;
    int32_t imm7=0;
    /*get funct3 and register*/
    uint32_t funct3 = (cmd>>12)&FUNCT3;
    uint32_t rs1 = (cmd>>15)&REGISTER;
    uint32_t rs2 = (cmd>>20)&REGISTER;
    /*get imm from two part of the cmd*/
    imm5 = (cmd>>7)&IMM5;
    imm7=((int32_t)cmd>>25)&SIGN_ALL;
    imm12 = (imm7<<5)|imm5;
    /*check whether the cmd is sw*/
    if(funct3 != 0x2)
        cmd_info->state = INCOMPRESSIBLE;
    /*check whether rs1 is between x8-x15*/
    else if(rs1>=16 || rs1<=7)
        cmd_info->state = INCOMPRESSIBLE;
    /*check whether rs2 is between x8-x15*/
    else if(rs2>=16 || rs2<=7)
        cmd_info->state = INCOMPRESSIBLE;
    /*incompressible if offset is negative*/
    else if(imm12 < 0)
        cmd_info->state = INCOMPRESSIBLE;
    /*incompressible if the least two significant bits is non-zero*/
    else if(imm12 % 4!=0)
        cmd_info->state = INCOMPRESSIBLE;
    /*incompressible if imm is larger than 64*/
    else if ((imm12>>2)>=32)
        cmd_info->state = INCOMPRESSIBLE;
    /*otherwise compressible*/
    else{
        cmd_info->state = COMPRESSIBLE;
        cmd_info ->c_format = CS_T1;
    }
}
/*bne,beq*/
void SB_check(cmd_info_t * cmd_info){
    /*get parameters*/
    uint32_t cmd = cmd_info->cmd;
    uint32_t funct3 = (cmd>>12)&REGISTER;
    /*get register*/
    uint32_t rs1 = (cmd>>15)&REGISTER;
    uint32_t rs2 = (cmd>>20)&REGISTER;
    /*funct3 larger than 2, incompressible*/
    if(funct3 >=2)
        cmd_info->state = B_J_INCOMPRESSIBLE;
    /*register not in the range*/
    else if(rs1>=16 || rs1<=7)
        cmd_info->state = B_J_INCOMPRESSIBLE;
    /*rs2 is not zero*/
    else if(rs2!=0)
        cmd_info->state = B_J_INCOMPRESSIBLE;
    /*otherwise, compressible*/
    else{
        cmd_info->state = COMPRESSIBLE;
        /*set compressed format*/
        cmd_info->c_format = CB_T1;
    }
}

/*check jalr*/
void JI_check(cmd_info_t * cmd_info){
    uint32_t cmd = cmd_info -> cmd;
    /*get parameters according to the format*/
    uint32_t funct3 = (cmd>>12)&FUNCT3;
    uint32_t rd = (cmd>>7)&REGISTER;
    /*get rs1 and imm12*/
    uint32_t rs1 = (cmd>>15)&REGISTER;
    uint32_t  imm12 = (cmd>>20)&IMM12;
    /*check funct3 to determine whether the command is jalr*/
    if(funct3 != 0x0)
        cmd_info->state = INCOMPRESSIBLE;
    /*incompressible if imm is non-zero or rs1 is x0 or rd is not x0 or x1 */
    else if(imm12 || !rs1 || rd>=2)
        cmd_info->state =INCOMPRESSIBLE;
    else{
        /*otherwise, compressible*/
        cmd_info->state= COMPRESSIBLE;
        cmd_info->c_format = CR;
    }
}

void U_check(cmd_info_t * cmd_info){
    /*since there two U_type commands and they are different in the opcode*/
    uint32_t cmd = cmd_info -> cmd;
    uint32_t rd = (cmd >>7)&REGISTER;
    /*get imm20*/
    int32_t imm20 = cmd&0xfffff000;
    /*highest bits are all one and 17 bit is 1*/
    if((imm20&0xfffc0000) == 0xfffc0000 && (imm20&0x00020000)!=0 && rd!=0 && rd!=2 && imm20){
        /*set c_format*/
        cmd_info->c_format = CI;
        /*compressible*/
        cmd_info->state = COMPRESSIBLE;
        return;
    }
    /*highest bits are all 0 and 17 bit is 0*/
    else if((imm20&0xfffc0000) == 0 && (imm20&0x00020000)==0 && rd!=0 && rd!=2 && imm20){
        /*set c_format*/
        cmd_info->c_format = CI;
        /*compressible*/
        cmd_info->state = COMPRESSIBLE;
        return;
    }
    /*otherwise, incompressible*/
    cmd_info -> state =INCOMPRESSIBLE;
}
/*jar*/
void UJ_check(cmd_info_t * cmd_info){
    /*get parameters*/
    uint32_t cmd = cmd_info->cmd;
    uint32_t rd = (cmd>>7)&REGISTER;
    /*if rd is not zero or x1*/
    if(rd>=2)
        cmd_info->state = B_J_INCOMPRESSIBLE;
    else{
        /*otherwise, compressible*/
        cmd_info->state = COMPRESSIBLE;
        /*set c_format*/
        cmd_info ->c_format = CJ;
    }
}

void handle_compressible(cmd_info_t * cmd_info){
    /*initilize */
    uint32_t i=0;
    for(i=0;i<cmd_num;i++){
        /*if incompressible, continue*/
        if(cmd_info[i].state != COMPRESSIBLE)
            continue;
        /*get c_format*/
        switch(cmd_info[i].c_format){
            /*CR*/
            case CR:
                CR_compress(&cmd_info[i]);
                break;
            /*CI*/
            case CI:
                CI_compress(&cmd_info[i]);
                break;
            /*CL*/
            case CL:
                CL_compress(&cmd_info[i]);
                break;
            /*CS*/
            case CS_T1: case CS_T2:
                CS_compress(&cmd_info[i]);
                break;
            /*CB_T2*/
            case CB_T2:
                CB_T2_Compress(&cmd_info[i]);
                break;
            /*otherwise*/
            default:
                continue;
        }
    }
}

void get_uj_offset(uint32_t imm20,int32_t * offset){
    *offset |=((imm20&0x7fe00)>>9); /*offset 10:1*/
    *offset |=((imm20&0x100)<<2); /*offset 11*/
    *offset |=((imm20&0xff)<<11); /*offset 19:12*/
    *offset |=((imm20&0x80000));/*offset 20*/
    /*set the rightest bit to 0*/
    *offset<<=1;
}

void handle_unsure(cmd_info_t * cmd_info){
    /*initialize*/
    uint32_t i=0;
    /*tranverse the cmd*/
    for(i=0;i<cmd_num;i++){
        uint32_t n_cmd = 0x1;
        uint32_t cmd = cmd_info[i].cmd;
        /*if incompressible, continue*/
        if(cmd_info[i].state==INCOMPRESSIBLE)
            continue;
        /*c.j and c.jal*/
        if(cmd_info[i].format == UJ){
            /*get rd and imm20*/
            uint32_t rd = (cmd>>7)&REGISTER;
            uint32_t imm20 = ((int32_t)cmd>>12)&SIGN_ALL;
            /*get offset*/
            int32_t offset = imm20&0xfff00000;
            /*get offset according to the green card*/
            get_uj_offset(imm20,&offset);
            /*if offset> 0,tranverse below*/
            if(offset>0){
                int32_t j=0;
                int32_t temp = offset;
                for(j=0;j<(temp/4);j++){
                    /*if compressible, offset minus 2*/
                    if(cmd_info[i+j].state == COMPRESSIBLE)
                        offset-=2;
                }
            }
            /*if offset> 0,tranverse above*/
            else if (offset<0){
                int32_t j=0;
                int32_t temp = offset;
                for(j=-1;-j<=-(temp/4);j--){
                    /*if compressible, offset plus 2*/
                    if(cmd_info[i+j].state == COMPRESSIBLE)
                        offset+=2;
                }
            }
            /*get the new cmd*/
            if(cmd_info[i].state == COMPRESSIBLE){
                /*if compressible*/
                n_cmd|=((offset&0x20)>>3);/*offset 5*/
                n_cmd|=((offset&0xe)<<2); /*offset 3:1*/
                n_cmd|=((offset&0x80)>>1); /*offset 7*/
                /*get new cmd*/
                n_cmd|=((offset&0x40)<<1);/*offset 6*/
                n_cmd|=((offset&0x400)>>2);/*offset 10*/
                n_cmd|=((offset&0x300)<<1);/*offset 9:8*/
                /*get new cmd*/
                n_cmd|=((offset&0x10)<<7);/*offset 4*/
                n_cmd|=((offset&0x800)<<1); /*offset 11*/
                /*c.jal*/
                if(rd)
                    n_cmd|=0x2000;
                /*c.j*/
                else
                    n_cmd|=0xa000;
                /*if compressible, set new cmd*/
                cmd_info[i].cmd = n_cmd;
            }
            else{
                /*set n_offset*/
                uint32_t n_offset=0;
                /*offset 19:12*/
                n_offset |=((offset&0xff000));
                /*offset 11*/
                n_offset |=((offset&800)<<9);
                /*offset 10:1*/
                n_offset |=((offset&0x7fe)<<20);
                /*offset 20*/
                n_offset |=((offset&80000)<<11);
                /*clear old offset*/
                cmd_info[i].cmd&=0xfff;
                /*set new offset*/
                cmd_info[i].cmd|=n_offset;
            }
        }
        /*c.beqz and c.bnez*/
        else if(cmd_info[i].format == SB){
            /*set rd to the required one*/
            uint32_t rd = ((cmd>>15)&REGISTER)-8;
            /*set funct3*/
            uint32_t funct3=(cmd>>12)&FUNCT3;
            /*save the 12 bit of offset and the signed bits*/
            int32_t offset =(((int32_t)cmd>>19)&0xfffff000);
            /*offset 4:1*/
            offset |=(((cmd>>8)&0xf)<<1);
            /*offset 10:5*/
            offset |=(((cmd>>25)&0x3f)<<5);
            /*offset 11*/
            offset |=(((cmd>>7)&0x1)<<11);
            /*if offset> 0,tranverse below*/
            if(offset>0){
                int32_t j=0;
                /*keep the condition unchanged*/
                int32_t temp = offset;
                for(j=0;j<(temp/4);j++){
                    /*if compressible, offset minus 2*/
                    if(cmd_info[i+j].state == COMPRESSIBLE)
                        offset-=2;
                }
            }
            /*if offset> 0,tranverse above*/
            else if (offset<0){
                int32_t j=0;
                /*keep the condition unchanged*/
                int32_t temp = offset;
                for(j=-1;-j<=-(temp/4);j--){
                    /*if compressible, offset plus 2*/
                    if(cmd_info[i+j].state == COMPRESSIBLE)
                        offset+=2;
                }
            }
            if(cmd_info[i].state == COMPRESSIBLE){
                /*get new cmd*/
                n_cmd |=((offset&0x20)>>3);
                n_cmd |=((offset&0x6)<<2);
                n_cmd |=((offset&0xc0)>>1);
                /*get new cmd*/
                n_cmd |=(rd<<7);
                n_cmd |=((offset&0x18)<<7);
                n_cmd |=((offset&0x100)<<4);
                /*c.beqz*/
                if(!funct3)
                    n_cmd |=0xc000;
                /*c.bnez*/
                else
                    n_cmd |=0xe000;
            /*if compressible, set the new cmd*/
                cmd_info[i].cmd = n_cmd;
            }
            else{
                uint32_t n_imm5=0;
                uint32_t n_imm7=0;
                /*clear origin imm*/
                cmd_info[i].cmd&=(~0xf80);
                cmd_info[i].cmd&=(~0xfe000000);
                /*offset 11*/
                n_imm5|=((offset>>11)&0x1);
                /*offset 4:1*/
                n_imm5|=(((offset>>1)&0xf)<<1);
                /*offset 10:5*/
                n_imm7|=((offset>>5)&0x3f);
                /*offset 12*/
                n_imm7|=(((offset>>12)&0x1)<<6);
                /*get new cmd*/
                cmd_info[i].cmd|=(n_imm5<<7);
                cmd_info[i].cmd|=(n_imm7<<25);
            }
        }
    }
}

void format_compressible_check(cmd_info_t * cmd_info,uint32_t * cmd_data){
    uint32_t i =0;
    /*transvse the whole cmd to check every one*/
    for(i=0;i<cmd_num;i++){
        cmd_info[i].cmd = cmd_data[i];
        switch(cmd_data[i] & OPCODE){
            /*R type, check whether compressible*/
            case R:
                cmd_info[i].format = R;
                R_check(&cmd_info[i]);
                break;
            /*I type, check whether compressible*/
            case I:
                cmd_info[i].format = I;
                I_check(&cmd_info[i]);
                break;
            /*LI type, check whether compressible*/
            case LI:
                cmd_info[i].format = LI;
                LI_check(&cmd_info[i]);
                break;
            /*U type, check whether compressible*/
            case U:
                cmd_info[i].format = U;
                U_check(&cmd_info[i]);
                break;
            /*S type, check whether compressible*/
            case S:
                cmd_info[i].format = S;
                S_check(&cmd_info[i]);
                break;
            /*SB type, check whether compressible*/
            case SB:
                cmd_info[i].format = SB;
                SB_check(&cmd_info[i]);
                break;
            /*UJ type, check whether compressible*/
            case UJ:
                cmd_info[i].format = UJ;
                UJ_check(&cmd_info[i]);
                break;
            /*JI type, check whether compressible*/
            case JI:
                cmd_info[i].format = JI;
                JI_check(&cmd_info[i]);
                break;
            /*otherwise*/
            default:
                cmd_info[i].format = UNVAILD;
                cmd_info[i].state = INCOMPRESSIBLE;
        }
    }
}

cmd_info_t * compress(uint32_t * cmd_data){
    /*get a cmd_info array to store all the information of the cmd*/
    cmd_info_t * cmd_info;
    cmd_info = (cmd_info_t *)calloc(sizeof(cmd_info_t),cmd_num);
    /*check their type and whether compressible*/
    format_compressible_check(cmd_info,cmd_data);
    /*unuse the cmd_data*/
    free(cmd_data);
    /*handle compressible*/
    handle_compressible(cmd_info);
    /*handle uj and sb*/
    handle_unsure(cmd_info);
    return cmd_info;
}

/* Your code here... */

