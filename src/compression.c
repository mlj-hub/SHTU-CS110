#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "compression.h"
#include "utils.h"

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
        uint32_t Old_rs2 = (cmd_info->cmd>20)&REGISTER;
        if (Old_rd == Old_rs1)
        {
            /*the cmd is [c.add]*/
            CR_funct4 = 0x9;
            CR_rdORrs1 = Old_rd;
            CR_rs2 = Old_rs2;
            CR_op = 0x2;

        }
        else
        {
            CR_funct4 = 0x8;
            CR_rdORrs1 = Old_rd;
            CR_rs2 = Old_rs2;
            CR_op = 0x2;
            /* the cmd is [c.mv]*/
        }
    }
    /* if the RISCV32 command is I,then in CR_compress the command is [jalr]*/
    else if(cmd_info->format == I)
    {
        uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
        if (Old_rd == 0x0)
        {
            /* the cmd is [c.jr] */
            CR_funct4 = 0x8;
            CR_rdORrs1 = Old_rd;
            CR_rs2 = 0;
            CR_op = 0x2;
        }
        else
        {
            CR_funct4 = 0x9;
            CR_rdORrs1 = Old_rd;
            CR_rs2 = 0;
            CR_op = 0x2;
            /* the cmd is [c.jalr] */
        }
    }
    /*change the 32-bit code in the cmd_info*/
    cmd_info->cmd = CR_op + (CR_rs2<<2) + (CR_rdORrs1<<7)+(CR_funct4<<12);
    return;

}

void CS_compress(cmd_info_t * cmd_info)
{
    /*if the RISCV32 command is S,then in CS_compress the command is [sw] */
    if (cmd_info->format == S)
    {
        /*parameters of the RISCV-16  Command CS_T1*/
        uint32_t CS_T1_funct3 = 0;
        uint32_t CS_T1_IMM3 = 0;
        uint32_t CS_T1_RS1 = 0;
        uint32_t CS_T1_IMM2 = 0;
        uint32_t CS_T1_RS2 = 0;
        uint32_t CS_T1_OP = 0;
        /*parameters of the RISCV-32 Command CS_T2*/
        uint32_t Old_imm12 = 0;
        uint32_t Old_rs1 = (cmd_info->cmd>>15)&REGISTER;
        uint32_t Old_rs2 = (cmd_info->cmd>>20)&REGISTER;
        /*get imm from two part of the cmd*/
        Old_imm12 |= (cmd_info->cmd>>7)&IMM5;
        Old_imm12 |= (cmd_info->cmd>>25)&IMM7;
        /* the  cmd is [c.sub]*/
        CS_T1_funct3 = 6;
        CS_T1_RS1 = Old_rs1-8;
        CS_T1_RS2 = Old_rs2-8;
        /*get the offset[5:3] of given imm*/
        CS_T1_IMM3 = (Old_imm12 >> 3) & 0x7;
        /*get the offset[2|6] of given imm*/
        CS_T1_IMM2 += ((Old_imm12 >> 2) & 0x1) <<1 ;
        CS_T1_IMM2 += (Old_imm12 >> 6) & 0x1;
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
        uint32_t Old_funct7 = (cmd_info->cmd>>25)&FUNCT7;

        uint32_t CS_T2_funct6 = 0x23;
        uint32_t CS_T2_RDorRS1 = Old_rs1 - 8;
        uint32_t CS_T2_funct2 = 0;
        uint32_t CS_T2_RS2 = Old_rs2 - 8;
        uint32_t CS_T2_OP = 1;

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
        return;
    }



}
/*check add,and,or,xor,sub*/
void CL_compress(cmd_info_t *cmd_info)
{
    /* the parameter of I*/
    uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
    uint32_t Old_rs1 = (cmd_info->cmd>>15)&REGISTER;
    int32_t  Old_imm12 = (cmd_info->cmd>>20)&IMM12;
    /* the parameter of CL*/
    uint32_t CL_funct3 = 2;
    uint32_t CL_IMM3 = (Old_imm12>>3) & 0x7;
    uint32_t CL_RS1 = Old_rs1 - 8;
    uint32_t CL_IMM2 = 0;
    uint32_t CL_RD= Old_rd - 8;
    uint32_t CL_OP = 0;
    CL_IMM2 += ((Old_imm12 >> 2) & 0x1) <<1 ;
    CL_IMM2 += (Old_imm12 >> 6) & 0x1;
    cmd_info->cmd = CL_OP + (CL_RD << 2) + (CL_IMM2<<5) +\
                    (CL_RS1 << 7) + (CL_IMM3 << 10) + (CL_funct3 << 13);
    return ;
}

void CI_compress(cmd_info_t *cmd_info)
{
    /*CI parameters*/
    uint32_t CI_funct3 = 0;
    uint32_t CI_IMM1 = 0;
    uint32_t CI_RS1orRD = 0 ;
    uint32_t CI_IMM5 = 0;
    uint32_t CI_OP = 0;
    /*the input command is lui,the old type is U*/
    if(cmd_info->format == U)
    {
        /*change the value of CI ,the commmand is [c.LUI]*/
        uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
        int32_t Old_imm20 = (cmd_info->cmd>>12)&IMM20;
        CI_funct3 = 3;
        CI_IMM1 = (Old_imm20 >> 17) & 0x1;
        CI_RS1orRD = Old_rd;
        CI_IMM5 = (Old_imm20 >> 12) & IMM5;
        CI_OP = 1;
    }
    /* the old type is I,the input is addi or slli*/
    else if (cmd_info->format == I)
    {
        uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
        uint32_t Old_rs1 = (cmd_info->cmd>>15)&REGISTER;
        uint32_t Old_funct3 = (cmd_info->cmd>>12)&FUNCT3;
        int32_t  Old_imm12 = (cmd_info->cmd>>20)&IMM12;
        if (Old_funct3 == 0)
        {
            /* the old cmd is addi */
            CI_IMM1 = (Old_imm12>>5) & 0x1;
            CI_RS1orRD = Old_rd;
            CI_IMM5 =  Old_imm12 & IMM5;
            CI_OP = 1;
            if (Old_rs1 == 0)
            {
                /* the commmand is [c.li]*/
                CI_funct3 = 2;
            }
            else
            {
                /* the command is [c.slli]*/
                CI_funct3 = 0;
            }
        }
        else
        {
            /* the old cmd is slli*/
            CI_funct3 = 0;
            /*shamt[5] should be zero, we directly assigned 0 to imm1*/
            CI_IMM1 = 0;
            CI_RS1orRD = Old_rd;
            CI_IMM5 =  Old_imm12 & IMM5;
            CI_OP = 2;

        }
    }
    /*changt the value of the command*/
    cmd_info->cmd = CI_OP + (CI_IMM5 <<2) + (CI_RS1orRD << 7) +\
                    (CI_IMM1<<12) + (CI_funct3 << 13);
    return;

}
void CB_T2_Compress(cmd_info_t *cmd_info)
{
        /*CI parameters*/
    uint32_t CB_T2_funct3 = 0x4;
    uint32_t CB_T2_IMM1 = 0;
    uint32_t CB_T2_funct2 = 0;
    uint32_t CB_T2_RS1orRD = 0 ;
    uint32_t CB_T2_IMM5 = 0;
    uint32_t CB_T2_OP = 1;
    /* given parameters,type I*/
    uint32_t Old_rd = (cmd_info->cmd>>7)&REGISTER;
    uint32_t Old_funct3 = (cmd_info->cmd>>12)&FUNCT3;
    uint32_t Old_imm12 = (cmd_info->cmd>>20)&IMM12;
    uint32_t Old_shamt = (cmd_info->cmd>>20)&REGISTER;
    uint32_t Old_funct7 = (cmd_info->cmd>>25)&FUNCT7;
    if (Old_funct3 == 0x5)
    {
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
    cmd_info->cmd = CB_T2_OP + (CB_T2_IMM5<<2) + (CB_T2_RS1orRD << 7) +\
                    +(CB_T2_funct2<< 10) + (CB_T2_IMM1<<12)+(CB_T2_funct3 << 13);

    return;
}

void R_check(cmd_info_t * cmd_info){
    uint32_t cmd = cmd_info->cmd;
    /*get parameters according to the format*/
    uint32_t rd = (cmd>>7)&REGISTER;
    uint32_t funct3 = (cmd>>12)&FUNCT3;
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
                    cmd_info->state = COMPRESSIBLE;
                    cmd_info->c_format = CR;
                    return ;
                }
            }
            /*c.sub*/
            else{
                if(rd==rs1 && rd>=8 && rd<=15){
                    cmd_info->state = COMPRESSIBLE;
                    cmd_info->c_format = CS_T2;
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
                return ;
            }
            cmd_info->state = INCOMPRESSIBLE;
            break;
        default:
            cmd_info->state = INCOMPRESSIBLE;
    }
}

void I_check(cmd_info_t * cmd_info){
    uint32_t cmd = cmd_info -> cmd;
    uint32_t rd = (cmd>>7)&REGISTER;
    uint32_t rs1 = (cmd>>15)&REGISTER;
    uint32_t funct3 = (cmd>>12)&FUNCT3;
    int32_t  imm12 = ((int32_t)cmd>>20)&SIGN_ALL;
    switch(funct3){
        /*addi*/
        case 0x0:
            /*c.li rd!=0,rs1=0,-32<=imm12<=31*/
            if(rd && !rs1 && imm12<=31 && imm12>=-32){
                cmd_info->state = COMPRESSIBLE;
                cmd_info->c_format = CI;
                return;
            }
            /*c.addi rd=rs1!=0, imm!=0*/
            else if(rd==rs1 && rd && imm12 && imm12<=31 && imm12>=-32){
                cmd_info->state = COMPRESSIBLE;
                cmd_info->c_format = CI;
                return;
            }
            /*otherwise it is incompressible*/
            cmd_info ->state = INCOMPRESSIBLE;
            break;
        /*slli*/
        case 0x1:
            if(((cmd>>20)&0x20) ==0 && rd==rs1 && !rd){
                cmd_info -> c_format = CI;
                cmd_info -> state = COMPRESSIBLE;
                return;
            }
            cmd_info->state = INCOMPRESSIBLE;
            break;
        /*srli or srai*/
        case 0x5:
            if(rd == rs1 && ((cmd>>20)&0x20) ==0 && rd>=8 && rd<=15){
                cmd_info -> c_format = CB_T2;
                cmd_info -> state = COMPRESSIBLE;
                return;
            }
            cmd_info->state = INCOMPRESSIBLE;
            break;
        /*andi*/
        case 0x7:
            /*c.andi rd=rs1, -32<=imm12<=31*/
            if(rd==rs1 && -32<=imm12 && imm12<=31 && rd>=8 && rd<=15){
                cmd_info->state = COMPRESSIBLE;
                cmd_info->c_format = CB_T2;
                return;
            }
            cmd_info -> state = INCOMPRESSIBLE;
            break;
        default:
            cmd_info -> state = INCOMPRESSIBLE;
    }
}

/*check lw*/
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
    else if ((imm12>>2)>=64)
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
    else if ((imm12>>2)>=64)
        cmd_info->state = INCOMPRESSIBLE;
    /*otherwise compressible*/
    else{
        cmd_info->state = COMPRESSIBLE;
        cmd_info ->c_format = CS_T1;
    }
}
/*bne,beq*/
void SB_check(cmd_info_t * cmd_info){
    uint32_t cmd = cmd_info->cmd;
    uint32_t funct3 = (cmd>>12)&REGISTER;
    uint32_t rs1 = (cmd>>15)&REGISTER;
    uint32_t rs2 = (cmd>>20)&REGISTER;
    if(funct3 >=2)
        cmd_info->state = INCOMPRESSIBLE;
    else if(rs1>=16 || rs1<=7)
        cmd_info->state = INCOMPRESSIBLE;
    else if(rs2!=0)
        cmd_info->state = INCOMPRESSIBLE;
    else{
        cmd_info->state = COMPRESSIBLE;
        cmd_info->c_format = CB_T1;
    }
}

/*check jalr*/
void JI_check(cmd_info_t * cmd_info){
    uint32_t cmd = cmd_info -> cmd;
    /*get parameters according to the format*/
    uint32_t funct3 = (cmd>>12)&FUNCT3;
    uint32_t rd = (cmd>>7)&REGISTER;
    uint32_t rs1 = (cmd>>15)&REGISTER;
    uint32_t  imm12 = (cmd>>20)&IMM12;
    /*check funct3 to determine whether the command is jalr*/
    if(funct3 != 0x0)
        cmd_info->state = INCOMPRESSIBLE;
    /*incompressible if imm is non-zero or rs1 is x0 or rd is not x0 or x1 */
    else if(imm12 || !rs1 || rd>=2)
        cmd_info->state =INCOMPRESSIBLE;
    else{
        cmd_info->state= COMPRESSIBLE;
        cmd_info->c_format = CR;
    }
}

void U_check(cmd_info_t * cmd_info){
    /*since there two U_type commands and they are different in the opcode*/
    uint32_t cmd = cmd_info -> cmd;
    uint32_t rd = (cmd >>7)&REGISTER;
    int32_t imm20 = (cmd>>12)&IMM20;
    if(rd !=0 && rd!=2 && !imm20 && imm20<=31 && imm20>=-32){
        cmd_info->c_format = CI;
        cmd_info->state = COMPRESSIBLE;
        return;
    }
    cmd_info -> state =INCOMPRESSIBLE;
}
/*jar*/
void UJ_check(cmd_info_t * cmd_info){
    uint32_t cmd = cmd_info->cmd;
    uint32_t rd = (cmd>>7)&REGISTER;
    if(rd>=2)
        cmd_info->state = INCOMPRESSIBLE;
    else{
        cmd_info->state = COMPRESSIBLE;
        cmd_info ->c_format = CJ;
    }
}

void handle_compressible(cmd_info_t * cmd_info){
    uint32_t i=0;
    for(i=0;i<cmd_num;i++){
        if(cmd_info[i].state != COMPRESSIBLE)
            continue;
        switch(cmd_info[i].c_format){
            case CR:
                CR_compress(&cmd_info[i]);
                break;
            case CI:
                CI_compress(&cmd_info[i]);
                break;
            case CL:
                CL_compress(&cmd_info[i]);
                break;
            case CS_T1: case CS_T2:
                CS_compress(&cmd_info[i]);
                break;
            case CB_T2:
                CB_T2_Compress(&cmd_info[i]);
                break;
            default:
                continue;
        }
    }
}

void handle_unsure(cmd_info_t * cmd_info){
    uint32_t i=0;
    uint32_t n_cmd = 0x1;
    for(i=0;i<cmd_num;i++){
        uint32_t cmd = cmd_info[i].cmd;
        /*c.j and c.jar*/
        if(cmd_info[i].c_format == CJ){
            uint32_t rd = (cmd>>7)&REGISTER;
            uint32_t imm20 = (cmd>>12)&IMM20;
            int32_t offset = 0;
            /*get offset according to the green card*/
            offset |=((imm20&0xff)<<12);
            offset |=((imm20&0x100)<<3);
            offset |=((imm20&0x7fe00)>>9);
            offset |=(imm20&0x80000);
            if(offset>0){
                int32_t j=0;
                for(j=0;j<(offset/4);j++){
                    if(cmd_info[i+j].state == COMPRESSIBLE)
                        offset-=2;
                }
            }
            else if (offset<0){
                int32_t j=0;
                for(j=-1;-j<-(offset/4);j--){
                    if(cmd_info[i+j].state == COMPRESSIBLE)
                        offset+=2;
                }
            }
            /*get the new cmd*/
            n_cmd|=((offset&0x20)>>3);
            n_cmd|=((offset&0xe)<<2);
            n_cmd|=((offset&0x80)>>1);
            n_cmd|=((offset&0x40)<<1);
            n_cmd|=((offset&0x400)>>2);
            n_cmd|=((offset&0x300)<<1);
            n_cmd|=((offset&0x10)<<6);
            n_cmd|=((offset&0x800)<<1);
            /*c.jal*/
            if(rd)
                n_cmd|=0x2000;
            /*c.j*/
            else
                n_cmd|=0xa000;
            cmd_info[i].cmd = n_cmd;
        }
        /*c.beqz and c.bnez*/
        else if(cmd_info[i].c_format == CB_T1){
            uint32_t rd = ((cmd>>15)&REGISTER)-8;
            int32_t offset =0;
            uint32_t funct3=(cmd>>12)&FUNCT3;
            offset |=((cmd&0x80)<<4);
            offset |=((cmd&0xf00)>>7);
            offset |=((cmd&0x7e000000)>>20);
            offset |=((cmd&0x80000000)>>19);
            if(offset>0){
                int32_t j=0;
                for(j=0;j<(offset/4);j++){
                    if(cmd_info[i+j].state == COMPRESSIBLE)
                        offset-=2;
                }
            }
            else if (offset<0){
                int32_t j=0;
                for(j=-1;-j<-(offset/4);j--){
                    if(cmd_info[i+j].state == COMPRESSIBLE)
                        offset+=2;
                }
            }
            n_cmd |=((offset&0x20)>>3);
            n_cmd |=((offset&0x6)<<2);
            n_cmd |=((offset&0xc0)>>1);
            n_cmd |=(rd<<7);
            n_cmd |=((offset&0x18)<<7);
            n_cmd |=((offset&0x100)<<4);
            /*c.beqz*/
            if(!funct3)
                n_cmd |=0xc000;
            /*c.bnez*/
            else
                n_cmd |=0xe000;
            cmd_info[i].cmd = n_cmd;
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
    return ;
}

cmd_info_t * compress(uint32_t * cmd_data){
    /*get a cmd_info array to store all the information of the cmd*/
    cmd_info_t * cmd_info;
    cmd_info = (cmd_info_t *)calloc(sizeof(cmd_info_t),cmd_num);
    /*check their type and whether compressible*/
    format_compressible_check(cmd_info,cmd_data);
    /*unuse the cmd_data*/
    free(cmd_data);
    handle_compressible(cmd_info);
    handle_unsure(cmd_info);
    return cmd_info;
}

/* Your code here... */

