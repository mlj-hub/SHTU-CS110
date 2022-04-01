#   Project 1.2 RVC instructions to RISC-V instructions in RISC-V

#   You can suppose you will get the number of machine codes you should translate from 
#   "lines_of_rvc_codes" and get the 16 or 32 bits-machine codes from "rvc_codes". 
#   Both of them are stored in file "input.S"

#   When we test your project, we will change the input files content.

# DO NOT REMOVE THIS LINE WHEN SUBMITTING
.import input.S

.text
# ================================================================================================================ #
# codes below are the main function

main:
    # load line_num to a1
    lw a1,lines_of_rvc_codes
    # load 4 to t1 for feature usage
    li t1,4
    # mul len by 4 to get the size of total cmd and store it in a1
    mul a1,a1,t1
    # malloc for data:
    # assign ID 9 to a0
    li a0,9 
    # malloc a1 bytes and store the pointer to a0
    ecall
    # move the data pointer into a2 
    mv a2,a0
    # malloc for flag
    li a0,9
    # malloc a1 bytes and store the pointer to a0
    ecall
    # move the flag pointer into a3
    mv a3,a0
    # load back line_num into a1
    lw a1,lines_of_rvc_codes

    # handle input
    jal input_handler
    # de compress
    jal de_compress
    # handle output
    jal output_handler
    
    # Exit with error code 0
    add	   a1, x0, x0
    addi   a0, x0, 17 
    ecall

# codes above are the main function
# ================================================================================================================ #

# ================================================================================================================ #
# codes below is used to handle input

input_handler:
    # move the stack three words below
    addi sp,sp,-16
    # store the return address
    sw ra,12(sp)
    # store the line num 
    sw a1,8(sp)
    # store the state array pointer
    sw a3,4(sp)
    # store the data array pointer
    sw a2,0(sp)
    # move line_num from a1 into t4
    mv t4,a1
    # load adress of input data to a1
    la a1,rvc_codes 
    loop_read_in:
        # /if already read in all commands, return
        beqz t4,return_input_handler
        # minus t4 by one to indicate that read in one cmd
        addi t4,t4,-1
        # load half word into t0
        lh t0,0(a1)
    # check if the cmd is 16 bits
        # get the rightest two bits of t0 and store it to t1
        andi t1,t0,3
        # load 3 into t2 for feather usage
        li t2,3
        # /if t1 and t2 is equal, it is 32 bits,branch to 32_bits_cmd
        beq t1,t2,32_bits_cmd
        16_bits_cmd:
            # left shift cmd by 16 bits to clear the sign bits
            slli t0,t0,16
            # right shift cmd by 16 bits to make the sign bits all zero
            srli t0,t0,16
            # store cmd into a2
            sw t0,0(a2)
            # store 1 into t3 for flag
            li t3, 1
            # store state: 1 for 16 bits
            sw t3,0(a3)
            # move a2(data pointer) one word forward
            addi a2,a2,4
            # move a3(flag pointer) one word forward
            addi a3,a3,4
            # move a1(input pointer) half word forward
            addi a1,a1,2
            # continue loop
            j loop_read_in
        32_bits_cmd:
            # load word into t0
            lw t0, 0(a1)
            # store cmd into a2(data pointer)
            sw t0, 0(a2)
            # store 1 into t3 for flag
            li t3, 0
            # store state: 1 for 16 bits
            sw t3,0(a3)
            # move a2(data pointer) one word forward
            addi a2,a2,4
            # move a3(state pointer) one word forward
            addi a3,a3,4
            # move a1(input pointer) one word forward
            addi a1,a1,4
            # continue loop
            j loop_read_in
return_input_handler:
    # load return address
    lw ra,12(sp)
    # load back the line num
    lw a1,8(sp)
    # load back the state pointer
    lw a3,4(sp)
    # load back the data pointer
    lw a2,0(sp)
    # reset sp
    addi sp,sp,16
    jr ra
# codes above is used to handle input data
# ================================================================================================================ #

# ================================================================================================================ #
# codes below are main funciton of de compressing

de_compress:
    # get three words in the stack
    addi sp,sp,-16
    # store return address
    sw ra,12(sp)
    # store len num
    sw a1,8(sp)
    # store data array
    sw a2,4(sp)
    # store state array
    sw a3,0(sp)
    # call check_type
    jal check_type
    # load back len num
    lw a1,8(sp)
    # load back data array
    lw a2,4(sp)
    # load back state array
    lw a3,0(sp)
    # call change_offset
    jal change_offset
    # return 
    j return_de_compress

return_de_compress:
    # load back return address
    lw ra,12(sp)
    # load back len num
    lw a1,8(sp)
    # load back data array
    lw a2,4(sp)
    # load back state array
    lw a3,0(sp)
    # reset stack
    addi sp,sp,16
    jr ra
# codes above are main funciton of de compressing
# ================================================================================================================ #

# ================================================================================================================ #
# codes below are main functions of check_type

check_type:
    # let t0 be i and initialize it to zero
    li t0,0
    # get one word on stack
    addi sp , sp,-4
    # store return address
    sw ra,0(sp)
    loop_check_type:
        # /if i== len num , return
        beq t0,a1,return_check_type
        # get one word on the stack
        addi sp,sp,-4
        # store t0 into stack
        sw t0,0(sp)
        # load cmd state into s2
        lw s2,0(a3)
        # get the address of cmd into a4
        mv a4,a2
        # load cmd into t1
        lw t1,0(a4)
        # cmd is 32bits, check whether if it is jump or branch
        beqz s2,b_j_check_32
        # get opcode into t2
        andi t2,t1,0x03
        # load 1 into t3 for opcode check
        li t3,1
        # check 00
        beq t2,zero,check_00
        # check 01
        beq t2,t3,check_01
        # check 10
        j check_10

check_00:
    # right shift the cmd by 15bits to get the highest bit
    srli t2,t1, 15
    # /if the highest bit is 0, the cmd is clw
    beqz t2,clw
    # /else, the cmd is csw, call csw_decompress
    csw: jal csw_decompress
    # continue
    j continue_check_type
    # call clw decompress
    clw: jal clw_decompress
    # continue 
    j continue_check_type

check_01:
    # right shift the cmd by 13 to get the highest thress bits into t2
    srli t2,t1,13
    # get 0x1 into s1 for type check
    li s1,0x1
    # get 0x2 into s2 for type check
    li s2,0x2
    # get 0x3 into s3 for type check
    li s3,0x3
    # get 0x5 into s5 for type check
    li s5,0x5
    # get 0x6 into s6 for type check
    li s6,0x6
    # get 0x7 into s7 for type check
    li s7,0x7
    # get 0x8 into s8 for type check
    li s8,0x8
    # /if funct3 is 000, the cmd is caddi
    beqz t2,caddi
    # /if funct3 is 001, the cmd is cjal
    beq t2,s1,cjal_cj
    # /if funct3 is 010, the cmd is cli
    beq t2,s2,cli
    # /if funct3 is 011, the cmd is clui
    beq t2,s3,clui
    # /if funct3 is 101, the cmd is cj
    beq t2,s5,cjal_cj
    # /if funct3 is 110, the cmd is cbeqz
    beq t2,s6,cbeqz_cbnez
    # /if funct3 is 111, the cmd is cbnez
    beq t2,s7,cbeqz_cbnez
    # /else, continue checking 
    # right shift the cmd by 10 bits to get the highest 6 bits into t2
    srli t2,t1,10
    # and t2 by 0x3 to get the lowese two bits 
    andi t2,t2,0x3
    # /if funct2 is 00, the cmd is csrli
    beqz t2,csrli
    # /if funct2 is 01, the cmd is csrai
    beq t2,s1,csrai
    # /if funct2 is 10, the cmd is candi
    beq t2,s2,candi
    # /else, continue checking
    # get the 6th to 5th bits into t2
    andi t2,t1,0x60
    # right shift t2 by 5 bits to get funct2
    srli t2,t2,5
    # /if funct2 is 00, the cmd is csub
    beqz t2, csub
    # /if funct2 is 01, the cmd is cxor
    beq t2,s1, cxor
    # /if funct2 is 10, the cmd is cor
    beq t2,s2, cor
    # /if funct2 is 11, the cmd is cand
    beq t2,s3,cand
    # /else, finish 01 check, just continue
    j continue_check_type
    # below are codes for calling decompress functions
    # call caddi decompress
    caddi: jal caddi_decompress
    # after decompress, continue
    j continue_check_type

    # handle cjal
    cjal_cj: 
    # set 5(16bits jump) to t2
    li t2,5
    # store t2(5) into state array
    sw t2,0(a3)
    # call cjal decompress
    jal cjal_cj_decompress
    # after decompress, continue
    j continue_check_type

    # call cli decompress
    cli: jal cli_decompress
    # after decompress, continue
    j continue_check_type

    # call clui decompress
    clui: jal clui_decompress
    # after decompress, continue
    j continue_check_type

    # handle cbeqz
    cbeqz_cbnez: 
    # set 3 (16bits branch) to t2
    li t2,3
    # store t2(3) into state array
    sw t2,0(a3)
    # call cbeqz decompress
    jal cbeqz_cbnez_decompress
    # after decompress, continue
    j continue_check_type

    # call cand decompress
    cand: jal cand_decompress
    # after decompress, continue
    j continue_check_type

    # call cor decompress
    cor: jal cor_decompress
    # after decompress, continue
    j continue_check_type

    # call xor decompress
    cxor: jal cxor_decompress
    # after decompress, continue
    j continue_check_type

    # call csub decompress
    csub: jal csub_decompress
    # after decompress, continue
    j continue_check_type

    # call csrli decompress
    csrli: jal csrli_decompress
    # after decompress, continue
    j continue_check_type

    # call csrai decompress
    csrai: jal csrai_decompress
    # after decompress, continue
    j continue_check_type

    # call candi decompress
    candi: jal candi_decompress
    # after decompress, continue
    j continue_check_type

check_10:
    # right shift the cmd by 15 to get the highest bit into t2
    srli t2,t1,15
    # /if the highest bit of the cmd is 0, the cmd is cslli
    beqz t2,cslli
    # get rs2 of the cmd into t2
    andi t2,t1,0x7c
    # /if rs2 is 0, the cmd is cjalr or cjr
    beqz t2,cjr_or_cjalr
    # /else , the cmd is cmv or cadd
    j cadd_or_cmv
    # check cjr and cjalr
    cjr_or_cjalr:
        # right shift the cmd by 12 to get the highest byte
        srli t2,t1,12
        # and t2 with 0x1 to get the lowest bit
        andi t2,t1,0x1
        # /if 13rd bit is 0, the cmd is cjr
        beqz t2,cjr
        # /else, the cmd is cjalr, call cjalr decompress
        cjalr: jal cjalr_decompress
        # continue
        j continue_check_type
        # call cjr decompress
        cjr: jal cjr_decompress
        # continue
        j continue_check_type
    # check cadd and cmv
    cadd_or_cmv:
        # right shift the cmd by 12 to get the highest byte
        srli t2,t1,12
        # and t2 with 0x1 to get the lowest bit
        andi t2,t2,0x1
        # /if 13rd bit is 0, the cmd is cmv
        beqz t2,cmv
        # /else, the cmd is cadd, call cadd decompress
        cadd: jal cadd_decompress
        j continue_check_type
        # call cmv decompress
        cmv: jal cmv_decompress
        j continue_check_type
    # /if the rs2 is 0, the cmd is cjr or cjalr
    cslli: jal cslli_decompress
    j continue_check_type

b_j_check_32:
    # load opcode into t2
    andi t2,t1,0x7f
    # load jump opcode into t3
    li t3,0x6f
    # load branch opcode into t4
    li t4,0x63
    # /if the cmd is jump, set the state array to 2
    beq t2,t3,j_state_set
    # /if the cmd is branch, set the state array to 2
    beq t2,t4,b_state_set
    # /else, leave the state array unchanged, which is 0
    j continue_check_type
    j_state_set:
        # load 4 (32bits jump) into t3
        li t3,4
        # change the state array
        sw t3,0(a3)
        # continue
        j continue_check_type
    b_state_set:
        # load 2 (32bits branch) into t3
        li t3,2
        # change the state array
        sw t3,0(a3)
        # continue
        j continue_check_type


continue_check_type:
    # load back t0
    lw t0,0(sp)
    # reset stack
    addi sp,sp,4
    # i++
    addi t0,t0,1
    # move state array one word forward
    addi a3,a3,4
    # move data array one word forward
    addi a2,a2,4
    # continue loop
    j loop_check_type

return_check_type:
    # load back return address
    lw ra,0(sp)
    # reset stack
    addi sp,sp,4
    # return 
    jr ra



# codes above are main functions of check_type
# ================================================================================================================ #

# ================================================================================================================ #
# codes below is used to decompress the data in a4

cadd_decompress:
    #reset the temp value  to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    # load the cmd from array pointer
    lw t0, 0(a4)
    #right move t0 2 bit and store its value in t2
    srli t2 t0 2
    #get the lowest 5 bit of t2 to have rs2 and store it in t2
    andi t2 t2 0x1f
    #right move t0 7 bit and store its value in t1
    srli t1 t0 7
    #get the lowest 5 bit of t1 to have rs1ORrd and store it in t1
    andi t1 t1 0x1f
    #suppose t3 is the value to put back to a4,add opcode to t3
    addi t3 t3 0x33
    #left move t1 7 bit to get it value to get the rd of the cmd in 32bit
    slli t1 t1 7
    #add t1 to t3 to give rd value to final cmd
    add t3 t3 t1
    #add func3 value to t3,since the func3 of add is zero so we add zero
    addi t3 t3 0
    #left move t1 8 bit more,move 15 bit in total get it value to get the rs1 of the cmd in 32bit
    slli t1 t1 8
    #add t1 to t3 to give rs1 value to final cmd
    add t3 t3 t1
    #left move t2 20 bit,get it value to get the rs2 value of the cmd in 32bit
    slli t2 t2 20
    #add t2 to t3 to give rs2 value to final cmd
    add t3 t3 t2
    #add func7 to t3,since func7 of add is zero,so we add zero to t3
    addi t3 t3 0
    #store the value of t3 to the register
    sw t3 0(a4)
    #return 
    jr ra

cmv_decompress:
    #reset the temp value  to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    # load the cmd from array pointer
    lw t0, 0(a4)
    #right move t0 2 bit and store its value in t2
    srli t2 t0 2
    #get the lowest 5 bit of t2 to have rs2 and store it in t2
    andi t2 t2 0x1f
    #right move t0 7 bit and store its value in t1
    srli t1 t0 7
    #get the lowest 5 bit of t1 to have rs1ORrd and store it in t1
    andi t1 t1 0x1f
    #suppose t3 is the value to put back to a4,add opcode to t3
    addi t3 t3 0x33
    #left move t1 7 bit to get it value to get the rd of the cmd in 32bit
    slli t1 t1 7
    #add t1 to t3 to give rd value to final cmd
    add t3 t3 t1
    #add func3 value to t3,since the func3 of add is zero so we add zero
    addi t3 t3 0
    #add rs1 value to t3,since rs1 is x0 ,so we add zero
    addi t3 t3 0
    #left move t2 20 bit, get it value to get the rs2 value of the cmd in 32bit
    slli t2 t2 20
    #add t2 to t3 to give rs2 value to final cmd
    add t3 t3 t2
    #add func7 to t3,since func7 of add is zero,so we add zero to t3
    addi t3 t3 0
    #store the value of t3 to the register
    sw t3 0(a4)
    #return 
    jr ra
cjr_decompress:
    #reset the temp value  to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load the word from the array pointer
    lw t0 0(a4)
    #add opcode  1100111 0x67 to the decompress cmd t3
    addi t3 t3 0x67
    #right move 7 bit to get the rs1ORrd of cmd
    srli t1 t0 7
    #use andi to get the last 5 bit of t1 to have rs1ORrd of cmd 
    andi t1 t1 0x1f
    #add rd to decompress cmd t3,since rd should be 0,we add 0
    addi t3 t3 0
    #add func3 to decompress cmd t3,since func3 should be 0,we add 0
    addi t3 t3 0
    #left move 15 bit to add the rs1 value to the decompress cmd
    slli t1 t1 15
    #add rs1 value to decompress cmd t3
    add t3 t3 t1 
    #add imm value to t3,since imm is zero,so we add zero
    addi t3 t3 0
    #store the decompress cmd to the address in a4 
    sw t3 0(a4)
    #return 
    jr ra
cjalr_decompress:
    #reset the temp value  to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load the word from the array pointer
    lw t0 0(a4)
    #add opcode  1100111 0x67 to the decompress cmd t3
    addi t3 t3 0x67
    #right move 7 bit to get the rs1ORrd of cmd
    srli t1 t0 7
    #use andi to get the last 5 bit of t1 to have rs1ORrd of cmd 
    andi t1 t1 0x1f
    #give value 1  to t2 since the rd value of jalr is 1
    li t2 1
    #left move t2 7 bit to add the value of rd to the decompress cmd
    slli t2 t2 7
    #add rd to decompress cmd t3
    add t3 t3 t2
    #add func3 to decompress cmd t3,since func3 should be 0,we add 0
    addi t3 t3 0
    #left move 10 bit to add the rs1 value to the decompress cmd
    slli t1 t1 15
    #add rs1 value to decompress cmd t3
    add t3 t3 t1 
    #add imm value to t3,since imm is zero,so we add zero
    addi t3 t3 0
    #store the decompress cmd to the address in a4 
    sw t3 0(a4)
    #return 
    jr ra
cli_decompress:
    #reset the temp value  to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move t0 2 bit to get imm[4:0] from the cli cmd
    srli t1 t0 2
    #get the last 5 bit of t1 to get imm[4:0]
    andi t1 t1 0x1f
    #right move t0 7 bit to get rd from the cli cmd
    srli t2 t0 7
    #get the last 5 bit of t2 to get rd
    andi t2 t2 0x1f
    #right move t0 12 bit to get imm[5] from the cli cmd
    srli t3 t0 12
    #get the last 1 bit of t1 to get imm[5]
    andi t3 t3 1
    #add op code to output cmd,the output cmd is t4,opcode is 0x13
    addi t4 t4 0x13
    #left move t2 7 bit to get the rd of the output cmd
    slli t2 t2 7
    #add rd to the output cmd
    add t4 t4 t2
    #add func3 to output cmd,since the opcode is 0,we add zero
    addi t4 t4 0
    #add rs1 to output cmd ,since the rs1 is x0 ,we add zero
    addi t4 t4 0
    #left move t1 20 bit to get imm[4:0]
    slli t1 t1 20
    #add imm[4:0] to the output cmd
    add t4 t4 t1
    #left move t3 25 bit to get imm[5]
    slli t3 t3 25
    #add imm[5] to the output cmd
    add t4 t4 t3
    #left move 6 bit to sign-extend
    slli t4 t4 6
    #right move 6 bit arithmatically to sign-extend
    srai t4 t4 6 
    #store the output cmd to the final cmd
    sw t4 0(a4)
    #return 
    jr ra
clui_decompress:
    #may have bugs in clui
    #reset the temp value  to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move t0 2 bit to get nzimm[16:12] from the clui cmd
    srli t1 t0 2
    #get the last 5 bit of t1 to get nzimm[16:12]
    andi t1 t1 0x1f
    #right move t0 7 bit to get rd from the cli cmd
    srli t2 t0 7
    #get the last 5 bit of t2 to get rd
    andi t2 t2 0x1f
    #right move t0 12 bit to get imm[5] from the cli cmd
    srli t3 t0 12
    #get the last 1 bit of t1 to get imm[5]
    andi t3 t3 1
    #add op code to the output cmd ,t4 is the output cmd ,the opcode is 0x37
    addi t4 t4 0x37
    #left move t2 7 bit to get the rd of the output cmd
    slli t2 t2 7
    #add rd to the output cmd
    add t4 t4 t2
    #left move t1 12 bit to get nzimm[16:12]
    slli t1 t1 12
    #add nzimm[16:12] to cmd 
    add t4 t4 t1
    #left move t3 17 bit to get nzimm[17]
    slli t3 t3 17
    #add nzimm[17] to cmd
    add t4 t4 t3
    #left move t4 14 vit to make nzimm[17] the top
    slli t4 t4 14
    #right move t4 14 bit arithmatically to make the nizmm[31:18] equal to nzimm[17]
    srai t4 t4 14
    #store the output cmd into the address 
    sw t4 0(a4)
    #return 
    jr ra

caddi_decompress:
    #reset the temp value  to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move t0 2 bit to get nzimm[4:0] from the caddi cmd
    srli t1 t0 2
    #get the last 5 bit of t1 to get nzimm[4:0]
    andi t1 t1 0x1f
    #right move t0 7 bit to get rd from the caddi cmd
    srli t2 t0 7
    #get the last 5 bit of t2 to get rd
    andi t2 t2 0x1f
    #right move t0 12 bit to get nzimm[5] from the caddi cmd
    srli t3 t0 12
    #get the last 1 bit of t1 to get nzimm[5]
    andi t3 t3 1
    #add op code to output cmd,the output cmd is t4,opcode is 0x13
    addi t4 t4 0x13
    #left move t2 7 bit to get the rd of the output cmd
    slli t2 t2 7
    #add rd to the output cmd
    add t4 t4 t2
    #add func3 to output cmd,since the opcode is 0,we add zero
    addi t4 t4 0
    #left move t2 8 bit more to get the rs1 of the output cmd 
    slli t2 t2 8
    #add rs1 to output cmd ,since the rs1 is rd ,we add t2 to t4
    add t4 t4 t2
    #left move t1 20 bit to get nzimm[4:0]
    slli t1 t1 20
    #add nzimm[4:0] to the output cmd
    add t4 t4 t1
    #left move t3 25 bit to get nzimm[5]
    slli t3 t3 25
    #add nzimm[5] to the output cmd
    add t4 t4 t3
    #left move 6 bit to sign-extend
    slli t4 t4 6
    #right move 6 bit arithmatically to sign-extend
    srai t4 t4 6 
    #store the output cmd to the final cmd
    sw t4 0(a4)
    #return 
    jr ra
cslli_decompress:
    #reset the temp value  to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move t0 2 bit to get shamt[4:0] from the cslli cmd
    srli t1 t0 2
    #get the last 5 bit of t1 to get shamt[4:0]
    andi t1 t1 0x1f
    #right move t0 7 bit to get rd from the cslli cmd
    srli t2 t0 7
    #get the last 5 bit of t2 to get rd
    andi t2 t2 0x1f
    #right move t0 12 bit to get shamt[5] from the cslli cmd
    srli t3 t0 12
    #get the last 1 bit of t1 to get shamt[5]
    andi t3 t3 1
    #add op code to output cmd,the output cmd is t4,opcode is 0x13
    addi t4 t4 0x13
    #left move t2 7 bit to get the rd of the output cmd
    slli t2 t2 7
    #add rd to the output cmd
    add t4 t4 t2
    #give t5 the value of func3
    li t5 1
    #left move t5 5 bit to add func3 to the output cmd
    slli t5 t5 12
    #add func3 to output cmd,since the opcode is 0,we add zero
    add t4 t4 t5
    #left move t2 8 bit more to get the rs1 of the output cmd 
    slli t2 t2 8
    #add rs1 to output cmd ,since the rs1 is rd ,we add t2 to t4
    add t4 t4 t2
    #left move t1 20 bit to get shamt[4:0]
    slli t1 t1 20
    #add nzimm[4:0] to the output cmd
    add t4 t4 t1
    #left move t3 25 bit to get shamt[5]
    slli t3 t3 25
    #add nzimm[5] to the output cmd
    add t4 t4 t3
    #store the output cmd to the final cmd
    sw t4 0(a4)
    #return 
    jr ra
clw_decompress:
    #reset the tmp value to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move t0 2 bit to get rd' and store it in t1
    srli t1 t0 2
    #get the last 3 bit of t1 as rd'
    andi t1 t1 0x7
    #add 8 to rd' to change it into rd 
    addi t1 t1 8
    #right move t0  6 bit to get offset[2]
    srli t2 t0 6
    #get the last bit of t2 to get offset[2]
    andi t2 t2 1
    #right move t0 5 bit to get offset[6]
    srli t6 t0 5
    #get the last bit of t6 to get offset[6]
    andi t6 t6 1
    #right move t0 7 bit to get rs1'
    srli t3 t0 7
    #get the last 3 bit of t3 to get rs1'
    andi t3 t3 0x7
    #add 8 to rs1' to change it into rd
    addi t3 t3 8
    #right move t0 5 bit to get offset[5:3]
    srli t4 t0 10
    #get the last 3 bit of t6 to get offset[5:3]
    andi t4 t4 0x7

    #add opcode to the output cmd,the opcode of lw is 0x03
    addi t5 t5 0x3
    #left move rd 7 bit to get rd of output cmd
    slli t1 t1 7
    #add rd to the output cmd
    add t5 t5 t1
    #give func3 value to the t0
    li t0 0x2
    #left move func3 12 bit to have func3 of the output cmd
    slli t0 t0 12
    #add func3 value to output cmd 
    add t5 t5 t0 
    #left move rs1 15 bit to have rs1 of the output cmd
    slli t3 t3 15
    #add rs1 to output cmd 
    add t5 t5 t3
    #left move offset[2] 22 bit to get offset[2] in output cmd 
    slli t6 t6 22
    #add offset[2] to output cmd
    add t5 t5 t6
    #left move offset[5:3] 23 bit to get offset[2] in output cmd 
    slli t4 t4 23
    #add offset[5:3] to output cmd
    add t5 t5 t4
    #left move offset[6] 25 bit to get offset[2] in output cmd 
    slli t6 t6 26
    #add offset[6] to output cmd
    add t5 t5 t6
    #store the value to the address of the cmd
    sw t5 0(a4)
    #return 
    jr ra
csw_decompress:
     #reset the tmp value to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move t0 2 bit to get rs2' and store it in t1
    srli t1 t0 2
    #get the last 3 bit of t1 as rs2'
    andi t1 t1 0x7
    #add 8 to rd' to change it into rs2 
    addi t1 t1 8
    #right move t0  6 bit to get offset[2]
    srli t2 t0 6
    #get the last bit of t2 to get offset[2]
    andi t2 t2 1
    #right move t0 5 bit to get offset[6]
    srli t6 t0 5
    #get the last bit of t6 to get offset[6]
    andi t6 t6 1
    #right move t0 7 bit to get rs1'
    srli t3 t0 7
    #get the last 3 bit of t3 to get rs1'
    andi t3 t3 0x7
    #add 8 to rs1' to change it into rs1
    addi t3 t3 8
    #right move t0 5 bit to get offset[5:3]
    srli t4 t0 10
    #get the last 3 bit of t6 to get offset[5:3]
    andi t4 t4 0x7
    #add opcode to the output cmd,since the output cmd is 
    addi t5 t5 0x23
    #give the value of offset[2] to t0
    mv t0 t2
    #left move t4 1 bit to have offset[5:3]
    slli t4 t4 1
    #add offset[5:3] to t0
    add t0 t0 t4
    #left move t6 4 bit to have offset[6]
    slli t6 t6 4
    #add offset[6] to t0
    add t0 t0 t6
    #put offset[6:5] into t4
    srli t4 t0 3
    #put offset[4:2] into t6
    andi t6 t0 3
    #left move t6 9 bit to add it to the output cmd
    slli t6 t6 9
    #add offser[4:2] to imm[4:0]
    add t5 t5 t6
    #give the func3 value to t0 
    li t0 0x2
    #left move li 12 bit to have the func3 of output
    slli t0 t0 12
    #add func3 to output cmd ,the output cmd is 
    add t5 t5 t0
    #left move t4 25 bit to add offset[6:5] to output cmd
    slli t4 t4 25
    #add offset[6:5] to output cmd
    add t5 t5 t4
    #left move t1 20 bit to have rs2 of the output cmd
    slli t1 t1 20
    #add rs2 to the output cmd
    add t5 t5 t1
    #left move t3 20 bit to have rs1 of the output cmd
    slli t3 t3 15
    #add rs1 ti the output cmd
    add t5 t5 t3 
    #store the output cmd into address a4
    sw t5 0(a4)
    #return 
    jr ra
cand_decompress:
    #reset the tmp value to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move t0 2 bit to get rs2'
    srli t2 t0 2
    #get the last 3 bit of t2 to have rs2'
    andi t2 t2 0x7
    #add 8 to t2  to change rs2' to rs2
    addi t2 t2 8
    #right move t0 7 bit to get rs1'
    srli t1 t0 7
    # get the last 3 bit of t2 to have rs1'
    andi t1 t1 0x7
    #add 8 to t1  to change rs1' to rs1
    addi t1 t1 8
    #add opcode to output cmd t3 ,opcode is 0x33
    addi t3 t3 0x33
    #left move t1 7 bit to get rd value
    slli t1 t1 7
    #add rd to output cmd
    add t3 t3 t1
    #give the func3 value to t0
    li t0 7
    #left move t0 12 bit to have the func3 of output
    slli t0 t0 12
    #add func3 to outputcmd t3
    add t3 t3 t0
    #left move t1 8 bit more to get rs1 value 
    slli t1 t1 8
    #add rs1 value to output cmd t3
    add t3 t3 t1
    #left move t2 20 bit to get rs2
    slli t2 t2 20
    #add  rs2 to the output cmd
    add t3 t3 t2
    #add func7 to output cmd ,since func 7 is zero,we add zero
    addi t3 t3 0
    #store the decompressed cmd to the address of a4
    sw t3 0(a4)
    #return 
    jr ra
cor_decompress:
    #reset the tmp value to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move t0 2 bit to get rs2'
    srli t2 t0 2
    #get the last 3 bit of t2 to have rs2'
    andi t2 t2 0x7
    #add 8 to t2  to change rs2' to rs2
    addi t2 t2 8
    #right move t0 7 bit to get rs1'
    srli t1 t0 7
    #get the last 3 bit of t2 to have rs1'
    andi t1 t1 0x7
    #add 8 to t1  to change rs1' to rs1
    addi t1 t1 8
    #add opcode to output cmd t3 ,opcode is 0x33
    addi t3 t3 0x33
    #left move t1 7 bit to get rd value
    slli t1 t1 7
    #add rd to output cmd
    add t3 t3 t1
    #give the func3 value to t0
    li t0 6
    #left move t0 12 bit to have the func3 of output
    slli t0 t0 12
    #add func3 to outputcmd t3,func3 is 0x6000
    add t3 t3 t0
    #left move t1 8 bit more to get rs1 value 
    slli t1 t1 8
    #add rs1 value to output cmd t3
    add t3 t3 t1
    #left move t2 20 bit to get rs2
    slli t2 t2 20
    #add  rs2 to the output cmd
    add t3 t3 t2
    #add func7 to output cmd ,since func 7 is zero,we add zero
    addi t3 t3 0
    #store the decompressed cmd to the address of a4
    sw t3 0(a4)
    #return
    jr ra
cxor_decompress:
    #reset the tmp value to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move t0 2 bit to get rs2'
    srli t2 t0 2
    #get the last 3 bit of t2 to have rs2'
    andi t2 t2 0x7
    #add 8 to t2  to change rs2' to rs2
    addi t2 t2 8
    #right move t0 7 bit to get rs1'
    srli t1 t0 7
    #get the last 3 bit of t2 to have rs1'
    andi t1 t1 0x7
    #add 8 to t1  to change rs1' to rs1
    addi t1 t1 8
    #add opcode to output cmd t3 ,opcode is 0x33
    addi t3 t3 0x33
    #left move t1 7 bit to get rd value
    slli t1 t1 7
    #add rd to output cmd
    add t3 t3 t1
    #give the func3 value to t0
    li t0 4
    #left move t0 12 bit to have the func3 of output
    slli t0 t0 12
    #add func3 to outputcmd t3,func3 is 0x4000
    add t3 t3 t0
    #left move t1 8 bit more to get rs1 value 
    slli t1 t1 8
    #add rs1 value to output cmd t3
    add t3 t3 t1
    #left move t2 20 bit to get rs2
    slli t2 t2 20
    #add  rs2 to the output cmd
    add t3 t3 t2
    #add func7 to output cmd ,since func 7 is zero,we add zero
    addi t3 t3 0
    #store the decompressed cmd to the address of a4
    sw t3 0(a4)
    #return
    jr ra
csub_decompress:
    #reset the tmp value to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move t0 2 bit to get rs2'
    srli t2 t0 2
    #get the last 3 bit of t2 to have rs2'
    andi t2 t2 0x7
    #add 8 to t2  to change rs2' to rs2
    addi t2 t2 8
    #right move t0 7 bit to get rs1'
    srli t1 t0 7
    #get the last 3 bit of t2 to have rs1'
    andi t1 t1 0x7
    #add 8 to t1  to change rs1' to rs1
    addi t1 t1 8
    #add opcode to output cmd t3 ,opcode is 0x33
    addi t3 t3 0x33
    #left move t1 7 bit to get rd value
    slli t1 t1 7
    #add rd to output cmd
    add t3 t3 t1
    #add func3 to outputcmd t3,func3 is 0x0000
    addi t3 t3 0x0000
    #left move t1 8 bit more to get rs1 value 
    slli t1 t1 8
    #add rs1 value to output cmd t3
    add t3 t3 t1
    #left move t2 20 bit to get rs2
    slli t2 t2 20
    #add  rs2 to the output cmd
    add t3 t3 t2
    #give the func7 value to t0
    li t0 1
    #left move t0 30 bit to have the func7 of output
    slli t0 t0 30
    #add func7 to output cmd ,since func 7 is 0x40000000,we add 0x40000000
    add t3 t3 t0
    #store the decompressed cmd to the address of a4
    sw t3 0(a4)
    #return
    jr ra
csrli_decompress:
    #reset the tmp value to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move 2 bit to get the shamt[4:0]
    srli t1 t0 2
    #get the last 5 bit of t1 to get shamt[4:0]
    andi t1 t1 0x1f
    #right move 7 bit to get rd'
    srli t2 t0 7
    #get the last 3 bit of t2 to get rd
    andi t2 t2 0x07
    #add 8 to change rd' into rd
    addi t2 t2 8
    #right move 12 bit to get shamt[5]
    srli t3 t0 12
    #get the last bit of t3 to get shamt[5]
    andi t3 t3 1
    #add opcode to the output cmd,opocode is 0x13
    addi t4 t4 0x13
    #left move rd 7 bit to get the rd of the output cmd
    slli t2 t2 7
    #add the rd to the output cmd
    add t4 t4 t2
    #give the func3 value to t0
    li t0 5
    #left move t0 12 bit to have the func3 of output
    slli t0 t0 12
    #add func3 to the output cmd,func3 is 0x5000
    add t4 t4 t0
    #left move rd 8 bit more to get the rs1 of the output cmd
    slli t2 t2 8
    #add the rs1 to the output cmd
    add t4 t4 t2
    #left move shamt[4:0] 20 bit to get the shamt[4:0] of the output
    slli t1 t1 20
    #add the shamt[4:0] to the output cmd
    add t4 t4 t1
    #left move shamt[5] 25 bit to get the shamt[5] of the output cmd
    slli t3 t3 25
    #add the shamt[5] to the output cmd
    add t4 t4 t3
    #store the output cmd to the address of the a4
    # the func7 of this cmd is zero,so we add 0
    sw t4 0(a4)
    #return
    jr ra
csrai_decompress:
    #reset the tmp value to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move 2 bit to get the shamt[4:0]
    srli t1 t0 2
    #get the last 5 bit of t1 to get shamt[4:0]
    andi t1 t1 0x1f
    #right move 7 bit to get rd
    srli t2 t0 7
    #get the last 3 bit of t2 to get rd
    andi t2 t2 0x07
    #add 8 to change rd' into rd
    addi t2 t2 8
    #right move 12 bit to get shamt[5]
    srli t3 t0 12
    #get the last bit of t3 to get shamt[5]
    andi t3 t3 1
    #add opcode to the output cmd,opocode is 0x13
    addi t4 t4 0x13
    #left move rd 7 bit to get the rd of the output cmd
    slli t2 t2 7
    #add the rd to the output cmd
    add t4 t4 t2
    #give the func3 value to t0
    li t0 5
    #left move t0 12 bit to have the func3 of output
    slli t0 t0 12
    #add func3 to the output cmd,func3 is 0x5000
    add t4 t4 t0
    #left move rd 8 bit more to get the rs1 of the output cmd
    slli t2 t2 8
    #add the rs1 to the output cmd
    add t4 t4 t2
    #left move shamt[4:0] 20 bit to get the shamt[4:0] of the output
    slli t1 t1 20
    #add the shamt[4:0] to the output cmd
    add t4 t4 t1
    #left move shamt[5] 25 bit to get the shamt[5] of the output cmd
    slli t3 t3 25
    #add the shamt[5] to the output cmd
    add t4 t4 t3
    #give the func7 value to t0
    li t0 1
    #left move t0 30 bit to have the func7 of output
    slli t0 t0 30
    #add func7 to the output cmd,func7 is 0x40000000
    add t4 t4 t0
    #store the output cmd to the address of the a4
    sw t4 0(a4)
    #return
    jr ra
candi_decompress:
    #reset the tmp value to zero
    li t0 0
    li t1 0
    li t2 0
    li t3 0
    li t4 0
    li t5 0
    li t6 0
    #load value to t0 from the address of a4
    lw t0 0(a4)
    #right move 2 bit to get the shamt[4:0]
    srli t1 t0 2
    #get the last 5 bit of t1 to get shamt[4:0]
    andi t1 t1 0x1f
    #right move 7 bit to get rd
    srli t2 t0 7
    #get the last 3 bit of t2 to get rd
    andi t2 t2 0x07
    #add 8 to change rd' into rd
    addi t2 t2 8
    #right move 12 bit to get shamt[5]
    srli t3 t0 12
    #get the last bit of t3 to get shamt[5]
    andi t3 t3 1
    #add opcode to the output cmd,opocode is 0x13
    addi t4 t4 0x13
    #left move rd 7 bit to get the rd of the output cmd
    slli t2 t2 7
    #add the rd to the output cmd
    add t4 t4 t2
    #give the func3 value to t0
    li t0 7
    #left move t0 12 bit to have the func3 of output
    slli t0 t0 12
    #add func3 to the output cmd,func3 is 0x7000
    add t4 t4 t0
    #left move rd 8 bit more to get the rs1 of the output cmd
    slli t2 t2 8
    #add the rs1 to the output cmd
    add t4 t4 t2
    #left move shamt[4:0] 20 bit to get the shamt[4:0] of the output
    slli t1 t1 20
    #add the shamt[4:0] to the output cmd
    add t4 t4 t1
    #left move shamt[5] 25 bit to get the shamt[5] of the output cmd
    slli t3 t3 25
    #add the shamt[5] to the output cmd
    add t4 t4 t3
    #left move 6 bit to sign-extend
    slli t4 t4 6
    #right move 6 bit arithmatically to sign-extend
    srai t4 t4 6 
    #store the output cmd to the address of the a4
    sw t4 0(a4)
    #return
    jr ra

cjal_cj_decompress:
    # get cmd into t1
    lw t1,0(a4)
    # get imm into t1
    srli t1,t1,2
    # clear highest three bits 
    andi t1,t1,0x7ff
    # now t1 contians imm11
    # store offset into t2 and use t3 to store part of the offset
    # initialize t2 to 0
    li t2,0
    # store imm(5) into t3
    andi t3,t1,0x1
    # move t3 to make imm(5) into its position
    slli t3,t3,5
    # merge t3 into t2
    or t2,t2,t3
    # get imm(3:1) into t3 and they are in their correct position
    andi t3,t1,0xe
    # merge t3 into t2
    or t2,t2,t3
    # get imm(7) into t3
    andi t3,t1,0x10
    # move t3 to make imm(7) into its position
    slli t3,t3,3
    # merge t3 into t2
    or t2,t2,t3
    # get imm(6) into t3
    andi t3,t1,0x20
    # move t3 to make imm(6) into its position
    slli t3,t3,1
    # merge t3 into t2
    or t2,t2,t3
    # get imm(10) into t3
    andi t3,t1,0x40
    # move t3 to make imm(10) into its position
    slli t3,t3,4
    # merge t3 into t2
    or t2,t2,t3
    # get imm(9:8) into t3
    andi t3,t1,0x180
    # move t3 to make imm(9:8) into its position
    slli t3,t3,1
    # merge t3 into t2
    or t2,t2,t3
    # get imm(4) into t3
    andi t3,t1,0x200
    # move t3 to make imm(4) into its position
    srli t3,t3,5
    # merge t3 into t2
    or t2,t2,t3
    # sign extended imm(11)
    slli t3,t1,21
    # sign extend
    srai t3,t3,31
    # move t3 to make imm(11) and sign bits into their position
    slli t3,t3,11
    # merge t3 into t2
    or t2,t2,t3
    # now, t2 get the correct offset
    # get opcode into t3
    li t3,0x6f
    # use t1 to store the 32bits command and initialzie it to 0
    li t1,0
    # merge t3 into t1
    or t1,t1,t3
    # get sign bit of offset into t3
    srai t3,t2,11
    # use t4 to store new imm and initialize it to 0
    li t4,0
    # get offset(11|19:12) which are all sign bits into t3
    andi t3,t3,0x1ff
    # merge t3 into t4
    or t4,t4,t3
    # get offset(10:1) int0 t3
    andi t3,t2,0x7fe
    # move t3 to get offset(10:1) its positon
    slli t3,t3,8
    # merge t3 into t4
    or t4,t3,t3
    # get offset(20) into t3
    srli t3,t2,20
    # clear other bits 
    andi t3,t3,0x1
    # move t3 to make offset(20) its position
    slli t3,t3,19
    # merge t3 into t4
    or t4,t3,t4
    # now t4 get the correct imm
    # move t4 to make imm its position
    slli t4,t4,12
    # merge t4 into t1
    or t1,t4,t1
    # get origin cmd  into t2
    lw t2, 0(a4)
    # get the highest bit to check what cnmd it is
    srli t2,t2,15
    # /if the highest bit is zero, the cmd is cjal
    beqz t2,cjal_decompress
    # /if the cmd is cj, the src is zero, just store the new cmd
    cj_decompress: j store_cmd
    cjal_decompress:
        # /if the cmd is cjal, the src is x1, change the new cmd
        # get 0x80 into t2
        li t2,0x80
        # merge t2 into t1
        or t1,t1,t2
        # store the new cmd
        j store_cmd
    # store t1 into data array
    store_cmd:  
    sw t1,0(a4)
    # return 
    jr ra

cbeqz_cbnez_decompress:
    # get cmd into t1
    lw t1,0(a4)
    # use t2 to store the offset and initialize it to zero
    li t2,0
    # use t3 to store part of the offset
    # get offset(5) into t3
    andi t3,t1,0x4
    # move t3 to get offset(5) to its position
    slli t3,t3,3
    # merge t3 into t2
    or t2,t2,t3
    # store offset(2:1) into t3
    andi t3,t1,0x18
    # move t3 to get offset(2:1) to its position
    srli t3,t3,2
    # merge t3 into t2
    or t2,t2,t3
    # get offset(7:6) into t3
    andi t3,t1,0x60
    # move t3 to get offset(7:6) to its position
    slli t3,t3,1
    # merge t3 into t2
    or t2,t2,t3
    # new t2 gets the correct imm5
    # move t1 to make imm3 to the rightest
    srli t1,t1,10
    # get offset(4:3) into t3
    andi t3,t1,0x3
    # move t3 to get offset(4:3) into its position
    slli t3,t3,3
    # merge t3 into t2
    or t2,t2,t3
    # get offset(8) into t3
    andi t3,t1,0x4
    # move t3 to get offset(8) into its position
    slli t3,t3,6
    # merge t3 into t2
    or t2,t2,t3
    # sign extended the offset
    slli t2,t2,23
    # sign extended the offset
    srai t2,t2,23
    # now t2 get the correct offset
    
    # use t4 to store the 32bits cmd and initialize it to opcode
    li t4,0x63
    # store imm(4:1|11) into t1 and initialize it to zero
    li t1,0
    # use t3 to store part of the offset
    # sotre offset(11) into t3
	# first right shift the offset and store it into t3
    srli t3,t2,11
    # then get the 11 bit
    andi t3,t3,0x1
    # merge t3 into t1
    or t1,t1,t3
    # get offset(4:1) into t3
    andi t3,t2,0x1e
    # merge t3 into t1
    or t1,t1,t3
    # move t1 to make imm5 to its position
    slli t1,t1,7
    # merge t1 into t4
    or t4,t4,t1

    # store imm(12|10:5) into t1 and initialize it to zero
    li t1,0
    # get offset(10:5) into t3
    andi t3,t2, 0x7e0
    # move t3 to get offset(10:5) into its position
    srli t3,t3,5
    # merge t3 into t1
    or t1,t1,t3
    # get offset(12) and the signal bits into t3
    srai t3,t2,12
    # move t3 to get the offset(12) into its position
    slli t3,t3,6
    # merge t3 into t1
    or t1,t1,t3
    # move t1 to make imm7 to its position
    slli t1,t1,25
    # merge t1 into t4
    or t4,t4,t1
    # now t4 has no func3, rs1 and rs2
    # load origin cmd into t1
    lw t1,0(a4)
    # get rs1 into t2
    # first shift the origin cmd to make lower bits zero
    srli t2,t1,7
    # then, get rs1
    andi t2,t2,0x7
    # change the popular register into origin register by adding it by 8
    addi t2,t2,8
    # move t2 to make rs1 its position
    slli t2,t2,15
    # merge t2 into t4
    or t4,t4,t2
    # now t4 has no funct3
    # get cmd(13) into t2 to check which compressed cmd it is
    # first right shift origin to make the lower bits zero
    srli t2,t1,13
    # then clear other bits
    andi t2,t2,0x1
    # /if cmd(13) is 0, the cmd is cbeqz
    beqz t2,cbeqz_decompress
    # /else, it is cbneqz
    cbnez_decompress: 
        # load 0x1 into t2
        li t2,0x1
        # move t2 to get func3 its positon
        slli t2,t2,12
        # merge t2 into t4
        or t4,t4,t2
        # then store the new cmd
        sw t4,0(a4)
        # return 
        jr ra
    cbeqz_decompress:
        # /if the cmd is cbeqz, funct3 is 0, need to do nothing, just store
        sw ,t4,0(a4)
        # return 
        jr ra

# codes above is used to decompress the data in a4
# ================================================================================================================ #

# ================================================================================================================ #

# ================================================================================================================ #
# codes below are for output handler
output_handler:
    # get 4 words in the stack
    addi sp,sp,-16
    # store return address
    sw ra,0(sp)
    # store the state array
    sw a3,4(sp)
    # store data array
    sw a2,8(sp)
    # store the num line
    sw a1,12(sp)
    # get i into t0
    li t0,0
    # move line num from a1 into t1
    mv t1,a1
    loop_output_cmd:
        # load 0x80000000 into t3 as a signal
        lui t3,0x80000
        # /if i is equal to the line num, return 
        beq t0,t1,return_output_cmd
        # get cmd into t2
        lw t2,0(a2)
        loop_print_cmd:
            # /if already print the cmd, return 
            beqz t3,continue_output_cmd
            # get (31-i)bit into t4
            and t4,t2,t3
            # /if t4 is 0, print zero
            beqz t4,print_zero
            # print one
            print_one: li a1,1
            # print 
            jal print_int
            # continue
            j continue_print_cmd
            # print zero
            print_zero: li a1,0
            # print
            jal print_int
            # continue
            j continue_print_cmd

continue_print_cmd:
    # right shift the signal
    srli t3,t3,1
    # jump back the loop
    j loop_print_cmd

continue_output_cmd:
    # i++
    addi t0,t0,1
    # move data array one word forward
    addi a2,a2,4
    # print a new line
    jal print_newline
    # jump back the loop
    j loop_output_cmd

return_output_cmd:
    # store return address
    lw ra,0(sp)
    # store the state array
    lw a3,4(sp)
    # store data array
    lw a2,8(sp)
    # store the num line
    lw a1,12(sp)
    # reset the stack
    addi sp,sp,16
    # return 
    jr ra
# codes above are for output handler
# ================================================================================================================ #

# ================================================================================================================ #
# codes below are for change offset

change_offset:
    # get 
    addi sp,sp,-16
    # store return address
    sw ra,0(sp)
    # store data array
    sw a2,4(sp)
    # store state array
    sw a3,8(sp)
    # store line num
    sw a1,12(sp)
    # let t0 be i and initialize it to 0
    li t0,0
    loop_change_offset:
        # /if i = line num,return
        beq t0,a1,return_change_offset
        # get state array into t1
        lw t1,0(a3)
        # get cmd address into a4
        mv a4,a2
        # get 2 into t2 for feather check
        li t2,2
        # get 3 into t3 for feather check
        li t3,3
        # get 4 into t4 for feather check
        li t4,4
        # get 5 into t5 for feather check
        li t5,5
        # /if the state is 2, the cmd is branch
        beq t1,t2,call_get_b_offset
        # /if the state is 3, the cmd is branch
        beq t1,t3,call_get_b_offset
        # /if the state is 4, the cmd is jump
        beq t1,t4,call_get_j_offset
        # /if the state is 5, the cmd is jump
        beq t1,t5,call_get_j_offset
        # /else continue
        j continue_change_offset

        # get b offset into a5
        call_get_b_offset: jal get_b_offset
        # get state array pointer into a6
        mv a6,a3
        # calc new offset and store it into a5
        jal calc_new_offset
        # get new branch cmd and store it into data array
        jal get_new_b_cmd
        # continue
        j continue_change_offset

        # get j offset into a5
        call_get_j_offset: jal get_j_offset
        # get state array pointer into a6
        mv a6,a3
        # calc new offset and store it into a5
        jal calc_new_offset
        # get new jump cmd and store it into data array
        jal get_new_j_cmd
        # continue 
        j continue_change_offset


get_new_j_cmd:
    # get cmd into t1
    lw t1,0(a4)
    # store correct offset into t2 and initialize it to zero
    li t2,0
    # use t3 to store part of the offset
    # get offset(19:12) into t3
    srli t3,a5,12
    # clear the higher bits
    andi t3,t3,0xff
    # merge t3 into t2
    or t2,t2,t3
    # get offset(11) into t3
    srli t3,a5,11
    # clear other bits
    andi t3,t3,0x1
    #move offset(11) into its position
    slli t3,t3,8
    # merge t3 into t2
    or t2,t2,t3
    # get offset(10:1) into t3
    srli t3,a5,1
    # clear other bits
    andi t3,t3,0x3ff
    # move offset(10:1) into their position
    slli t3,t3,9
    # merge t3 into t2
    or t2,t2,t3
    # get offset(20) and the signal bits into t3
    srai t3,a5,20
    # move offset(20) into its position
    slli t3,t3,19
    # merge t3 into t2
    or t2,t3,t2
    # now t2 get the correct imm
    # left shift t1(cmd) by 20 bits to clear the imm
    slli t1,t1,20
    # right shift t1(cmd) by 20 bits to clear the imm
    srli t1,t1,20
    # move t2 to make the imm into its position
    slli t2,t2,12
    # merge t2 into t1
    or t1,t1,t2
    # store the new cmd into a4
    sw t1,0(a4)
    # return 
    jr ra

get_new_b_cmd:
    # get cmd into t1
    lw t1,0(a4)
    # store imm(4:1|11) into t2 and initialize it to zero
    li t2,0
    # use t3 to store part of the offset
    # sotre offset(11) into t3
	# first right shift the offset and store it into t3
    srli t3,a5,11
    # then get the 11 bit
    andi t3,t3,0x1
    # merge t3 into t2
    or t2,t2,t3
    # get offset(4:1) into t3
    andi t3,a5,0x1e
    # merge t3 into t2
    or t2,t2,t3
    # load clear signal 0xfffff07f into t4
    # first load high bits
    lui t4, 0xfffff
    # then load 0x7f into t5
    li t5,0x7f
	# finally, merge t5 into t4
    or t4,t4,t5
    # clear coresponding positon of origin cmd
    and t1,t1,t4
    # move t2 to make imm5 to its position
    slli t2,t2,7
    # merge t2 into t1
    or t1,t2,t1
    # store imm(12|10:5) into t2 and initialize it to zero
    li t2,0
    # get offset(10:5) into t3
    andi t3,a5, 0x7e0
    # move t3 to get offset(10:5) into its position
    srli t3,t3,5
    # merge t3 into t2
    or t2,t2,t3
    # get offset(12) and the signal bits into t3
    srai t3,a5,12
    # move t3 to get the offset(12) into its position
    slli t3,t3,6
    # merge t3 into t2
    or t2,t2,t3
    # move t2 to make imm7 to its position
    slli t2,t2,25
    # clear the imm7 of origin cmd
    slli t1,t1,7
    # clear the imm7 of origin cmd
    srli t1,t1,7
    # merge t2 into t1
    or t1,t1,t2
    # store t1 into data array
    sw t1,0(a4)
    # return 
    jr ra


calc_new_offset:
    # let t2 be the new offset and initialize it to 0
    li t2,0
    # load 2 into t2 for feather check
    li t3,2
    # load 4 into t4 for feather check
    li t4,4
    # /if the offset is more than 0, loop down
    bgt a5,zero,loop_down_calc_new_offset
    # /else, loop up
    # move the state array one word after since jump up donnot consider the cmd it
    addi a6,a6,-4
    # loop up
    j loop_up_calc_new_offset
    loop_up_calc_new_offset:
        # /if offset is 0, return
        beqz a5,return_calc_new_offset
        # add new offset by -4
        addi t2,t2,-4
        # get state array into t5
        lw t5,0(a6)
        # /if the state is 0, increase the offset by 4
        beqz t5,increase_4
        # /if the state is 2, increase the offset by 4
        beq t5,t3,increase_4
        # /if the state is 4, increase the offset by 4
        beq t5,t4,increase_4
        j increase_2
        
    loop_down_calc_new_offset: 
        # /if offset is 0, return
        beqz a5,return_calc_new_offset
        # add new offset by 4
        addi t2,t2,4
        # get state array into t5
        lw t5,0(a6)
        # /if the state is 0, decrease the offset by 4
        beqz t5,decrease_4
        # /if the state is 2, decrease the offset by 4
        beq t5,t3,decrease_4
        # /if the state is 4, decrease the offset by 4
        beq t5,t4,decrease_4
        j decrease_2

decrease_2:
    # decrease offset by 2
    addi a5,a5,-2
    # continue
    j continue_down_calc_new_offset

decrease_4:
    # decrease offset by 4
    addi a5,a5,-4
    # continue
    j continue_down_calc_new_offset

continue_down_calc_new_offset:
    # move state array one word forword
    addi a6,a6,4
    # continue
    j loop_down_calc_new_offset

increase_2:
    # increase offset by 2
    addi a5,a5,2
    # continue
    j continue_up_calc_new_offset

increase_4:
    # increase offset by 4
    addi a5,a5,4
    # continue
    j continue_up_calc_new_offset

continue_up_calc_new_offset:
    # move state array one word afterword
    addi a6,a6,-4
    # continue
    j loop_up_calc_new_offset


return_calc_new_offset:
    # move new_offset into a5
    mv a5,t2
    # return 
    jr ra
continue_change_offset:
    # i++
    addi t0,t0,1
    # move data array one word forward
    addi a2,a2,4
    # move state array one word forward
    addi a3,a3,4
    # continue
    j loop_change_offset
get_j_offset:
    # get cmd into a5
    lw a5,0(a4)
    # get all imm into t1 by right shift the cmd by 12 bits(since the signal bit is the first,keep it)
    srai t1,a5,12
    # use t2 to store the offset and t3 to store the part of the offset
    # store imm(20|10:1) into t2 by right shift the imm by 9 bits
    srli t2,t1,9
    # get imm(10:1) into t2 by and a signal on t2
    andi t2,t2,0x3ff
    # left shift t2 by 1 bit to get the least significant bit zero,now the t2 get imm(10:0s)
    slli t2,t2,1
    # get 11 bit into t3 by right shift the imm by 8 bits
    srli t3,t1,8
    # make other bits zero
    andi t3,t3,0x1
    # left shift t3 by 11 bits to make 11 bit inot its position
    slli t3,t3,11
    # merge it into t2 and now t2 get imm(11:0)
    or t2,t2,t3
    # get imm(19:12) by and a signal onto t1 and store it into t3
    andi t3,t1,0xff
    # left shift t3 by 12 bits to make bit 19:12 their position
    slli t3,t3,12
    # merge t3 into t2 and now t2 get imm(19:0)
    or t2,t2,t3
    # get the signal bits and imm (20) into t3 by right shifting t1 by 19 bits
    srai t3,t1,19
    # left shift t3 to make bit 20 its position
    slli t3,t3,19
    # merge t3 into t2
    or t2,t2,t3
    # now t2 get the offset, store it into a5
    mv a5,t2
    # return 
    jr ra

get_b_offset:
    # get cmd into a5
    lw a5,0(a4)
    # get imm(12|10:5) into t1 with sign bits
    srai t1,a5,25
    # left shift t1 by 5 bits to save space for imm(4:1|11)
    slli t1,t1,5
    # get imm(4:1|11) into t2
    srli t2,a5,7
    # and t2 with a signal to clear high bits
    andi t2,t2, 0x1f
    # merge t2 with t1 into t1 to make t1 the whole imm
    or t1,t1,t2
    # use t2 for offset and t3 for part of the imm
    # get imm(4:1) into t2 by and a signal 
    andi t2,t1,0x1e
    # get imm(10:5) into t3
    andi t3,t1,0x7e0
    # merge t3 into t2
    or t2,t2,t3
    # get imm(11) into t3
    andi t3,t1,0x1
    # left shift t3 by 11 bits to get 11th bit into right position
    slli t3,t3,11
    # merge t3 into t2
    or t2,t2,t3
    # get imm(31:12) into t3
    srai t3,t1,11
    # left shift t3 by 12 bits to get 11th bit into right position
    slli t3,t3,12
    # merge t3 into t2
    or t2,t2,t3
    # now t2 get the correct offset
    # move t2 into a5
    mv a5,t2
    # return 
    jr ra


return_change_offset:
    # load back return address
    lw ra,0(sp)
    # load back data array
    lw a2,4(sp)
    # load back state array
    lw a3,8(sp)
    # load back line num
    lw a1,12(sp)
    # reset stack
    addi sp,sp,16
    # return 
    jr ra

# codes above are for change offset
# ================================================================================================================ #

# ================================================================================================================ #
#codes below are some helpful functions

print_int:
    # ID for print int 
    li a0,1
    # print 
    ecall
    # return 
    jr ra

print_newline:
    # load a1 with '\n'
    li a1,'\n'
    # ID for print char
    li a0,11
    # print
    ecall 
    # return
    jr ra

print_ints:
    # move sp one word below
    addi sp,sp,-4
    # store the return register
    sw ra,0(sp)
    # print variables in a loop
    print_loop:
        # /if print all variables, return
        beqz a3,return_print_ints
        # decrease the length by one
        addi a3,a3,-1
        # load the variable which is need printing into a1
        lw a1,0(a2)
        # call print
        jal print_int
        # print a new line
        jal print_newline
        # move a2(data pointer) one word forward
        addi a2,a2,4
        # continue
        j print_loop
return_print_ints:
    # get the return adress from stack
    lw ra,0(sp)
    # reset sp
    addi sp,sp,4
    # return 
    jr ra

# codes above are some helpful functions
# ================================================================================================================ #

