#   Project 1.2 RVC instructions to RISC-V instructions in RISC-V

#   You can suppose you will get the number of machine codes you should translate from 
#   "lines_of_rvc_codes" and get the 16 or 32 bits-machine codes from "rvc_codes". 
#   Both of them are stored in file "input.S"

#   When we test your project, we will change the input files content.

# DO NOT REMOVE THIS LINE WHEN SUBMITTING
.import input.S

.text
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

# ================================================================#
# codes below is used to handle input
input_handler:
    # move the stack three words below
    addi sp,sp,-12
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
read_in:
    # /if already read in all commands, return
    beqz t4,return_read_in
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
    j read_in
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
    j read_in
return_read_in:
    # load back the line num
    lw a1,8(sp)
    # load back the state pointer
    lw a3,4(sp)
    # load back the data pointer
    lw a2,0(sp)
    # reset sp
    addi sp,sp,12
# codes above is used to handle input data
# =================================================================#


de_compress:
    # let t0 be i and initialize it to zero
    li t0,0
    # get three words in the stack
    addi sp,sp,-12
    # store len num
    sw a1,8(sp)
    # store data array
    sw a2,4(sp)
    # store state array
    sw a3,0(sp)
# in this function, loop through all cmds to check their types and de-compress them
check_type:
    # /if i== len num , return
    beq t0,a1,return_check_type
    # move state array one word forward
    addi a3,a3,4
    # move data array one word forward
    addi a2,a2,4
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
return_check_type:
    # load back len num
    lw a1,8(sp)
    # load back data array
    lw a2,4(sp)
    # load back state array
    lw a3,0(sp)
    # reset stack
    addi sp,sp,12
    # Exit with error code 0
    add	   a1, x0, x0
    addi   a0, x0, 17 
    ecall

#==================================================================#
#codes below are some helpful functions

# this function is used to print an int variable in a1 
print_int:
    # ID for print int 
    li a0,1
    # print 
    ecall
    # return 
    jr ra

# this function is used to print a chr variable in a1
print_char:
    # ID for print char
    li a0,11
    # print
    ecall 
    # return
    jr ra

# this function is used to print all int variables whose head is in a2 with length in a3
# after each number, there is a new line
print_ints:
    # move sp one word below
    addi sp,sp,-4
    # store the return register
    sw ra,0(sp)
    # print variables in a loop
    print_loop:
        # if print all variables, return
        beqz a3,return_print_ints
        # decrease the length by one
        addi a3,a3,-1
        # load the variable which is need printing into a1
        lw a1,0(a2)
        # call print
        jal print_int
        # print a new line
        li a1,'\n'
        # call print
        jal print_char
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
# ==================================================================#