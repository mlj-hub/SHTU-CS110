# code architecture

input data -> input_handler ->  de-compress -> output_handler -> ouput file

# input_handler
* store all commands in a array(**data array**) whose head pointer is stored in **a2**
* another array is used to indicate cmds' length: 1 for 16 bits and 0 for 32 bits. The head pointer of flag array(**state array**) is in **a3**
* The data length(**line_num**) is in **a1**

# de_compress
## NOTE:
In this part, there are three steps: **check type** , **de_compress** and **change offset**. The main fuction is **de_compress** and it will call **check_type** first, in which we will call sub-functions to **de_compress** all commands. Finally we call **change_offset**.

### check_type
* In this part, we need to check all cmds and change the value in **state array** to remark them:
  * 0: 32bits command which is not jump or branch  
  * 1: 16bits command which is not jump or branch
  * 2: 32bits command which is jump or branch
  * 3: 16bits command which is jump or branch
* The main function is check_type. This function will handle the check loop and call check functions according to opcode. For 32bits cmds, just check whether it is jump or branch and change the **state array**, whic will use **j_b_check_32**. For 16bits cmds, call the following functions(in these functions, **state array** will be changed):
  * check_00: for those cmds whose opcode is 00
  * check_10: for those cmds whose opcode is 10
  * check_01: for those cmds whose opcode is 01
* Once determining it's corresponding command in 32 bits, we de-compress part of them by calling de-compress functions seperately in the check funcitons, which follows the rules below:
  * for commands that are not in branch or jump, we directly de-compress it or leave it unchanged (for 32bits commands)
  * for jump and branch commands, leave it unchanged
  * need to change the **state array** for all commands

### de-compress
* In this part, we de-compress the cmds and change the value in **data array**
* we use a set of sub-functions to de-compress all commands. They will be called in **check_type** and the parameter is the **address** of the cmd which is stored in **a4**. These functions need to de-compress the value and store it back into the address. Some of the functions are listed below:
  * cadd_decompress
  * cmv_decompress
  * cjr_decompress
  * cjalr_decomress
  * ......

### change_offset
* In this part, we change offset of the jump and branch commands
* We use two sub-functions in this part :  **change_offset_j** and **change_offset_b**. The parameters are: **state array**(**a3**),**data array**(**a2**), **line_num**(**a1**). They need to change the offset and store it back into the address. They will be called in function **change_offset**. 

# output_handler
* output all values in the **data array**