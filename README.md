input_file --> input_handler --> compress --> output_handler --> output_file

**translator**  : check input and handle with the file opening and closing 

**utils**       : handle data input and data output 

**compression** : compress input commands and store compressed commands

main step of compression:
 * check command type 
 * check whether it can be compressed and mark them as compressible, icompressible, unsure
 * compress commands which are compressible
 * compress commands which are unsure 
