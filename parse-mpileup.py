#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_


# Version information START --------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard

Version-01:
  2019-05-24 解析mpileup文件；

Version-02:
  2019-06-04 增加多核运算功能

Version-03:
  2019-08-01 修复bug，当^X出现时会导致程序崩溃

Version-04:
  2019-08-02 修复bug，对insertion, deletion的信号能够准确识别

E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------

# Learning Part START --------------------------------------------------------
LEARNING_PART = \
"""
Input samtools mpileup file
"""
# Learning Part END-----------------------------------------------------------

import argparse
import gzip
import sys
import os 
import string
import random
import time
from multiprocessing import Process

###############################################################################
# function part 
###############################################################################
def split_file_and_make_temp(input_filename, threads_num=1,temp_dir=None):
    """
    """
    #--------------------------------------------------->>>>>>
    # set temp dir 
    #--------------------------------------------------->>>>>>        
    if temp_dir is None:
        temp_dir = os.path.dirname(input_filename)
    
    #--------------------------------------------------->>>>>>
    # get input basename 
    #--------------------------------------------------->>>>>>    
    input_file_basename = os.path.basename(input_filename)

    #--------------------------------------------------->>>>>>
    # creat temp file name and open file 
    #--------------------------------------------------->>>>>>    
    #>>>>>> LOG 
    sys.stderr.write("Making temp files... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    temp_filename_list = []
    temp_file_list = []
    
    for index in range(threads_num):
        # make filename 
        temp_file_basename = "temp_" + input_file_basename + "." + str(index) + "." +  "".join(random.sample(string.ascii_letters + string.digits, 16))
        temp_file_name = os.path.join(temp_dir , temp_file_basename)
        temp_filename_list.append(temp_file_name)
        
        # open file
        temp_file = open(temp_file_name, "w") 
        temp_file_list.append(temp_file)
    
    #--------------------------------------------------->>>>>>
    # counting input file line number
    #--------------------------------------------------->>>>>>   
    #>>>>>> LOG 
    sys.stderr.write("Counting input file... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    
    total_input_line_num = 0
    
    if ".gz" == input_filename[-3:]:
        input_file = gzip.open(input_filename,"r")
    else:
        input_file = open(input_filename,"r")

    for line in input_file:
        total_input_line_num += 1 

    if (total_input_line_num % threads_num == 0):
        each_file_line_num = (total_input_line_num // threads_num)
    else:
        each_file_line_num = (total_input_line_num // threads_num) + 1

    input_file.close() 

    #>>>>>> LOG 
    sys.stderr.write("Done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))  

    #--------------------------------------------------->>>>>>
    # split temp files
    #--------------------------------------------------->>>>>> 
    # write into output filename 
    if ".gz" in input_filename:
        input_file = gzip.open(input_filename,"r")
    else:
        input_file = open(input_filename,"r")

    for index, line in enumerate(input_file):
        file_index = index // each_file_line_num 
        temp_file_list[file_index].write(line)

    # close output filename 
    input_file.close()
    [temp_file.close() for temp_file in temp_file_list]

    #>>>>>> LOG     
    sys.stderr.write("Make temp files done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    
    return(temp_filename_list)



def run_part(INPUT_FILE_PATH, OUTPUT_FILE_PATH="Stdout", MUT_NUM=0):
    """
    input_filename
    output_filename
    """
    
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # open input file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ".gz" == INPUT_FILE_PATH[-3:]:
        INPUT_FILE = gzip.open(INPUT_FILE_PATH,"r")
    else:
        INPUT_FILE = open(INPUT_FILE_PATH,"r")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # open output file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if OUTPUT_FILE_PATH == "Stdout":
        OUTPUT_FILE = sys.stdout
    elif ".gz" == OUTPUT_FILE_PATH[-3:]:
        OUTPUT_FILE = gzip.open(OUTPUT_FILE_PATH,"w")
    else:
        OUTPUT_FILE = open(OUTPUT_FILE_PATH,"w")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # main part
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # parse mpileup file
    for line in INPUT_FILE:
        data = line.strip().split('\t')
        chr_name = data[0]        
        bp = data[1]
        bases = data[4].upper()
        ref = data[2].upper()

        types = {'A':0,'G':0,'C':0,'T':0,'-':[],'+':[],'Not':[]}

        # coverage > 0 
        i = 0
        while i < len(bases):
            base = bases[i]

            if base == '^' :
                i += 2
            
            elif base == '$':
                i += 1 

            elif base == '-':
                i += 1
                del_str_num = ""
                
                while bases[i].isdigit():
                    del_str_num += bases[i]
                    i += 1
                
                delNum = int(del_str_num)
                delSeq = ""
                for a in range(delNum):
                    delSeq += bases[i]
                    i += 1
                types['-'].append(delSeq)

            elif base == '*':
                types['-'].append(bases[i])
                i += 1

            elif base == '+':
                i += 1
                add_str_num = ""
                
                while bases[i].isdigit():
                    add_str_num += bases[i]
                    i += 1
                    
                addNum = int(add_str_num)
                addSeq = ""
                for a in range(addNum):
                    addSeq += bases[i]
                    i += 1
                types['+'].append(addSeq)
                
            elif (base == '.') or (base == ','):
                types[ref] += 1
                i += 1

            else:
                if types.has_key(base):
                    types[base] += 1
                else:
                    types['Not'].append(base)
                i += 1

        adds = '.'
        adds_count = 0
        if len(types['+']) > 0:
            adds = ','.join(types['+'])
            adds_count = len(types['+'])
            
        dels = "."
        dels_count = 0
        if len(types['-']) > 0:
            dels = ",".join(types['-'])
            dels_count = len(types['-'])

        amb = '.'
        amb_count = 0
        if len(types['Not']) > 0:
            amb = ','.join(types['Not'])
            amb_count = len(types['Not'])

        # get other mutation number 
        if types.has_key(ref):
            mut_num = types["A"] + types["T"] + types["G"] + types["C"] - types[ref]
        else:
            mut_num = 0

        out_list = [
            chr_name,
            bp,
            ref,
            types['A'],
            types['G'],
            types['C'],
            types['T'],
            dels_count,
            adds_count,
            amb_count,
            dels,
            adds,
            amb,
            mut_num
        ]

        # make a filter 
        if mut_num >= MUT_NUM:
            OUTPUT_FILE.write("\t".join(map(str,out_list)) + "\n")


    INPUT_FILE.close()
    OUTPUT_FILE.close()


def merge_out_bmat(temp_filename_list, output_filename, remove_temp=True):
    """
    """
    # open output file
    if output_filename == "Stdout":
        output_file = sys.stdout
    elif ".gz" == OUTPUT_FILE_PATH[-3:]:
        output_file = gzip.open(output_filename,"w")
    else:
        output_file = open(output_filename,"w")
    
    # write header  
    header = [
        "chr_name",
        "chr_index",
        "ref_base",
        "A",
        "G",
        "C",
        "T",
        "del_count",
        "insert_count",
        "ambiguous_count",
        "deletion",
        "insertion",
        "ambiguous",
        "mut_num"
    ]
    output_file.write("\t".join(header) + "\n")    


    # open temp bmat file 
    for temp_filename in temp_filename_list:
        # write file 
        with open(temp_filename + ".bmat","r") as temp_file:
            for line in temp_file:
                output_file.write(line)
    
    # remove temp files
    if remove_temp:
        for temp_filename in temp_filename_list:
            os.remove(temp_filename)
            os.remove(temp_filename + ".bmat")

###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="convert mpileup file to info file")

    parser.add_argument("-i", "--Input",
        help="samtools mpileup format file",required=True)

    parser.add_argument("-o", "--Output",
        help="Output parsed file",default="Stdout")

    parser.add_argument("-p", "--Threads",
        help="Multiple threads number, default=1",default="1")

    parser.add_argument("-n", "--MutNum",
        help="Only contain mutation info go to the output, set 0 mean output all site, default=0",default="0")

    parser.add_argument("--TempDir",
        help="Where to keep temp files, default is the same dir with --Input",default=None)


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load the paramters
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ARGS = parser.parse_args()

    INPUT_FILE_PATH = ARGS.Input
    OUTPUT_FILE_PATH = ARGS.Output

    THREADS_NUM = int(ARGS.Threads)
    MUT_NUM = int(ARGS.MutNum)

    TEMP_DIR = ARGS.TempDir

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # make temp files
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    temp_filename_list = split_file_and_make_temp(input_filename=INPUT_FILE_PATH, threads_num=THREADS_NUM,temp_dir=TEMP_DIR)
    
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # multiple processing part 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sys.stderr.write("Parsing files... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    procs_list = []

    for index,temp_filename in enumerate(temp_filename_list):
        out_filename =  temp_filename + ".bmat"
        
        # add multiprocess
        sub_proc = Process(target=run_part, args=(temp_filename,out_filename,MUT_NUM,))
        procs_list.append(sub_proc)
        sub_proc.start()

    for sub_proc in procs_list:
        sub_proc.join()

    sys.stderr.write("Done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))        

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # merge output files
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    #>>>>>> LOG  
    sys.stderr.write("Merging files... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    merge_out_bmat(temp_filename_list=temp_filename_list, output_filename=OUTPUT_FILE_PATH, remove_temp=True)

    #>>>>>> LOG  
    sys.stderr.write("Done!... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))



    

    









# 2019-08-02