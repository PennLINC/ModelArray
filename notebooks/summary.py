'''
This is for summarizing the results
'''

import argparse  
import os 

if __name__ == '__main__'   # main function
    parser = argparse.ArgumentParser(description="parse args")
    parser.add_argument('--folder',type=str,help='folder of outputs')

    args=parser.parse_args()

    # re-define the argument names:
    folder = args.folder

    # folder exists?
    if os.path.isdir(folder) == FALSE:
        print("the folder does not exist!")
        exit(1)

    # output.txt exists in this folder?
    fn_output_txt = folder + "/" + "output.txt"
    if os.path.isfile(fn_output_txt) == FALSE:
        print("the output.txt does not exist in this folder!")
        exit(1)
    else:
        # call function to extract related values:
        
