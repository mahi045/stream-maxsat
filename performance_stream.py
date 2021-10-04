import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--i', help='Description', required=True)
parser.add_argument('-a','--a', help='Description', required=True)
args = parser.parse_args()

file_name = args.i 
algorithm = args.a
file_name_open_wbo = 'assignment_' + args.i
if algorithm == "streaming":
    file_name_open_wbo = "result_streaming_" + file_name
elif algorithm == "sampling":
    file_name_open_wbo = "result_sampled_" + file_name


assignment = None
file_name_open_wbo_f = open(file_name_open_wbo, 'r')
assignment_line = None
for line in file_name_open_wbo_f:
    if line.startswith("v "):
        assignment_line = line

assignment = assignment_line.split()[1:]

input_file_f = open(file_name, 'r')
first_line = input_file_f.readline()
hard_clause_sign = int(first_line.split()[-1])

os.system('rm -f assess_{0}'.format(file_name))
os.system('cp {0} assess_{0}'.format(file_name))
write_file = open('assess_{0}'.format(file_name), 'a')
for var in assignment:
    write_file.write('{0} {1} 0 \n'.format(hard_clause_sign, var))

write_file.close()
os.system('./open-wbo_static {0} -no-print-model -cpu-lim={1}'.format('assess_{0}'.format(file_name),
    10000))

os.system('rm -f assess_*')

