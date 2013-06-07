#!/usr/bin/python

# Suppose you got the results from bench.cc
# We here parse the result file and produce a PDF table
# in which the results are summarized

import sys
from subprocess import call

if(len(sys.argv) != 2):
    print "Usage : ", sys.argv[0], " <data file>"
    sys.exit()

file = open(sys.argv[1],'r')
latex_filename = sys.argv[1]+".tex"
latex_file = open(latex_filename,'w')

# Header of the latex file
latex_file.write('\\documentclass[a4paper,12pt]{article}\n')
latex_file.write('\\textwidth17cm \\textheight24.0cm \\topmargin-2.1cm\n')
latex_file.write('\\evensidemargin-0.5cm \\oddsidemargin-0.5cm \\parindent0.0cm\n')
latex_file.write('\\unitlength1.0mm\n')
latex_file.write('\\usepackage{graphicx}\n')
latex_file.write('\\usepackage{longtable}\n')
latex_file.write('\\usepackage{fancyhdr}\n')
latex_file.write('\\usepackage{graphics}\n')
latex_file.write('\\setlength{\\headheight}{15.2pt}\n\n')
latex_file.write('\\begin{document}\n')

column_formats = ['%s','%s','%s', '%s','%d','%.2e','%.2e','%.2e','%.2e','%.2e','%s']
column_names = ['RNG','Function', 'Dim', 'Algorithm', 'Nb trials', 'Error (mean)' ,'Error (std)', 'FE (mean)', 'Log progress','Best min','Success rate']

latex_file.write('\\begin{longtable}{')
for i in range(len(column_names)):
    latex_file.write('c')
latex_file.write('}\n')

for i in range(len(column_names)-1):
    latex_file.write(column_names[i]+" & ")
latex_file.write(column_names[len(column_names)-1]+" \\\\ \n")

for l in file.readlines():
    #First we split the line with the ";" separator
    contents = l.split(';')
    for i in range(len(column_names)-2):
        if(column_formats[i] == '%d' or column_formats[i] == '%.2e'):
            latex_file.write(column_formats[i] % float(contents[i].split("=")[1])+"& ")
        else:
            latex_file.write(column_formats[i] % (contents[i].split("=")[1])+"& ")

    latex_file.write(contents[len(column_names)-2].split("=")[1]+"\\\\ \n")
    
latex_file.write('\\end{longtable}\n')

latex_file.write('\\end{document}\n')
latex_file.close()


print "Compiling latex .... "
call(["pdflatex",latex_filename])
print "\n\n"
print "Done ... file %s.pdf generated"%sys.argv[1]


# We may also draw the table as ....
# import matplotlib.pylab as plt


# fig = plt.figure()
# ax=plt.gca()
# y=[1,2,3,4,5,4,3,2,1,1,1,1,1,1,1,1]
# col_labels=['col1','col2','col3']
# row_labels=None #['row1','row2','row3']
# table_vals=[[11,12,13],[21,22,23],[31,32,33]]
# # the rectangle is where I want to place the table
# the_table = plt.table(cellText=table_vals,
#                   colWidths = None,
#                   rowLabels=row_labels,
#                   colLabels=col_labels,
#                   loc='center right')
# ax.add_table(the_table)
# plt.axis('off')
# plt.savefig("test.pdf",bbox_inches='tight')
# plt.show()

# Or generate an HTML file with the results
