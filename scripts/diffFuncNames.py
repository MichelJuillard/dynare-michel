import os
import string

for dirname, dirnames, filenames in os.walk('../matlab'):
    for filename in filenames:
        filename = string.strip(filename)

        if filename[-2:] != '.m' or filename == 'msstart2.m' or filename == 'msstart_setup.m' or filename == 'qmc_sequence.m':
            continue

        fullfilename = os.path.join(dirname, filename)
        f = open(fullfilename, 'r')
        funcDef = ''
        inComment = False
        while True:
            funcDef += f.read(1)
            if funcDef[-1:] == '%':
                inComment = True

            if inComment:
                if funcDef[-1:] == '\n':
                    inComment = False
            else:
                if funcDef[-1:] == '(':
                    break
        f.close()

        spliteq = string.rsplit(funcDef, '=')
        if len(spliteq) == 1:
            spliteq = string.rsplit(funcDef, 'function ')

        spliteq = spliteq.pop()
        spliteq = string.strip(spliteq, '. ')
        spliteq = string.strip(spliteq, '\n ')
        spliteq = string.strip(spliteq, '( ')

        if filename[:-2] != spliteq:
            print fullfilename + ': ' + spliteq

