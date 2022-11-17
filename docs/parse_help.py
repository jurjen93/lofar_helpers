
f = open('../h5_merger.py', 'r')
lines = f.readlines()

output = {}

for l in lines:
    if 'help' in l:
        llist = l.split(',')
        print(llist)

        for ll in llist:
            # param, help = None, None
            if '--' in ll:
                # print(ll)
                if ' ' not in ll:
                    param = ll
                    # print(param)
            elif 'parser.add_argument' in ll:
                param = ll.replace('parser.add_argument(', '')
                # print(param)
            if 'help' in ll:
                help = l.replace(')','').replace('(','').replace('help=','')
                # print(help)

            # if param is not None and help is not None:
            #     output.update({ll:help})

# print(output)



