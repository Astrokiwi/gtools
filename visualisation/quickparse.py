# plot_parameters = list()
#allthings = list()


with open("toparse.dat") as f:
#     current_parameters = None
#     col = None
    for line in f:
        if ( line[0]!="#"):
            if ( "plot_thing" in line ):
                subStart = line.index("'")
                subEnd = line.index("'",subStart+1)
                plotLabel = line[subStart+1:subEnd]
    #             print('plotLabel["{}"] = '.format(line[subStart+1:subEnd]))
    #             subStart = line.index("'")
    #             subEnd = line.index("'",subStart+1)
                col = 1
    #             current_parameters = dict()
    #             current_parameters["plot_name"] = line[subStart+1:subEnd]
    #             plot_parameters.append(current_parameters)
            elif ( "cols==2" in line ):
                col = 2
            elif ( "makesph_plot" in line and col==1 ):
                subStart = line.index('"')
                subEnd = line.index('"',subStart+1)
                plotFormat = line[subStart-1:subEnd+1]
                print('plotLabel["{}"] = {}'.format(plotLabel,plotFormat))
            
    #             splut = line.split(",")
    #             print(splut)