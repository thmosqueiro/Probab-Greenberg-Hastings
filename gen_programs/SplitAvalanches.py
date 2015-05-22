import os


for input_name in os.listdir('./'):
    
    if input_name.endswith(".data"):
        
        print " Separando avalanches de " + input_name
        
        inputfile = open(input_name, 'r')
        output_T_file = open('Time_' + input_name, 'w')
        output_S_file = open('Size_' + input_name, 'w')
        
        for line in inputfile:
            
            Sspted = line.split("Size:")    
            Tspted = line.split("Tempo:")
            
            if ( len(Tspted) > 1 ) :
                output_T_file.write( str( int(Tspted[1]) ) + "\n")
            if ( len(Sspted) > 1 ) :
                output_S_file.write( str( int(Sspted[1]) ) + "\n")
        
        
        inputfile.close()
        output_T_file.close()
        output_S_file.close()
