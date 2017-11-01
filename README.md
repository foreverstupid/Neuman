# Neuman
Numerical solution of the equlibrium integral equation
    
## Compiling flags:
    -DSHOUT - turn off full info and turn on progress bar
    -DBAR_WIDTH=50 - width in characters of progress bar (default value is 70)
    
## The list of values you must give the program in command line
    - The number of iterations
    - The number of nodes in grid
    - The limit of integration
    - The path to store data file
    So executing of this programm looks like:
        ./neuman 1000 5001 20.0 calculated_data.plt
        
    Nota bene: you can miss last values. So, programm will use default values

## Default values:
    5001 - nodes
    1000 - iterations
    20.0 - integration limit
    "graph.plt" - file with stored data (x, accurate solution, calculated solution, relative error)
    
**Author isn't responsible for your psycology helth that can be damaged during reading this code**
