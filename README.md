# Neuman
Numerical solution of the equlibrium integral equation

## Default values:
    5001 - nodes
    1000 - iterations
    "graph.plt" - file with stored data (x, accurate solution, calculated solution, relative error)
    70 - width of the progress bar in characters
 For each iteration the programm shows a C1-norm of difference between the last and current solution iteration
    
## The list of flags, that you must add during compilation to change parameters of program:
    -DSHOUT - off full info about each iteration and on the progress bar
    -DPATH=\"new_path.plt\" - change file for storing data
    -DBAR_WIDTH=80 - change the progress bar width
    
## The list of values you must give the program in command line
    - The number of iterations
    - The number of nodes in greed
    So executing of this programm looks like:
        * ./neuman 1000 5001 *
    
**Author isn't responsible for your psycology helth that can be damaged during reading this code**
