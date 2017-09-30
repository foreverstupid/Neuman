# Neuman
Numerical solution of integral equation

## Default values:
    5001 - nodes
    1000 - iterations
    "graph.plt" - file with stored data (x, accurate solution, calculated solution, relative error)
    70 - width of progress bar in characters
 For each iteration programm shows C1-norm of difference between last and current solution iteration
    
## The list of flags, that you must add during compilation to change parameters of program:
    -DSHOUT - off full info about each iteration and on progress bar
    -DITERS=20001 - change iteration count
    -DN_COUNT=42 - change nodes count
    -DPATH=\"new_path.plt\" - change file for storing data
    -DBAR_WIDTH=80 - change progress bar width
    
**Author isn't responsible for your psycology helth that can be damaged during reading this code**
