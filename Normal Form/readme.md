# Normal Form Code in Julia
All normal form implementation is done in the Julia programming language for performance reasons. To run the Julia code, follow the steps below.

1. Install julia from https://julialang.org/downloads/
2. Open a terminal window at `NFQPOs`.
3. In the terminal, run 
    ``` zsh
    julia --project
    ] activate .
    ] instantiate
    ``` 
> 📘  Note
>
> If working in the terminal is unappealing, then Julia can also be used through VSCode. In this case, just make sure the REPL is at the right path, which can be ensured by using `] activate` when in the `NFQPOs` folder

4. The data files required to actually execute the normal form analysis are included. If lost, simply run `./data/create_required_data` to regenerate them. Alternatively, each data file can be manually recreated using `./data/create_store`

5. The rest of the code should be fairly easy to follow. The utilites used to compute and operate with the normal form can be found in `Subroutines/NFCR3BP.jl`